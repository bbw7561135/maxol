#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "coordinate.h"
#include "variable.h"
#include "../config/comp_param.h"
#include "../config/phys_param.h"

extern double dt;

/*
 * Discretize the integral form of the Ampere's equation which is like below:
 *   _             __                ->
 *  / ->  ->      //  ->            @E    ->
 *  O B *dr = \mu || (j + \epsilon ----)*dS
 * _/A           _//A               @t
 *
 * I integrate electric field in a surface element like the square
 * EFGH; and countour-integrate magnetic flux densities along the edge of it.
 *
 * The surface element is not along computational grid, but it is along
 * contravariant basic vectors of the coordinate system.
 * It is vertical to the grid. In the figure, AB, DC are vertical with EH, FG;
 * and AD, BC are vertical with EF, HG. Additionally, EFGH is vertical with the
 * p-axis.
 *
 * The edges of surface element is on the points: I, J, K, L, where electric
 * field vectors are put.
 *
 *             q=j-1/2
 *     H |._     /                     q=j+1/2
 *       |  `-._/                        /
 *       |     /`-._                    /
 *       |  D /     `-._  K            /
 *       |   +----------`+._----------+------- r=k-1/2
 *       |  /               `-._  G  / C
 *       | /                    `-. /
 *       |/                       |/
 *     L +                        + J
 *      /|                       /|
 *     / `-._                   / |
 *    /  E   `-._              /  |
 * A +-----------`+._---------+---|--- r=k+1/2
 *               I   `-._      B  |
 *                       `-._     |
 *                           `-._ |
 *                               `| F
 *
 * Viewing from r-axis, the surface element is like [. The --> allow is the
 * magnetic flux density vector to be evoluted. The * points and the --> allow
 * mean electric field vectors to be countour-integrated.
 *
 * On the * points, electric fields vector is not defined. I substitute
 * electric fields on + points.
 *
 *                 p=i
 *                  /
 *      -------*---+--- q=j-1/2
 *             [  /
 *            ^[ /
 *            |[/
 *             [-->
 *            /[
 *           / [
 *          /  [
 *      ---+---*------- q=j+1/2
 *        /
 *
 * ----------------------------------------------------------------------------
 *
 * As evoluting E(i',j,k), definition of variable is below.
 *
 *                      lenA
 *                 <----------->
 *                       B1(i',j,k-1/2)
 *                       +---->
 *           ^     +-----------+     ^
 *           |     |           |     |
 *       lenC|    +|    dS     |+    |lenD
 *           |    ||           ||    |
 *           v    |+-----------+|    v
 * B2(i',j-1/2,k) v      +----> v B2(i',j+1/2,k)
 *                       B1(i',j,k+1/2)
 *                 <----------->
 *                      lenB
 *
 * Range of (i',j,k) is
 *
 *   i': 1/2 ~ N0-3/2
 *   j :  1  ~  N1-2
 *   k :  1  ~  N2-2
 *
 *         p1=j
 * p0=i'     |
 *     |     |    |
 *  ---+-----==>--+---
 *     |          |
 *     |          |  ____ p2=k
 *    ||    /     ||
 *    v|   L      |v
 *  ---+-----==>--+---
 *     |          |
 *
 * Calculation is this order. Thus I reused the values on the right edge.
 *
 *   j 1   2   3
 * k +---+---+---+
 * 1 | 1-------> |
 *   +---+---+---+
 * 3 | 2-------> |
 *   +---+---+---+
 * 4 | 3-------> |
 *   +---+---+---+
 */

/*
 * contravariant basic vectors normalized to the same size (Jacobian) as
 * covariant basic vectors.
 */
template <int c, int d>
static inline struct vector normalized_contravariant_basic_vector(struct vector position)
{
	const double coefficient = pow(jacobian<d>(position), 2.0/3.0);
	return vector_multiple(
				contravariant_basic_vector<c, d>(position), coefficient);
}

template <int c, int N0, int N1, int N2>
static double __evolute_elect(double *E, const double *B1, const double *B2)
{
	// total energy about a component of electric field
	double eng = 0.0;

	// Extract for-loop for vectorization.
	// Don't extract inner loop because reuse a value: intD.
	// for i: 0 ~ N0-2, for k: 1 ~ N2-2
	for (int ik = 0; ik < (N2-2)*(N0-1); ik++) {
		const int i = ik / (N2-2);
		const int k = ik % (N2-2) + 1;

		// (p0, p1, p2) is the center of the surface element.
		// p1 is declared after.
		const double p0 = (double)i + 0.5;
		const double p2 = (double)k;

		// Calculate intD previously because cannot reuse it at the first loop.

		// length of the edge (i',3/2,k)
		const double lenD = vector_length(
								covariant_basic_vector<2, c>({p0, 1.5, p2}));
		double intD = lenD * (*magnt_flux<2, N0, N1, N2>(i, 1, k, B2));

		for (int j = 1; j < N1-1; j++) {
			const double p1 = (double)j;

			/*
			 * The volume equals to Jacobian.
			 *   +-----+    ^
			 *  / \   / \   | length of covariant basic vector
			 * +-----+   \  |
			 *  \   +-\---+ v
			 *   \ / dS\ /
			 *    +-----+
			 */
			const double jcb = jacobian<c>({p0, p1, p2});
			const double dS = jcb / vector_length(
					covariant_basic_vector<c, c>({p0, p1, p2}));

			// length of the edge (i',j,k-1/2)
			const double lenA = vector_length(
				normalized_contravariant_basic_vector<1, c>({p0, p1, p2-0.5}));
			// length of the edge (i',j,k+1/2)
			const double lenB = vector_length(
				normalized_contravariant_basic_vector<1, c>({p0, p1, p2+0.5}));
			// length of the edge (i',j+1/2,k)
			const double lenD = vector_length(
				normalized_contravariant_basic_vector<2, c>({p0, p1+0.5, p2}));

			// integral of magnetic flux density along each edge
			const double intA = lenA *
					(*magnt_flux<1, N0, N1, N2>(i, j, k-1, B1));
			const double intB = lenB *
					(*magnt_flux<1, N0, N1, N2>(i, j, k  , B1));
			const double intC = intD; // reuse from previous loop
			intD = lenD * (*magnt_flux<2, N0, N1, N2>(i, j, k, B2));

			// contour integral of magnetic flux density around surface element
			const double oint = - intA + intB - intC + intD;

			// E at (i',j,k)
			double *e = elect_field<N1, N2>(i, j, k, E);
			// Ampere's circuital law. Omitting electric current term.
			const double E_new = *e +
					1.0 / (magnetic_permeability * electric_permittivity) *
					dt * oint / dS;
			*e = E_new;

			// Sum up energy in the grid
			eng += 0.5 * electric_permittivity * E_new * E_new * jcb;
		}
	}

	return eng;
}

/*
 * Evolute electric fields with leap-frog method
 * following Ampere's circuital law.
 * NOTICE:  Now I omit electric current term.
 */
double evolute_elect(
		// electric fields at nt
		double *Ep, double *Eq, double *Er,
		 // magnetic flux densities at nt + 0.5
		const double *BP, const double *BQ, const double *BR)
{
	double eng_p = __evolute_elect<0, NP, NQ, NR>(Ep, BQ, BR);
	double eng_q = __evolute_elect<1, NQ, NR, NP>(Eq, BR, BP);
	double eng_r = __evolute_elect<2, NR, NP, NQ>(Er, BP, BQ);
	return (eng_p + eng_q + eng_r)/3.0;
}
