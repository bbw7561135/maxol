#include <assert.h>
#include <math.h>

#include "coordinate.h"
#include "../config/comp_param.h"
#include "../config/phys_param.h"

extern double dt;

/*
 * Discretize the integral form of the Faraday's equation which is like below:
 *   _           __  ->
 *  / ->  ->    //  @E   ->
 *  O E *dr = - || ----*dS
 * _/A         _//A @t
 *
 * I integrate magnetic flux density in a surface element like ###;
 * and countour-integrate electric fields along the edge of it.
 * Magnetic flux density is represented by the value at the center of
 * surface element.
 * Electric fields are represented by the center of each edge.
 *
 *            q=j-1/2   q=j+1/2
 *              /         /
 *      -------+---------+--- r=k-1/2
 *            / # # # # /
 *           / # # # # /
 *          / # # # # /
 *      ---+---------+------- r=k+1/2
 *        /         /
 *
 * Vewing from r-axis, the surface element is like ===. The allow vertical to
 * the surface element is the magnetic flux density vector to be evoluted. The
 * * points and the --> allow mean electric field vectors to be countour-
 * integrated.
 *
 *           q=j-1        q=j        q=j+1
 *            /           /           /
 *           +-----------+-----------+-- p=i-1/2
 *          /         ^ /           /
 *         /          |/           /
 *        /     *===========*     /
 *       /           /-->        /
 *      /           /           /
 *     +-----------+-----------+-- p=i+1/2
 *
 * ----------------------------------------------------------------------------
 *
 * As evoluting B(i,j',k'), definition of variable is below.
 *
 *                       lenA
 *                  <----------->
 *                        E1(i,j',k'-1/2)
 *                        +---->
 *            ^     +-----------+     ^
 *            |     |           |     |
 *        lenC|    +|    dS     |+    |lenD
 *            |    ||           ||    |
 *            v    |+-----------+|    v
 * E2(i,j'-1/2,k') v      +----> v E2(i,j'+1/2,k')
 *                        E1(i,j',k'+1/2)
 *                  <----------->
 *                       lenB
 *
 * Range of (i,j',k') is
 *
 *   i :  1  ~  NP-2
 *   k': 1/2 ~ NR-3/2
 *   j': 1/2 ~ NQ-3/2
 *
 *                     /
 *                    /| p1 = j'
 *                   / |
 *                  /  |
 *                 /   |
 *                /    |
 * p2 = k'        |  ---->
 *        /-------|----+------------/
 *       /        |   /|           /
 *      /     .   |  /.|   .      /___ p0 = i
 *     /      |   | // |   |     /
 *    /       v   |/L  |   v    /
 *   /------------+------------/
 *                |  ---->
 *                |    /
 *                |   /
 *                |  /
 *                | /
 *                |/
 *
 * Calculation is this order. Thus I reused the values on the right edge.
 *
 *   p1 1   2   3
 * p2 +---+---+---+
 *  1 | 1-------> |
 *    +---+---+---+
 *  3 | 2-------> |
 *    +---+---+---+
 *  4 | 3-------> |
 *    +---+---+---+
 */

template <int c, int N0, int N1, int N2>
static double __evolute_magnt(double *B, const double *E1, const double *E2)
{
	// total energy about a component of magnetic field
	double eng = 0.0;

	// Extract for-loop for vectorization.
	// Don't extract inner loop because reuse a value: intD.
	for (int ik = 0; ik < (N0-2)*(N2-1); ik++) {
		const int i = ik / (N2-1) + 1; // for i: 1 ~ N0-2
		const int k = ik % (N2-1);     // for k: 0 ~ N2-2, k': 0.5 ~ N2-1.5

		// (p0, p1, p2) is the center of the surface element.
		const double p0 = (double)i;
		const double p2 = (double)k + 0.5;

		// Calculate intD previously because cannot reuse it at the first loop.

		// length of the edge (i,1,k')
		const double lenD = vector_length(
				covariant_basic_vector<2, c>({p0, 1.0, p2}));
		const int lD = i + N0*(1 + N1*k);
		double intD = E2[lD]*lenD;

		// for j': 0.5 ~ N1-1.5
		for (int j = 0; j < N1-1; j++) {
			const double p1 = (double)j + 0.5;

			/*
			 * The volume equals to Jacobian.
			 *   +-----+    ^
			 *  / \   / \   | 1.0 / (length of contravariant basic vector)
			 * +-----+   \  |
			 *  \   +-\---+ v
			 *   \ / dS\ /
			 *    +-----+
			 */
			const double jcb = jacobian<c>({p0, p1, p2});
			const double dS = jcb * vector_length(
					contravariant_basic_vector<c, c>({p0, p1, p2}));

			// length of the edge (i,j',k'-1/2)
			const double lenA = vector_length(
					covariant_basic_vector<1, c>({p0, p1, p2-0.5}));
			// length of the edge (i,j',k'+1/2)
			const double lenB = vector_length(
					covariant_basic_vector<1, c>({p0, p1, p2+0.5}));
			// length of the edge (i,j'+1/2,k')
			const double lenD = vector_length(
					covariant_basic_vector<2, c>({p0, p1+0.5, p2}));

			// position of E1(i,j',k'-1/2)
			const int lA = k-1 + N2*(i+ N0*j);
			// position of E1(i,j',k'+1/2)
			const int lB = lA + 1;
			// position of E2(i,j'+1/2,k')
			const int lD = i + N0*(j + N1*k);

			// integral of electric field along each edge
			const double intA = E1[lA]*lenA;
			const double intB = E1[lB]*lenB;
			const double intC = intD; // reuse from previous loop
			intD = E2[lD]*lenD;

			// contour integral of electric field around surface element
			const double oint = - intA + intB - intC +  intD;

			// position of B(i,j',k')
			const int l = j + (N1-1)*(k + (N2-1)*i);
			// Faraday's equation
			const double B_new = B[l] - dt * oint / dS;
			B[l] = B_new;

			// Sum up energy in the grid
			eng += 0.5 / magnetic_permeability * B_new * B_new * jcb;
		}
	}

	return eng;
}

/*
 * Evolute electric fields with leap-frog method following Faraday's law.
 */
double evolute_magnt(
		// magnetic flux densities at nt
		double *BP, double *BQ, double *BR,
		// electric fields at nt + 0.5
		const double *Ep, const double *Eq, const double *Er)
{
	double eng_P = __evolute_magnt<0, NP, NQ, NR>(BP, Eq, Er);
	double eng_Q = __evolute_magnt<1, NQ, NR, NP>(BQ, Er, Ep);
	double eng_R = __evolute_magnt<2, NR, NP, NQ>(BR, Ep, Eq);
	return (eng_P + eng_Q + eng_R)/3.0;
}
