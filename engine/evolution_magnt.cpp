#include <assert.h>
#include <math.h>

#include "coordinate.h"
#include "../config/comp_param.h"
#include "../config/phys_param.h"

extern double dt;

/*
 * I integrate magnetic flux density in a surface element like the square
 * EFGH; and countour-integrate electric field along the edge of it.
 *
 * The surface element is not along the computational grid because it is along
 * contravariance basic vectors of the coodinate system.
 * It is vertical to the grid. In the figure, AB, DC are vertical with EH, FG;
 * and AD, BC are vertical with EF, HG. Additionaly, EFGH is vertical with the
 * r-axis.
 *
 * The edges of surface element is on the points: I, J, K, L, where electric
 * field vectors are put.
 *
 *              p-1/2
 *     H |._     /                      p+1/2
 *       |  `-._/                        /
 *       |     /`-._                    /
 *       |  D /     `-._  K            /
 *       |   +----------`+._----------+------- q-1/2
 *       |  /               `-._  G  / C
 *       | /                    `-. /
 *       |/                       |/
 *     L +                        + J
 *      /|                       /|
 *     / `-._                   / |
 *    /  E   `-._              /  |
 * A +-----------`+._---------+---|--- q+1/2
 *               I   `-._      B  |
 *                       `-._     |
 *                           `-._ |
 *                               `| F
 *
 * As evoluting BP(i,j',k'), definition of variable is below.
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

template <int c, int N0, int N1, int N2>
static void __evolute_magnt(double *B, const double *E1, const double *E2)
{
#ifdef _OPENMP
#pragma omp for
#endif
	// I extract for-loop for vectorization.
	// I don't extract inner loop because I reuse a value: intD.
	// for i: 1 ~ N0-2, for k: 0 ~ N2-2
	for (int ik = 0; ik < (N2-1)*(N0-2); ik++) {
		const int i = ik / (N2-1) + 1;
		const int k = ik % (N2-1);

		// (p0, p1, p2) is the center of the surface element.
		const double p0 = (double)i;
		const double p2 = (double)k + 0.5;

		// Calculate intD previously
		// because we cannot reuse it at the first loop.
		double intD;
		{
			// length of the edge (i',1,k')
			const double lenD =
					len_of_normalized_contravariant_basic_vector<2, c>(
							p0, 1.5, p2);
			const int lD = i + N0*(1 + N1*k);
			intD = E2[lD]*lenD;
		}

		for (int j = 0; j < N1-1; j++) {
			const double p1 = (double)j + 0.5;

			/*
			 * The volume equals to Jacobian.
			 *   +-----+    ^
			 *  / \   / \   | 1 / (length of covariant basic vector)
			 * +-----+   \  |
			 *  \   +-\---+ v
			 *   \ / dS\ /
			 *    +-----+
			 */
			const double dS = jacobian<c>(p0, p1, p2) *
					len_of_covariant_basic_vector<c, c>(p0, p1, p2);
			assert(dS != 0.0);

			// length of the edge (i,j',k'-1/2)
			const double lenA =
					len_of_normalized_contravariant_basic_vector<1, c>(
							p0, p1, p2-0.5);
			// length of the edge (i,j',k'+1/2)
			const double lenB =
					len_of_normalized_contravariant_basic_vector<1, c>(
							p0, p1, p2+0.5);
			// length of the edge (i,j'+1/2,k')
			const double lenD =
					len_of_normalized_contravariant_basic_vector<2, c>(
							p0, p1+0.5, p2);

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
			const int l = j + N1*(k + N2*i);
			// Maxwell–Faraday equation
			B[l] -= dt*oint/dS;
		}
	}
}

/*
 * Evolute electric fields with leap-frog method following Faraday's law.
 */
void evolute_magnt(
		// magnetic flux densities at nt
		double *BP, double *BQ, double *BR,
		// electric fields at nt + 0.5
		const double *Ep, const double *Eq, const double *Er)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
	{
		__evolute_magnt<0, NP, NQ, NR>(BP, Eq, Er);
		__evolute_magnt<1, NQ, NR, NP>(BQ, Er, Ep);
		__evolute_magnt<2, NR, NP, NQ>(BR, Ep, Eq);
	}
}
