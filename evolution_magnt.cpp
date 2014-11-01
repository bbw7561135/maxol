#include <assert.h>
#include <math.h>

#include "comp_param.h"
#include "coordinate.h"
#include "phys_param.h"

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
 *                       dtA
 *                  <----------->
 *                        E1(i,j',k'-1/2)
 *                        +---->
 *            ^     +-----------+     ^
 *            |     |           |     |
 *         dtC|    +|    dS     |+    |dtD
 *            |    ||           ||    |
 *            v    |+-----------+|    v
 * E2(i,j'-1/2,k') v      +----> v E2(i,j'+1/2,k')
 *                        E1(i,j',k'+1/2)
 *                  <----------->
 *                       dtB
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
static void __evolute_magnetic(double *B, const double *E1, const double *E2)
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
			const double dtD =
					length_of_normalized_contravariant_basic_vector<2, c>(
							p0, 1.5, p2);
			const int lD = i + N0*(1 + N1*k);
			intD = E2[lD]*dtD;
		}

		for (int j = 0; j < N1-1; j++) {
			const double p1 = (double)j + 0.5;

			// components of contravariant basic vectors along surface element
			const double dx_dP1 =
				normalized_contravariant_basic_vector<0, 1, c>(p0, p1, p2);
			const double dy_dP1 =
				normalized_contravariant_basic_vector<1, 1, c>(p0, p1, p2);
			const double dz_dP1 =
				normalized_contravariant_basic_vector<2, 1, c>(p0, p1, p2);

			const double dx_dP2 =
				normalized_contravariant_basic_vector<0, 2, c>(p0, p1, p2);
			const double dy_dP2 =
				normalized_contravariant_basic_vector<1, 2, c>(p0, p1, p2);
			const double dz_dP2 =
				normalized_contravariant_basic_vector<2, 2, c>(p0, p1, p2);

			// cross product vector
			const double n0 = dy_dP1*dz_dP2 - dy_dP2*dz_dP1;
			const double n1 = dz_dP1*dx_dP2 - dz_dP2*dx_dP1;
			const double n2 = dx_dP1*dy_dP2 - dx_dP2*dy_dP1;

			// area of the surface element
			const double dS = sqrt(n0*n0 + n1*n1 + n2*n2);
			assert(dS != 0.0);

			// length of the edge (i,j',k'-1/2)
			const double dtA =
					length_of_normalized_contravariant_basic_vector<1, c>(
							p0, p1, p2-0.5);
			// length of the edge (i,j',k'+1/2)
			const double dtB =
					length_of_normalized_contravariant_basic_vector<1, c>(
							p0, p1, p2+0.5);
			// length of the edge (i,j'+1/2,k')
			const double dtD =
					length_of_normalized_contravariant_basic_vector<2, c>(
							p0, p1+0.5, p2);

			// position of E1(i,j',k'-1/2)
			const int lA = k-1 + N2*(i+ N0*j);
			// position of E1(i,j',k'+1/2)
			const int lB = lA + 1;
			// position of E2(i,j'+1/2,k')
			const int lD = i + N0*(j + N1*k);

			// integral of electric field along each edge
			const double intA = E1[lA]*dtA;
			const double intB = E1[lB]*dtB;
			const double intC = intD; // reuse from previous loop
			intD = E2[lD]*dtD;

			// contour integral of electric field around surface element
			const double oint = - intA + intB - intC + intD;

			// position of B(i,j',k')
			const int l = j + N1*(k + N2*i);
			// Maxwell–Faraday equation
			B[l] -= DT*oint/dS;
		}
	}
}

/*
 * Evolute electric fields with leap-frog method following Faraday's law.
 */
void evolute_magnetic(
		// magnetic flux densities at nt
		double *BP, double *BQ, double *BR,
		// electric fields at nt + 0.5
		const double *Ep, const double *Eq, const double *Er)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
	{
		__evolute_magnetic<0, NP, NQ, NR>(BP, Eq, Er);
		__evolute_magnetic<1, NQ, NR, NP>(BQ, Er, Ep);
		__evolute_magnetic<2, NR, NP, NQ>(BR, Ep, Eq);
	}
}
