#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "coordinate.h"
#include "../config/comp_param.h"
#include "../config/phys_param.h"

// TODO: To fix some bug (buffer overrun?)

/*
 * I integrate electric field in a surface element like ###;
 * and countour-integrate magnetic flux density along the edge of it.
 * Electric field is represented by the value at the center of
 * surface element.
 * Magnetic flux densities are represented by the center of each edge.
 *
 *             p-1/2     p+1/2
 *              /         /
 *      -------+---------+--- q-1/2
 *            / # # # # /
 *           / # # # # /
 *          / # # # # /
 *      ---+---------+------- q+1/2
 *        /         /
 *
 * As evoluting E(i',j,k), definition of variable is below.
 *
 *                      dtA
 *                 <----------->
 *                       B1(i',j,k-1/2)
 *                       +---->
 *           ^     +-----------+     ^
 *           |     |           |     |
 *        dtC|    +|    dS     |+    |dtD
 *           |    ||           ||    |
 *           v    |+-----------+|    v
 * B2(i',j-1/2,k) v      +----> v B2(i',j+1/2,k)
 *                       B1(i',j,k+1/2)
 *                 <----------->
 *                      dtB
 *
 * Range of (i',j,k) is
 *
 *   i': 1/2 ~ N0-3/2
 *   j :  1  ~  N1-2
 *   k :  1  ~  N2-2
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
static void __evolute_elect(double *E, const double *B1, const double *B2)
{
#ifdef _OPENMP
#pragma omp for
#endif
	// Extract for-loop for vectorization.
	// Don't extract inner loop because a value intD is reused.
	// for i: 0 ~ N0-2, for k: 1 ~ N2-2
	for (int ik = 0; ik < (N2-2)*(N0-1); ik++) {
		const int i = ik / (N2-2);
		const int k = ik % (N2-2) + 1;

		// (p0, p1, p2) is the center of the surface element.
		// p1 is declared after.
		const double p0 = (double)i + 0.5;
		const double p2 = (double)k;

		// Calculate intD previously because we cannot reuse it at the first loop.

		// length of the edge (i',3/2,k)
		const double dtD = len_of_covariant_basic_vector<2, c>(p0, 1.5, p2);
		const int lD = i + (N0-1)*(1 + (N1-1)*k);
		double intD = B2[lD]*dtD;

		for (int j = 1; j < N1-1; j++) {
			const double p1 = (double)j;

#if 0 /* Calculation with jacobian is more simple */
			// components of covariant basic vectors along surface element
			const double dx_dp1 = covariant_basic_vector<0, 1, c>(p0, p1, p2);
			const double dy_dp1 = covariant_basic_vector<1, 1, c>(p0, p1, p2);
			const double dz_dp1 = covariant_basic_vector<2, 1, c>(p0, p1, p2);

			const double dx_dp2 = covariant_basic_vector<0, 2, c>(p0, p1, p2);
			const double dy_dp2 = covariant_basic_vector<1, 2, c>(p0, p1, p2);
			const double dz_dp2 = covariant_basic_vector<2, 2, c>(p0, p1, p2);

			// cross product vector
			const double n0 = dy_dp1 * dz_dp2 - dy_dp2 * dz_dp1;
			const double n1 = dz_dp1 * dx_dp2 - dz_dp2 * dx_dp1;
			const double n2 = dx_dp1 * dy_dp2 - dx_dp2 * dy_dp1;

			// area of the surface element
			const double dS = sqrt(n0*n0 + n1*n1 + n2*n2);
#endif
			const double dS = jacobian<c>(p0, p1, p2) *
						len_of_contravariant_basic_vector<c, c>(p0, p1, p2);
			assert(dS != 0.0);

			// length of the edge (i',j,k-1/2)
			const double dtA =
					len_of_covariant_basic_vector<1, c>(p0, p1, p2-0.5);
			// length of the edge (i',j,k+1/2)
			const double dtB =
					len_of_covariant_basic_vector<1, c>(p0, p1, p2+0.5);
			// length of the edge (i',j+1/2,k)
			const double dtD =
					len_of_covariant_basic_vector<2, c>(p0, p1+0.5, p2);

			// position of B1(i',j,k-1/2)
			const int lA = k-1 + (N2-1)*(i + (N0-1)*j);
			// position of B1(i',j,k+1/2)
			const int lB = lA + 1;
			// position of B2(i',j+1/2,k)
			const int lD = i + (N0-1)*(j + (N1-1)*k);

			// integral of magnetic flux density along each edge
			const double intA = B1[lA]*dtA;
			const double intB = B1[lB]*dtB;
			const double intC = intD; // reuse from previous loop
			intD = B2[lD]*dtD;

			// contour integral of magnetic flux density around surface element
			const double oint = - intA + intB - intC + intD;

			// position of E(i',j,k)
			const int l = j + N1*(k + N2*i);
			// Ampère's circuital law. Omitting electric current term.
			E[l] += DT/(permeability*permittivity)*oint/dS;
		}
	}
}

/*
 * Evolute electric fields with leap-frog method
 * following Ampère's circuital law.
 * NOTICE:  Now I omit electric current term.
 */
void evolute_elect(
		// electric fields at nt
		double *Ep, double *Eq, double *Er,
		 // magnetic flux densities at nt + 0.5
		const double *BP, const double *BQ, const double *BR)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
	{
		__evolute_elect<0, NP, NQ, NR>(Ep, BQ, BR);
		__evolute_elect<1, NQ, NR, NP>(Eq, BR, BP);
		__evolute_elect<2, NR, NP, NQ>(Er, BP, BQ);
	}
}
