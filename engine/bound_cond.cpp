#include "../config/comp_param.h"

/*
 *                / 6
 *           3   L
 *        +--|---- +
 *       /   v    /|
 *      /        / |
 * 5   +--------+ <--2    ^ k
 * --> |        |  +      |
 *     |   1    | /       +---> j
 *     |        |/       /
 *     +--------+       L i
 *           ^
 *           | 4
 *
 * View from direction 1.
 */

template <int c, int N0, int N1, int N2>
static void __bound_cond_elect(double *E, double t)
{
	// tangential components to boundary is zero

#ifdef _OPENMP
#pragma omp for
#endif
	for (int i = 0; i < N0; i++){
		for (int j = 1; j < N1-1; j++) {
			// Bound 4
			const int k0 = 0;
			const int l0 = j + N1*(k0 + N2*i);
			E[l0] = 0.0;

			// Bound 3
			const int k1 = N2-1;
			const int l1 = j + N1*(k1 + N2*i);
			E[l1] = 0.0;
		}
		// boundary j=0 and j=N1-1
		for (int k = 0; k < N2; k++) {
			// Bound 5
			const int j0 = 0;
			const int l0 = j0 + N1*(k + N2*i);
			E[l0] = 0.0;

			// Bound 2
			const int j1 = N1-1;
			const int l1 = j1 + N1*(k + N2*i);
			E[l1] = 0.0;
		}
	}
}

void bound_cond_elect(double *Ep, double *Eq, double *Er, double t)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
	{
		__bound_cond_elect<0, NP-1, NQ, NR>(Ep, t);
		__bound_cond_elect<1, NQ-1, NR, NP>(Eq, t);
		__bound_cond_elect<2, NR-1, NP, NQ>(Er, t);
	}
}

// Only vertical components face on the boundaries. However, they are not
// affect because light is transverse wave. Certainly, they not used for
// computing. I don't think this condition might be needed. I don't know why
// boundary condition in magnetic field is not required.

template <int c, int N0, int N1, int N2>
static void __bound_cond_magnt(double *B, double t)
{
	// vertical component to boundary is zero
#ifdef _OPENMP
#pragma omp for
#endif
	for (int k = 0; k < N2; k++)
		for (int j = 0; j < N1; j++) {
			// Bound 6
			const int i0 = 0;
			const int l0 = j + N1*(k + N2*i0);
			B[l0] = 0.0;

			// Bound 1
			const int i1 = N0-1;
			const int l1 = j + N1*(k + N2*i1);
			B[l1] = 0.0;
		}
}

void bound_cond_magnt(double *BP, double *BQ, double *BR, double t)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
	{
		__bound_cond_magnt<0, NP, NQ-1, NR-1>(BP, t);
		__bound_cond_magnt<1, NQ, NR-1, NP-1>(BQ, t);
		__bound_cond_magnt<2, NR, NP-1, NQ-1>(BR, t);
	}
}
