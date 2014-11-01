#include "../config/comp_param.h"

template <int c, int N0, int N1, int N2>
static void __bound_cond_elect(double *E, double t)
{
	// tangential components to boundary is zero

#ifdef _OPENMP
#pragma omp for
#endif
	for (int i = 0; i < N0-1; i++){
		// boundary k=0 and k=N2-1
		for (int j = 1; j < N1-1; j++) {
			const int k0 = 0;
			const int l0 = j + N1*(k0 + N2*i);
			E[l0] = 0.0;

			const int k1 = N2-1;
			const int l1 = j + N1*(k1 + N2*i);
			E[l1] = 0.0;
		}
		// boundary j=0 and j=N1-1
		for (int k = 1; k < N2-1; k++) {
			const int j0 = 0;
			const int l0 = j0 + N1*(k + N2*i);
			E[l0] = 0.0;

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
		__bound_cond_elect<0, NP, NQ, NR>(Ep, t);
		__bound_cond_elect<1, NQ, NR, NP>(Eq, t);
		__bound_cond_elect<2, NR, NP, NQ>(Er, t);
	}
}

// Only vertical components face on the boundaries.
// However, they are not affect because light is transverse wave.
// Certainly, they not used for computing.
// I think they might not be needed.
// I don't know why boundary condition in magnetic field is not required.

template <int c, int N0, int N1, int N2>
static void __bound_cond_magnt(double *B, double t)
{
	// vertical component to boundary is zero
#ifdef _OPENMP
#pragma omp for
#endif
	for (int k = 0; k < N2-1; k++)
		for (int j = 0; j < N1-1; j++) {
			const int i0 = 0;
			const int l0 = j + (N1-1)*(k + (N2-1)*i0);
			B[l0] = 0.0;

			const int i1 = N0-1;
			const int l1 = j + (N1-1)*(k + (N2-1)*i1);
			B[l1] = 0.0;
		}
}

void bound_cond_magnt(double *BP, double *BQ, double *BR, double t)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
	{
		__bound_cond_magnt<0, NP, NQ, NR>(BP, t);
		__bound_cond_magnt<1, NQ, NR, NP>(BQ, t);
		__bound_cond_magnt<2, NR, NP, NQ>(BR, t);
	}
}
