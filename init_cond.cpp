#include "maxol.h"

template <int c, int N0, int N1, int N2>
static void __set_init_cond(double *A)
{
	for (int i = 0; i < N0; i++)
		for (int k = 0; k < N2; k++)
			for (int j = 0; j < N1; j++) {
					int l = i + N0*(j + N1*k);
					// the smallest positive value and the biggest negative value.
					A[l] = l%2 ? 0x0000000000000001 : 0x8000000000000001;
				}
}

void set_init_cond(double *Ep, double *Eq, double *Er,
                   double *BP, double *BQ, double *BR)
{
	__set_init_cond<0, NP-1, NQ, NR>(Ep);
	__set_init_cond<1, NQ-1, NR, NP>(Eq);
	__set_init_cond<2, NR-1, NP, NQ>(Er);
	__set_init_cond<3, NP, NQ-1, NR-1>(BP);
	__set_init_cond<4, NQ, NR-1, NP-1>(BQ);
	__set_init_cond<5, NR, NP-1, NQ-1>(BR);
}
