#include "variable.h"
#include "vector.h"
#include "../config/comp_param.h"
#include "../config/my_bound_cond.h"

/*
 *                / A
 *           F   L
 *        +--|---- +
 *       /   v    /|
 *      /        / |
 * B   +--------+ <--E    ^ k
 * --> |        |  +      |
 *     |   D    | /       +---> j
 *     |        |/       /
 *     +--------+       L i
 *           ^
 *           | C
 *
 */

#define BOUND_A 0
#define BOUND_B 1
#define BOUND_C 2
#define BOUND_D 3
#define BOUND_E 4
#define BOUND_F 5

template <int c, int b>
static inline double __my_bound_cond_elect(int p0, int p1, int p2, double *E)
{
	int p, q, r;
	switch (c) {
	case 0: p = p0, q = p1, r = p2; break;
	case 1: p = p2, q = p0, r = p1; break;
	case 2: p = p1, q = p2, r = p0; break;
	}
	// case c == 0: bound(A,B,C,D,E,F) -> bound(0,1,2,3,4,5)
	// case c == 1: bound(A,B,C,D,E,F) -> bound(1,2,0,4,5,3)
	// case c == 2: bound(A,B,C,D,E,F) -> bound(2,0,1,5,3,4)
	const int bound = b < 3 ? (b + c) % 3 : (b + c) % 3 + 3;
	return my_bound_cond_elect<c, bound>(p, q, r, E);
}

template <int c, int N0, int N1, int N2>
static void __bound_cond_elect(double *E)
{
	for (int i = 0; i < N0-1; i++) {
		for (int j = 1; j < N1-1; j++) {
			// Bound C
			const int k0 = 0;
			*elect_field<N1, N2>(i, j, k0, E) =
					__my_bound_cond_elect<c, BOUND_C>(i, j, k0, E);

			// Bound F
			const int k1 = N2-1;
			*elect_field<N1, N2>(i, j, k1, E) =
					__my_bound_cond_elect<c, BOUND_F>(i, j, k1, E);
		}
		for (int k = 0; k < N2; k++) {
			// Bound B
			const int j0 = 0;
			*elect_field<N1, N2>(i, j0, k, E) =
					__my_bound_cond_elect<c, BOUND_B>(i, j0, k, E);

			// Bound E
			const int j1 = N1-1;
			*elect_field<N1, N2>(i, j1, k, E) =
					__my_bound_cond_elect<c, BOUND_E>(i, j1, k, E);
		}
	}
}

void bound_cond_elect(double *Ep, double *Eq, double *Er)
{
	__bound_cond_elect<0, NP, NQ, NR>(Ep);
	__bound_cond_elect<1, NQ, NR, NP>(Eq);
	__bound_cond_elect<2, NR, NP, NQ>(Er);
}

// Only vertical components face on the boundaries. However, they are not
// affect because light is transverse wave. Certainly, they not used for
// computing. I don't think this condition might be needed. I don't know why
// boundary condition in magnetic field is not required.

template <int c, int N0, int N1, int N2>
static void __bound_cond_magnt(double *B)
{
	// vertical component to boundary is zero
	for (int k = 0; k < N2-1; k++)
		for (int j = 0; j < N1-1; j++) {
			// Bound A
			const int i0 = 0;
			*magnt_flux<N1, N2>(i0, j, k, B) = 0.0;

			// Bound D
			const int i1 = N0-1;
			*magnt_flux<N1, N2>(i1, j, k, B) = 0.0;
		}
}

void bound_cond_magnt(double *BP, double *BQ, double *BR)
{
		__bound_cond_magnt<0, NP, NQ, NR>(BP);
		__bound_cond_magnt<1, NQ, NR, NP>(BQ);
		__bound_cond_magnt<2, NR, NP, NQ>(BR);
}
