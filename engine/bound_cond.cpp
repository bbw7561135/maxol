#include "vector.h"
#include "../config/comp_param.h"
#include "../config/my_bound_cond.h"

// TODO: Some bugs might be in boundary condition.

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
 * case c == 0: bound(A,B,C,D,E,F) -> bound(0,1,2,3,4,5)
 * case c == 1: bound(A,B,C,D,E,F) -> bound(1,2,0,4,5,3)
 * case c == 2: bound(A,B,C,D,E,F) -> bound(2,0,1,5,3,4)
 */
template <int c>
static inline double elect_val_on_bound_B(struct vector pos, const double *E)
{
	pos = swap<c>(pos);

	switch(c%3){
	case 0:  return bound_cond_elect_p1(pos._2, pos._0, E);
	case 1:  return bound_cond_elect_q2(pos._0, pos._1, E);
	case 2:  return bound_cond_elect_r3(pos._1, pos._2, E);
	}
}

template <int c>
static inline double elect_val_on_bound_C(struct vector pos, const double *E)
{
	pos = swap<c>(pos);

	switch(c%3){
	case 0:  return bound_cond_elect_p2(pos._0, pos._1, E);
	case 1:  return bound_cond_elect_q0(pos._1, pos._2, E);
	case 2:  return bound_cond_elect_r1(pos._2, pos._0, E);
	}
}

template <int c>
static inline double elect_val_on_bound_E(struct vector pos, const double *E)
{
	pos = swap<c>(pos);

	switch(c%3){
	case 0: return bound_cond_elect_p4(pos._2, pos._0, E);
	case 1: return bound_cond_elect_q5(pos._0, pos._1, E);
	case 2: return bound_cond_elect_r3(pos._1, pos._2, E);
	}
}

template <int c>
static inline double elect_val_on_bound_F(struct vector pos, const double *E)
{
	pos = swap<c>(pos);

	switch(c%3){
	case 0: return bound_cond_elect_p5(pos._0, pos._1, E);
	case 1: return bound_cond_elect_q3(pos._1, pos._2, E);
	case 2: return bound_cond_elect_r4(pos._2, pos._0, E);
	}
}

template <int c, int N0, int N1, int N2>
static void __bound_cond_elect(double *E)
{
	// tangential components to boundary is zero

	for (int i = 0; i < N0; i++) {
		for (int j = 1; j < N1-1; j++) {
			// Bound C
			const int k0 = 0;
			const int l0 = j + N1*(k0 + N2*i);
			E[l0] = elect_val_on_bound_C<c>(
					{(double)i, (double)j, (double)k0}, E);

			// Bound F
			const int k1 = N2-1;
			const int l1 = j + N1*(k1 + N2*i);
			E[l1] = elect_val_on_bound_F<c>(
					{(double)i, (double)j, (double)k1}, E);
		}
		for (int k = 0; k < N2; k++) {
			// Bound B
			const int j0 = 0;
			const int l0 = j0 + N1*(k + N2*i);
			E[l0] = elect_val_on_bound_B<c>(
					{(double)i, (double)j0, (double)k}, E);

			// Bound E
			const int j1 = N1-1;
			const int l1 = j1 + N1*(k + N2*i);
			E[l1] = elect_val_on_bound_E<c>(
					{(double)i, (double)j1, (double)k}, E);
		}
	}
}

void bound_cond_elect(double *Ep, double *Eq, double *Er)
{
	__bound_cond_elect<0, NP-1, NQ, NR>(Ep);
	__bound_cond_elect<1, NQ-1, NR, NP>(Eq);
	__bound_cond_elect<2, NR-1, NP, NQ>(Er);
}

// Only vertical components face on the boundaries. However, they are not
// affect because light is transverse wave. Certainly, they not used for
// computing. I don't think this condition might be needed. I don't know why
// boundary condition in magnetic field is not required.

template <int c, int N0, int N1, int N2>
static void __bound_cond_magnt(double *B)
{
	// vertical component to boundary is zero
	for (int k = 0; k < N2; k++)
		for (int j = 0; j < N1; j++) {
			// Bound A
			const int i0 = 0;
			const int l0 = j + N1*(k + N2*i0);
			B[l0] = 0.0;

			// Bound D
			const int i1 = N0-1;
			const int l1 = j + N1*(k + N2*i1);
			B[l1] = 0.0;
		}
}

void bound_cond_magnt(double *BP, double *BQ, double *BR)
{
		__bound_cond_magnt<0, NP, NQ-1, NR-1>(BP);
		__bound_cond_magnt<1, NQ, NR-1, NP-1>(BQ);
		__bound_cond_magnt<2, NR, NP-1, NQ-1>(BR);
}
