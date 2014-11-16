/*
 * Edit as you like.
 */

#ifndef MY_BOUND_COND_H_
#define MY_BOUND_COND_H_

#include "../engine/coordinate.h"

/*
 *                / 0
 *           5   L
 *        +--|---- +
 *       /   v    /|
 *      /        / |
 * 1   +--------+ <--4    ^ r
 * --> |        |  +      |
 *     |   3    | /       +---> q
 *     |        |/       /
 *     +--------+       L p
 *           ^
 *           | 2
 */

// Bound 0 - q & r component

static inline double bound_cond_elect_q0(int q, int r, const double *E)
{
	const int l = r + NR*(0 + NP*q);
/*
	// Dirichlet
	return 0.0;

	// Naumann. Assume grid is vertical to boundary.
	const double grad = 0.0;
	return (- 2.0*grad/len_of_contravariant_basic_vector<0, 0>(0, q, r) +
			4*E[l+NR] - E[l+2*NR]) * (1.0/3.0);
*/
	// Periodic
	return E[l + NR*(NP-2)];
}

static inline double bound_cond_elect_r0(int q, int r, const double *E)
{
	const int l = 0 + NP*(q + NQ*r);
/*
	// Dirichlet
	return 0.0;

	// Naumann. Assume grid is vertical to boundary.
	const double grad = 0.0;
	return (- 2.0*grad/len_of_contravariant_basic_vector<0, 0>(0, q, r) +
			4*E[l+1] - E[l+2]) * (1.0/3.0);
*/
	// Periodic
	return E[l + (NP-2)];
}

// Bound 3 - q & r component

static inline double bound_cond_elect_q3(int q, int r, const double *E)
{
	const int l = r + NR*((NP-1) + NP*q);
/*
	// Dirichlet
	return 0.0;

	// Naumann. Assume grid is vertical to boundary.
	const double grad = 0.0;
	return (2.0*grad/len_of_contravariant_basic_vector<0, 0>(NP-1, q, r) +
			4*E[l-NR] - E[l-2*NR]) * (1.0/3.0);
*/
	// Periodic
	return E[l - NR*(NP-2)];
}

static inline double bound_cond_elect_r3(int q, int r, const double *E)
{
	const int l = (NP-1) + NP*(q + NQ*r);
/*
	// Dirichlet
	return 0.0;

	// Naumann. Assume grid is vertical to boundary.
	const double grad = 0.0;
	return (2.0*grad/len_of_contravariant_basic_vector<0, 0>(NP-1, q, r) +
			4*E[l-1] - E[l-2]) * (1.0/3.0);
*/
	// Periodic
	return E[l - (NP-2)];
}

// Bound 1 - r & p component

static inline double bound_cond_elect_r1(int r, int p, const double *E)
{
	const int l = p + NP*(0 + NQ*r);
/*
	// Dirichlet
	return 0.0;

	// Naumann. Assume grid is vertical to boundary.
	const double grad = 0.0;
	return (- 2.0*grad/len_of_contravariant_basic_vector<1, 0>(p, 0, r) +
			4*E[l+NP] - E[l+2*NP]) * (1.0/3.0);
*/
	// Periodic
	return E[l + NP*(NQ-2)];
}

static inline double bound_cond_elect_p1(int r, int p, const double *E)
{
	const int l = 0 + NQ*(r + NR*p);
/*
	// Dirichlet
	return 0.0;

	// Naumann. Assume grid is vertical to boundary.
	const double grad = 0.0;
	return (- 2.0*grad/len_of_contravariant_basic_vector<1, 0>(p, 0, r) +
			4*E[l+1] - E[l+2]) * (1.0/3.0);
*/
	// Periodic
	return E[l + (NQ-2)];
}

// Bound 4 - r & p component

static inline double bound_cond_elect_r4(int r, int p, const double *E)
{
	const int l = p + NP*((NQ-1) + NQ*r);
/*
	// Dirichlet
	return 0.0;

	// Naumann. Assume grid is vertical to boundary.
	const double grad = 0.0;
	return (2.0*grad/len_of_contravariant_basic_vector<1, 0>(p, NQ-1, r) +
			4*E[l-NP] - E[l-2*NP]) * (1.0/3.0);
*/
	// Periodic
	return E[l - NP*(NQ-2)];
}

static inline double bound_cond_elect_p4(int r, int p, const double *E)
{
	const int l = (NQ-1) + NQ*(r + NR*p);
/*
	// Dirichlet
	return 0.0;

	// Naumann. Assume grid is vertical to boundary.
	const double grad = 0.0;
	return (2.0*grad/len_of_contravariant_basic_vector<1, 0>(p, NQ-1, r) +
			4*E[l-1] - E[l-2]) * (1.0/3.0);
*/
	// Periodic
	return E[l - (NQ-2)];
}

// Bound 2 - p & q component

static inline double bound_cond_elect_p2(int p, int q, const double *E)
{
	const int l = q + NQ*(0 + NR*p);
/*
	// Dirichlet
	return 0.0;

	// Naumann. Assume grid is vertical to boundary.
	const double grad = 0.0;
	return (- 2.0*grad/len_of_contravariant_basic_vector<2, 0>(p, q, 0) +
			4*E[l+NQ] - E[l+2*NQ]) * (1.0/3.0);
*/
	// Periodic
	return E[l + NQ*(NR-2)];
}

static inline double bound_cond_elect_q2(int p, int q, const double *E)
{
	const int l = 0 + NR*(p + NP*q);
/*
	// Dirichlet
	return 0.0;

	// Naumann. Assume grid is vertical to boundary.
	const double grad = 0.0;
	return (- 2.0*grad/len_of_contravariant_basic_vector<2, 0>(p, q, 0) +
			4*E[l+1] - E[l+2]) * (1.0/3.0);
*/
	// Periodic
	return E[l + (NR-2)];
}

// Bound 5 - p & q component

static inline double bound_cond_elect_p5(int p, int q, const double *E)
{
	const int l = q + NQ*(NR-1 + NR*p);
/*
	// Dirichlet
	return 0.0;

	// Naumann. Assume grid is vertical to boundary.
	const double grad = 0.0;
	return (2.0*grad/len_of_contravariant_basic_vector<2, 0>(p, q, NR-1) +
			4*E[l-NQ] - E[l-2*NQ]) * (1.0/3.0);
*/
	// Periodic
	return E[l - NQ*(NR-2)];
}

static inline double bound_cond_elect_q5(int p, int q, const double *E)
{
	const int l = (NR-1) + NR*(p + NP*q);
/*
	// Dirichlet
	return 0.0;

	// Naumann. Assume grid is vertical to boundary.
	const double grad = 0.0;
	return (2.0*grad/len_of_contravariant_basic_vector<2, 0>(p, q, NR-1) +
			4*E[l-1] - E[l-2]) * (1.0/3.0);
*/
	// Periodic
	return E[l - (NR-2)];
}

#endif /* MY_BOUND_COND_H_ */
