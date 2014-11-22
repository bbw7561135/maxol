/*
 * Edit as you like.
 */

#ifndef MY_BOUND_COND_H_
#define MY_BOUND_COND_H_

#include <assert.h>

#include "../engine/coordinate.h"

#define P_UP   0
#define Q_UP   1
#define R_UP   2
#define P_DOWN 3
#define Q_DOWN 4
#define R_DOWN 5

#define CMPO_P (0 << 3)
#define CMPO_Q (1 << 3)
#define CMPO_R (2 << 3)

/*
 *                 / P_UP
 *         R_DOWN L
 *         +--|---- +
 *        /   v    /|
 *       /        / |
 * Q_UP +--------+ <--Q_DOWN ^ r
 *  --> |        |  +        |
 *      | P_DOWN | /         +---> q
 *      |        |/         /
 *      +--------+         L p
 *            ^
 *            | R_UP
 */
#if 0
template <int c, int b>
static inline double my_bound_cond_elect(int p, int q, int r, double *E)
{
	// All boundaries are zero
	switch (c << 3 | b) {

	case  CMPO_P | P_UP: break;
	case  CMPO_Q | P_UP: return 0.0;
	case  CMPO_R | P_UP: return 0.0;

	case  CMPO_P | P_DOWN: break;
	case  CMPO_Q | P_DOWN: return 0.0;
	case  CMPO_R | P_DOWN: return 0.0;

	case  CMPO_P | Q_UP: return 0.0;
	case  CMPO_Q | Q_UP: break;
	case  CMPO_R | Q_UP: return 0.0;

	case  CMPO_P | Q_DOWN: return 0.0;
	case  CMPO_Q | Q_DOWN: break;
	case  CMPO_R | Q_DOWN: return 0.0;

	case  CMPO_P | R_UP: return 0.0;
	case  CMPO_Q | R_UP: return 0.0;
	case  CMPO_R | R_UP: break;

	case  CMPO_P | R_DOWN: return 0.0;
	case  CMPO_Q | R_DOWN: return 0.0;
	case  CMPO_R | R_DOWN: break;
	}
	assert(0);
}
#endif

#if 1
template <int c, int b>
static inline double my_bound_cond_elect(int p, int q, int r, double *E)
{
	// All boundaries are periodic.
	switch (c << 3 | b) {

	case  CMPO_P | P_UP: break;
	case  CMPO_Q | P_UP: return *elect_field<1, NP, NQ, NR>(NP-2, q, r, E);
	case  CMPO_R | P_UP: return *elect_field<2, NP, NQ, NR>(NP-2, q, r, E);

	case  CMPO_P | P_DOWN: break;
	case  CMPO_Q | P_DOWN: return *elect_field<1, NP, NQ, NR>(1, q, r, E);
	case  CMPO_R | P_DOWN: return *elect_field<2, NP, NQ, NR>(1, q, r, E);

	case  CMPO_P | Q_UP: return *elect_field<0, NP, NQ, NR>(p, NQ-2, r, E);
	case  CMPO_Q | Q_UP: break;
	case  CMPO_R | Q_UP: return *elect_field<2, NP, NQ, NR>(p, NQ-2, r, E);

	case  CMPO_P | Q_DOWN: return *elect_field<0, NP, NQ, NR>(p, 1, r, E);
	case  CMPO_Q | Q_DOWN: break;
	case  CMPO_R | Q_DOWN: return *elect_field<2, NP, NQ, NR>(p, 1, r, E);

	case  CMPO_P | R_UP: return *elect_field<0, NP, NQ, NR>(p, q, NR-2, E);
	case  CMPO_Q | R_UP: return *elect_field<1, NP, NQ, NR>(p, q, NR-2, E);
	case  CMPO_R | R_UP: break;

	case  CMPO_P | R_DOWN: return *elect_field<0, NP, NQ, NR>(p, q, 1, E);
	case  CMPO_Q | R_DOWN: return *elect_field<0, NP, NQ, NR>(p, q, 1, E);
	case  CMPO_R | R_DOWN: break;
	}
	assert(0);
}
#endif

#endif /* MY_BOUND_COND_H_ */
