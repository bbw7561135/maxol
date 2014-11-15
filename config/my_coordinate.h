/*
 * Rewrite for the coordinate system you want to use.
 */

#ifndef MY_COORDINATE_H_
#define MY_COORDINATE_H_

#include "comp_param.h"

// Get reciprocal coordinates member.

// p, q, r -> x, y, z
static inline struct vector xyz(struct vector pqr)
{
	double x = (1.0/(double)(NP-2)) * pqr._0 + (1.0/(double)(NQ-2)) * pqr._1;
	double y = (1.0/(double)(NQ-2)) * pqr._1 - (1.0/(double)(NP-2)) * pqr._0;
	double z = (1.0/(double)(NR-2)) * pqr._2;

	return {x, y, z};
}

// covariant basic vectors

// p, q, r -> dx/dp, dy/dp, dz/dp
static inline struct vector dxyz_dp(struct vector pqr)
{
	return {1.0/(double)(NP-2),
			- (1.0/(double)(NP-2)),
			0.0};
}

// p, q, r -> dx/dq, dy/dq, dz/dq
static inline struct vector dxyz_dq(struct vector pqr)
{
	return {1.0/(double)(NQ-2),
			1.0/(double)(NQ-2),
			0.0};
}

// p, q, r -> dx/dr, dy/dr, dz/dr
static inline struct vector dxyz_dr(struct vector pqr)
{
	return {0.0,
			0.0,
			1.0/(double)(NR-2)};
}

#endif /* MY_COORDINATE_H_ */
