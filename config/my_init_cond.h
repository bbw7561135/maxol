/*
 * Rewrite for your system.
 */

#ifndef MY_INIT_COND_H_
#define MY_INIT_COND_H_

#include <math.h>

#include "comp_param.h"
#include "phys_param.h"

extern double dt;

// At nt = 0

static struct vector init_cond_elect(struct vector position)
{
	double ex = sin(- 2.0 * M_PI * position._2);
	double ey = 0.0;
	double ez = 0.0;

	return {ex, ey, ez};
}

// At nt = 0.5

static struct vector init_cond_magnt(struct vector position)
{
	double bx = 0.0;
	double by = sin(2.0 * M_PI * (0.5 * dt - position._2));
	double bz = 0.0;

	return {bx, by, bz};
}

#endif /* MY_INIT_COND_H_ */
