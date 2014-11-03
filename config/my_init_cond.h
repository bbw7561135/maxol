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

static double init_cond_elect_x(double x, double y, double z)
{
	return sin(- 2.0 * M_PI * z );
}

static double init_cond_elect_y(double x, double y, double z)
{
	return 0.0;
}

static double init_cond_elect_z(double x, double y, double z)
{
	return 0.0;
}

// At nt = 0.5

static double init_cond_magnt_x(double x, double y, double z)
{
	return 0.0;
}

static double init_cond_magnt_y(double x, double y, double z)
{
	return sin(2.0 * M_PI * (0.5 * dt - z));
}

static double init_cond_magnt_z(double x, double y, double z)
{
	return 0.0;
}

#endif /* MY_INIT_COND_H_ */
