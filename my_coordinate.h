/*
 * Rewrite for the coordinate system you want to use.
 */

#ifndef MY_COORDINATE_H_
#define MY_COORDINATE_H_

// Get reciprocal coordinates member.

// p, q, r -> x
static inline double x(double p, double q, double r)
{
	return p;
}

// p, q, r -> y
static inline double y(double p, double q, double r)
{
	return q;
}

// p, q, r -> z
static inline double z(double p, double q, double r)
{
	return r;
}

// covariant basic vectors

// p, q, r -> dx/dp
static inline double dx_dp(double p, double q, double r)
{
	return 1.0;
}

// p, q, r -> dx/dq
static inline double dx_dq(double p, double q, double r)
{
	return 0.0;
}

// p, q, r -> dx/dr
static inline double dx_dr(double p, double q, double r)
{
	return 0.0;
}

// p, q, r -> dy/dp
static inline double dy_dp(double p, double q, double r)
{
	return 0.0;
}

// p, q, r -> dy/dq
static inline double dy_dq(double p, double q, double r)
{
	return 1.0;
}

// p, q, r -> dy/dr
static inline double dy_dr(double p, double q, double r)
{
	return 0.0;
}

// p, q, r -> dz/dp
static inline double dz_dp(double p, double q, double r)
{
	return 0.0;
}

// p, q, r -> dz/dq
static inline double dz_dq(double p, double q, double r)
{
	return 0.0;
}

// p, q, r -> dz/dr
static inline double dz_dr(double p, double q, double r)
{
	return 1.0;
}

#endif /* MY_COORDINATE_H_ */
