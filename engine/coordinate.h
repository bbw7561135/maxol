#ifndef COORDINATE_H_
#define COORDINATE_H_

#include <assert.h>
#include <math.h>

#include "../config/my_coordinate.h"
#include "../config/comp_param.h"

// Reset the components to viewing from d=0
template <int d, typename T>
static inline void swap(T *p, T *q, T *r)
{
	T tmp;

	switch (d % 3) {
	case 0: // (p, q, r) -> (p, q, r)
		break;
	case 1: // (q, r, p) -> (p, q, r)
		tmp = *p; *p = *r; *r = *q; *q = tmp;
		break;
	case 2: // (r, p, q) -> (p, q, r)
		tmp = *p; *p = *q; *q = *r; *r = tmp;
		break;
	}
}

/*
 * Return othnormal position vector at (p, q, r)
 *
 * c: component
 * d: direction
 *
 * d means viewing direction. (p, q, r) are components when viewing from d.
 * If d == 0, the view is normal.
 */
template <int c, int d, typename T> double othnormal_position(T p, T q, T r)
{
	swap<d>(&p, &q, &r);
	switch (c%3) {
	case 0: return x((double)p, (double)q, (double)r);
	case 1: return y((double)p, (double)q, (double)r);
	case 2: return z((double)p, (double)q, (double)r);
	}
}

// also called jacobi matrix
// XXX: memoize
template <int cx, int cp, int d>
double covariant_basic_vector(double p, double q, double r)
{
	swap<d>(&p, &q, &r);
	switch (cx%3 + 3*(cp%3)) {
	// p
	case 0: return dx_dp(p, q, r);
	case 1: return dy_dp(p, q, r);
	case 2: return dz_dp(p, q, r);
	// q
	case 3: return dx_dq(p, q, r);
	case 4: return dy_dq(p, q, r);
	case 5: return dz_dq(p, q, r);
	// r
	case 6: return dx_dr(p, q, r);
	case 7: return dy_dr(p, q, r);
	case 8: return dz_dr(p, q, r);
	}
}

template <int ci, int cj, int d>
static inline double covariant_metric_tensor(double p, double q, double r)
{
	swap<d>(&p, &q, &r);

	if (ci == cj) {
		const double v0 = covariant_basic_vector<0, ci, 0>(p, q, r);
		const double v1 = covariant_basic_vector<1, ci, 0>(p, q, r);
		const double v2 = covariant_basic_vector<2, ci, 0>(p, q, r);
		return v0*v0 + v1*v1 + v2*v2;
	}

	return covariant_basic_vector<0, ci, 0>(p, q, r) *
		   covariant_basic_vector<0, cj, 0>(p, q, r) +
	       covariant_basic_vector<1, ci, 0>(p, q, r) *
		   covariant_basic_vector<1, cj, 0>(p, q, r) +
	       covariant_basic_vector<2, ci, 0>(p, q, r) *
		   covariant_basic_vector<2, cj, 0>(p, q, r);
}

// The length of the "c" component of covariant basic vector.
template <int c, int d>
double len_of_covariant_basic_vector(double p, double q, double r)
{
	return sqrt(covariant_metric_tensor<c, c, d>(p, q, r));
}

template <int d>
double jacobian(double p, double q, double r)
{
	swap<d>(&p, &q, &r);
	return dx_dp(p, q, r)*dy_dq(p, q, r)*dz_dr(p, q, r)
	     + dx_dq(p, q, r)*dy_dr(p, q, r)*dz_dp(p, q, r)
	     + dx_dr(p, q, r)*dy_dp(p, q, r)*dz_dq(p, q, r)
	     - dx_dr(p, q, r)*dy_dq(p, q, r)*dz_dp(p, q, r)
	     - dx_dq(p, q, r)*dy_dp(p, q, r)*dz_dr(p, q, r)
	     - dx_dp(p, q, r)*dy_dr(p, q, r)*dz_dq(p, q, r);
}

// contravariant basic vectors

// p, q, r -> dp/dx
static inline double dp_dx(double p, double q, double r)
{
	return (dy_dq(p, q, r)*dz_dr(p, q, r) -
			dy_dr(p, q, r)*dz_dq(p, q, r)) /
			jacobian<0>(p, q, r);
}

// p, q, r -> dq/dx
static inline double dq_dx(double p, double q, double r)
{
	return (dy_dr(p, q, r)*dz_dp(p, q, r) -
			dy_dp(p, q, r)*dz_dr(p, q, r)) /
			jacobian<0>(p, q, r);
}

// p, q, r -> dr/dx
static inline double dr_dx(double p, double q, double r)
{
	return (dy_dp(p, q, r)*dz_dq(p, q, r) -
			dy_dq(p, q, r)*dz_dp(p, q, r)) /
			jacobian<0>(p, q, r);
}

// p, q, r -> dp/dy
static inline double dp_dy(double p, double q, double r)
{
	return (dz_dq(p, q, r)*dx_dr(p, q, r) -
			dz_dr(p, q, r)*dx_dq(p, q, r)) /
			jacobian<0>(p, q, r);
}

// p, q, r -> dq/dy
static inline double dq_dy(double p, double q, double r)
{
	return (dz_dr(p, q, r)*dx_dp(p, q, r) -
			dz_dp(p, q, r)*dx_dr(p, q, r)) /
			jacobian<0>(p, q, r);
}

// p, q, r -> dr/dy
static inline double dr_dy(double p, double q, double r)
{
	return (dz_dp(p, q, r)*dx_dq(p, q, r) -
			dz_dq(p, q, r)*dx_dp(p, q, r)) /
			jacobian<0>(p, q, r);
}

// p, q, r -> dp/dz
static inline double dp_dz(double p, double q, double r)
{
	return (dx_dq(p, q, r)*dy_dr(p, q, r) -
			dx_dr(p, q, r)*dy_dq(p, q, r)) /
			jacobian<0>(p, q, r);
}

// p, q, r -> dq/dz
static inline double dq_dz(double p, double q, double r)
{
	return (dx_dr(p, q, r)*dy_dp(p, q, r) -
			dx_dp(p, q, r)*dy_dr(p, q, r)) /
			jacobian<0>(p, q, r);
}

// p, q, r -> dr/dz
static inline double dr_dz(double p, double q, double r)
{
	return (dx_dp(p, q, r)*dy_dq(p, q, r) -
			dx_dq(p, q, r)*dy_dp(p, q, r)) /
			jacobian<0>(p, q, r);
}

// also called inverse jacobi matrix
template <int cx, int cp, int d>
double contravariant_basic_vector(double p, double q, double r)
{
	swap<d>(&p, &q, &r);
	switch (cx%3 + 3*(cp%3)) {
	// p
	case 0: return dp_dx(p, q, r);
	case 1: return dp_dy(p, q, r);
	case 2: return dp_dz(p, q, r);
	// q
	case 3: return dq_dx(p, q, r);
	case 4: return dq_dy(p, q, r);
	case 5: return dq_dz(p, q, r);
	// r
	case 6: return dr_dx(p, q, r);
	case 7: return dr_dy(p, q, r);
	case 8: return dr_dz(p, q, r);
	}
}

template <int ci, int cj, int d>
static inline double contravariant_metric_tensor(double p, double q, double r)
{
	swap<d>(&p, &q, &r);

	if (ci == cj) {
		const double v0 = contravariant_basic_vector<0, ci, 0>(p, q, r);
		const double v1 = contravariant_basic_vector<1, ci, 0>(p, q, r);
		const double v2 = contravariant_basic_vector<2, ci, 0>(p, q, r);
		return v0*v0 + v1*v1 + v2*v2;
	}

	return contravariant_basic_vector<0, ci, 0>(p, q, r) *
	       contravariant_basic_vector<0, cj, 0>(p, q, r) +
	       contravariant_basic_vector<1, ci, 0>(p, q, r) *
	       contravariant_basic_vector<1, cj, 0>(p, q, r) +
	       contravariant_basic_vector<2, ci, 0>(p, q, r) *
	       contravariant_basic_vector<2, cj, 0>(p, q, r);
}

// The length of the "c" component of contravariant basic vector.
template <int c, int d>
double len_of_contravariant_basic_vector(double p, double q, double r)
{
	return sqrt(contravariant_metric_tensor<c, c, d>(p, q, r));
}

// p, q, r -> dx/dP
static inline double dx_dP(double p, double q, double r)
{
	return pow(jacobian<0>(p, q, r), 2.0/3.0) * dp_dx(p, q, r);
}

// p, q, r -> dy/dP
static inline double dy_dP(double p, double q, double r)
{
	return pow(jacobian<0>(p, q, r), 2.0/3.0) * dp_dy(p, q, r);
}

// p, q, r -> dz/dP
static inline double dz_dP(double p, double q, double r)
{
	return pow(jacobian<0>(p, q, r), 2.0/3.0) * dp_dz(p, q, r);
}

// p, q, r -> dx/dQ
static inline double dx_dQ(double p, double q, double r)
{
	return pow(jacobian<0>(p, q, r), 2.0/3.0) * dq_dx(p, q, r);
}

// p, q, r -> dy/dQ
static inline double dy_dQ(double p, double q, double r)
{
	return pow(jacobian<0>(p, q, r), 2.0/3.0) * dq_dy(p, q, r);
}

// p, q, r -> dz/dQ
static inline double dz_dQ(double p, double q, double r)
{
	return pow(jacobian<0>(p, q, r), 2.0/3.0) * dq_dz(p, q, r);
}

// p, q, r -> dx/dR
static inline double dx_dR(double p, double q, double r)
{
	return pow(jacobian<0>(p, q, r), 2.0/3.0) * dr_dx(p, q, r);
}

// p, q, r -> dy/dR
static inline double dy_dR(double p, double q, double r)
{
	return pow(jacobian<0>(p, q, r), 2.0/3.0) * dr_dy(p, q, r);
}

// p, q, r -> dz/dR
static inline double dz_dR(double p, double q, double r)
{
	return pow(jacobian<0>(p, q, r), 2.0/3.0) * dr_dz(p, q, r);
}
/*
 * contravariant basic vectors normalized to the same size as
 * covariant basic vectors:
 * dx->/dp * (dx->/dQ X dx->/dR) == dx->/dp * (dx->/dq X dx->/dr)
 * Although this function is heavy, I use it for simplification the code.
 * XXX: memoize
 */
template <int cx, int cp, int d>
double normalized_contravariant_basic_vector(double p, double q, double r)
{
	swap<d>(&p, &q, &r);
	switch (cx%3 + 3*(cp%3)) {
	// P
	case 0: return dx_dP(p, q, r);
	case 1: return dy_dP(p, q, r);
	case 2: return dz_dP(p, q, r);
	// Q
	case 3: return dx_dQ(p, q, r);
	case 4: return dy_dQ(p, q, r);
	case 5: return dz_dQ(p, q, r);
	// R
	case 6: return dx_dR(p, q, r);
	case 7: return dy_dR(p, q, r);
	case 8: return dz_dR(p, q, r);
	}
}

template <int c, int d>
double len_of_normalized_contravariant_basic_vector(
		double p, double q, double r)
{
	swap<d>(&p, &q, &r);

	const double dx = normalized_contravariant_basic_vector<0, c, 0>(p, q, r);
	const double dy = normalized_contravariant_basic_vector<1, c, 0>(p, q, r);
	const double dz = normalized_contravariant_basic_vector<2, c, 0>(p, q, r);

	return sqrt(dx*dx + dy*dy + dz*dz);
}

#endif /* COORDINATE_H_ */
