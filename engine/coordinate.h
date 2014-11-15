#ifndef COORDINATE_H_
#define COORDINATE_H_

#include <assert.h>
#include <math.h>

#include "vector.h"
#include "../config/my_coordinate.h"

// Reset the components to viewing from d=0
template <int d>
static inline struct vector swap(struct vector position)
{
	switch (d % 3) {
	// (p, q, r) -> (p, q, r)
	case 0: return position;
	// (q, r, p) -> (p, q, r)
	case 1: return {position._2, position._0, position._1};
	// (r, p, q) -> (p, q, r)
	case 2: return {position._1, position._2, position._0};
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
template <int d>
vector othnormal_position(struct vector position)
{
	position = swap<d>(position);
	return xyz(position);
}

// also called jacobi matrix
// XXX: memoize
template <int c, int d>
struct vector covariant_basic_vector(struct vector position)
{
	position = swap<d>(position);
	switch (c%3) {
	case 0: return dxyz_dp(position);
	case 1: return dxyz_dq(position);
	case 2: return dxyz_dr(position);
	}
}

template <int ci, int cj, int d>
static inline double covariant_metric_tensor(struct vector position)
{
	struct vector gi = covariant_basic_vector<ci, d>(position);
	if (ci == cj) return inner_product(gi, gi);
	struct vector gj = covariant_basic_vector<cj, d>(position);
	return inner_product(gi, gj);
}

template <int d>
double jacobian(struct vector position)
{
	position = swap<d>(position);
	struct vector dp = dxyz_dp(position);
	struct vector dq = dxyz_dq(position);
	struct vector dr = dxyz_dr(position);

	return dp._0 * dq._1 * dr._2
		 + dq._0 * dr._1 * dp._2
		 + dr._0 * dp._1 * dq._2
		 - dr._0 * dq._1 * dp._2
		 - dq._0 * dp._1 * dr._2
		 - dp._0 * dr._1 * dq._2;
}

// also called inverse jacobi matrix
template <int c, int d>
struct vector contravariant_basic_vector(struct vector position)
{
	return vector_multiple(
			exterior_product(covariant_basic_vector<c+1, d>(position),
							 covariant_basic_vector<c+2, d>(position)),
			1.0/jacobian<d>(position));
}

template <int ci, int cj, int d>
static inline double contravariant_metric_tensor(struct vector position)
{
	struct vector gi = contravariant_basic_vector<ci, d>(position);
	if (ci == cj) return inner_product(gi, gi);
	struct vector gj = contravariant_basic_vector<cj, d>(position);
	return inner_product(gi, gj);
}

#endif /* COORDINATE_H_ */
