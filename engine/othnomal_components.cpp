#include "coordinate.h"
#include "vector.h"

// Certainly a conditional branch in for loop is not fast,
// but use it for readablity.

// Linear-extrapolate if on boundary else linear-interpolate.
template <int N, int NS>
static double interpolate2(const double *A, int i, int l)
{
	switch (i) {
	// boundary
	case 0:   return   1.5 * A[l+NS] - 0.5 * A[l+2*NS];
	case N-1: return - 0.5 * A[l-NS] + 1.5 * A[l];
	// non-boundary
	default:  return   0.5 * A[l]    + 0.5 * A[l+NS];
	}
}

// Linear-extrapolate if on boundary else linear-interpolate.
template <int N1, int N2, int NS1, int NS2>
static double interpolate4(const double *A, int i1, int i2, int l)
{
	double a, b;
	switch(i1) {
	// boundary
	case 0:
		a = interpolate2<N2, NS2>(A, i2, l+NS1);
		b = interpolate2<N2, NS2>(A, i2, l+2*NS1);
		return   1.5*a - 0.5*b;
	case N1-1:
		a = interpolate2<N2, NS2>(A, i2, l-NS1);
		b = interpolate2<N2, NS2>(A, i2, l);
		return - 0.5*a + 1.5*b;
	// non-boundary
	default:
		a = interpolate2<N2, NS2>(A, i2, l);
		b = interpolate2<N2, NS2>(A, i2, l+NS1);
		return   0.5*a + 0.5*b;
	}
}

// This function is usually called in a for-loop. Thus, the overhead of calling
// might not be negligible small.
struct vector othnomal_electric_field(int i, int j, int k,
		const double *Ep, const double *Eq, const double *Er)
{
	/*
	 *         p   +---- Calculate electric field vector here.
	 *   +-----+--/--+
	 *   |     | /   |
	 *   |     v/    |
	 * q +-->  +     +-->
	 *   |           |
	 *   |           |
	 *   +-----+-----+
	 *         |
	 *         v
	 */
	const struct vector pos = {(double)i, (double)j, (double)k};

	const int l0 = j + NQ*(k + NR*(i-1));
	const int l1 = k + NR*(i + NP*(j-1));
	const int l2 = i + NP*(j + NQ*(k-1));

	struct vector e;

	// physical components along to the grids
	// interpolate with p direction
	e._0 = interpolate2<NP, NQ*NR>(Ep, i, l0);
	// interpolate with q direction
	e._1 = interpolate2<NQ, NR*NP>(Eq, j, l1);
	// interpolate with r direction
	e._2 = interpolate2<NR, NP*NQ>(Er, k, l2);

	const struct vector dp = covariant_basic_vector<0, 0>(pos);
	const struct vector dq = covariant_basic_vector<1, 0>(pos);
	const struct vector dr = covariant_basic_vector<2, 0>(pos);

	// Convert to unphysical component.
	e._0 /= vector_length(dp);
	e._1 /= vector_length(dq);
	e._2 /= vector_length(dr);

	const struct vector dx = {dp._0, dq._0, dr._0};
	const struct vector dy = {dp._1, dq._1, dr._1};
	const struct vector dz = {dp._2, dq._2, dr._2};

	// coordinate transformation
	return {inner_product(e, dx), inner_product(e, dy), inner_product(e, dz)};
}

struct vector othnomal_magnetic_flux(int i, int j, int k,
		const double *BP, const double *BQ, const double *BR)
{
	/*
	 *         p
	 *         +-->
	 *   +-----------+
	 *   |       +------ Calculate magnetic flux density vector here.
	 *   |      /    |
	 * q+|     +     |+
	 *  ||           ||
	 *  v|           |v
	 *   +-----------+
	 *         +-->
	 */
	const struct vector pos = {(double)i, (double)j, (double)k};

	const int l0 = j-1 + (NQ-1)*(k-1 + (NR-1)*i);
	const int l1 = k-1 + (NR-1)*(i-1 + (NP-1)*j);
	const int l2 = i-1 + (NP-1)*(j-1 + (NQ-1)*k);

	struct vector b;
	// physical components vertical with grids
	// interpolate with q-r surface
	b._0 = interpolate4<NQ, NR, 1, NQ-1>(BP, j, k, l0);
	// interpolate with r-p surface
	b._1 = interpolate4<NR, NP, 1, NR-1>(BQ, k, i, l1);
	// interpolate with p-q surface
	b._2 = interpolate4<NP, NQ, 1, NP-1>(BR, i, j, l2);

	const struct vector dp = contravariant_basic_vector<0, 0>(pos);
	const struct vector dq = contravariant_basic_vector<1, 0>(pos);
	const struct vector dr = contravariant_basic_vector<2, 0>(pos);

	// Convert to unphysical component.
	b._0 /= vector_length(dp);
	b._1 /= vector_length(dq);
	b._2 /= vector_length(dr);

	const struct vector dx = {dp._0, dq._0, dr._0};
	const struct vector dy = {dp._1, dq._1, dr._1};
	const struct vector dz = {dp._2, dq._2, dr._2};

	// coordinate transformation
	return {inner_product(b, dx), inner_product(b, dy), inner_product(b, dz)};
}
