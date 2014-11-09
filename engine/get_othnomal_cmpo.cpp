#include "config/comp_param.h"

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

void get_othnomal_cmpo_elect(
		float *EX, float *EY, float *EZ,
		const double *EP, const double *EQ, const double *ER)
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
	for (int k = 0; k < NR; k++)
	for (int j = 0; j < NQ; j++)
	for (int i = 0; i < NP; i++) {
		const double p = (double)i;
		const double q = (double)j;
		const double r = (double)k;

		const int l0 = j + NQ*(k + NR*(i-1));
		const int l1 = k + NR*(i + NP*(j-1));
		const int l2 = i + NP*(j + NQ*(k-1));

		// physical components along to the grids
		// interpolate with p direction
		double e0 = interpolate2<NP, NQ*NR>(EP, i, l0);
		// interpolate with q direction
		double e1 = interpolate2<NQ, NR*NP>(EQ, j, l1);
		// interpolate with r direction
		double e2 = interpolate2<NR, NP*NQ>(ER, k, l2);

		// Convert to unphysical component.
		e0 /= len_of_covariant_basic_vector<0, 0>(p, q, r);
		e1 /= len_of_covariant_basic_vector<1, 0>(p, q, r);
		e2 /= len_of_covariant_basic_vector<2, 0>(p, q, r);

		const int l = i + NP*(j + NQ*k);

		EX[l] = (float)(
				covariant_basic_vector<0, 0, 0>(p, q, r) * e0 +
				covariant_basic_vector<0, 1, 0>(p, q, r) * e1 +
				covariant_basic_vector<0, 2, 0>(p, q, r) * e2);

		EY[l] = (float)(
				covariant_basic_vector<1, 0, 0>(p, q, r) * e0 +
				covariant_basic_vector<1, 1, 0>(p, q, r) * e1 +
				covariant_basic_vector<1, 2, 0>(p, q, r) * e2);

		EZ[l] = (float)(
				covariant_basic_vector<2, 0, 0>(p, q, r) * e0 +
				covariant_basic_vector<2, 1, 0>(p, q, r) * e1 +
				covariant_basic_vector<2, 2, 0>(p, q, r) * e2);
	}
}

void get_othnomal_cmpo_magnt(
		float *BX, float *BY, float *BZ,
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
	for (int k = 0; k < NR; k++)
	for (int j = 0; j < NQ; j++)
	for (int i = 0; i < NP; i++) {
		const double p = (double)i;
		const double q = (double)j;
		const double r = (double)k;

		const int l0 = j-1 + (NQ-1)*(k-1 + (NR-1)*i);
		const int l1 = k-1 + (NR-1)*(i-1 + (NP-1)*j);
		const int l2 = i-1 + (NP-1)*(j-1 + (NQ-1)*k);

		// physical components vertical with grids
		// interpolate with q-r surface
		double b0 = interpolate4<NQ, NR, 1, NQ-1>(BP, j, k, l0);
		// interpolate with r-p surface
		double b1 = interpolate4<NR, NP, 1, NR-1>(BQ, k, i, l1);
		// interpolate with p-q surface
		double b2 = interpolate4<NP, NQ, 1, NP-1>(BR, i, j, l2);

		// Convert to unphysical component.
		b0 /= len_of_contravariant_basic_vector<0, 0>(p, q, r);
		b1 /= len_of_contravariant_basic_vector<1, 0>(p, q, r);
		b2 /= len_of_contravariant_basic_vector<2, 0>(p, q, r);

		const int l = i + NP*(j + NQ*k);

		BX[l] = (float)(
				contravariant_basic_vector<0, 0, 0>(p, q, r) * b0 +
				contravariant_basic_vector<0, 1, 0>(p, q, r) * b1 +
				contravariant_basic_vector<0, 2, 0>(p, q, r) * b2);

		BY[l] = (float)(
				contravariant_basic_vector<1, 0, 0>(p, q, r) * b0 +
				contravariant_basic_vector<1, 1, 0>(p, q, r) * b1 +
				contravariant_basic_vector<1, 2, 0>(p, q, r) * b2);

		BZ[l] = (float)(
				contravariant_basic_vector<2, 0, 0>(p, q, r) * b0 +
				contravariant_basic_vector<2, 1, 0>(p, q, r) * b1 +
				contravariant_basic_vector<2, 2, 0>(p, q, r) * b2);
	}
}
