#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>

#include "coordinate.h"
#include "../config/comp_param.h"

// Linear-extrapolate if on boundary else linear-interpolate.
template <int N, int NS>
double interpolate_line(const double *A, int i, int l)
{
	switch (i) {
	case 0:   return   1.5 * A[l+NS] - 0.5 * A[l+2*NS];
	case N-1: return - 0.5 * A[l-NS] + 1.5 * A[l];
	default:  return   0.5 * A[l]    + 0.5 * A[l+NS];
	}
}

static void get_othnomal_cmpo_elect(
		float *EX, float *EY, float *EZ,
		const double *EP, const double *EQ, const double *ER)
{
	for (int k = 0; k < NR; k++)
	for (int j = 0; j < NQ; j++)
	for (int i = 0; i < NP; i++) {
		const int l0 = j + NQ*(k + NR*(i-1));
		const int l1 = k + NR*(i + NP*(j-1));
		const int l2 = i + NP*(j + NQ*(k-1));

		const double p = (double)i;
		const double q = (double)j;
		const double r = (double)k;

		// physical components along to the grids
		double e0 = interpolate_line<NP, NQ*NR>(EP, i, l0);
		double e1 = interpolate_line<NQ, NR*NP>(EQ, j, l1);
		double e2 = interpolate_line<NR, NP*NQ>(ER, k, l2);

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

// Linear-extrapolate if on boundary else linear-interpolate.
template <int N1, int N2, int NS1, int NS2>
double interpolate_plane(const double *A, int i, int j, int l)
{
	double a, b;

	switch(i) {
	case 0:
		a = interpolate_line<N2, NS2>(A, j, l+  NS1);
		b = interpolate_line<N2, NS2>(A, j, l+2*NS1);
		return   1.5*a - 0.5*b;
	case N1-1:
		a = interpolate_line<N2, NS2>(A, j, l-NS1);
		b = interpolate_line<N2, NS2>(A, j, l);
		return - 0.5*a + 1.5*b;
	default:
		a = interpolate_line<N2, NS2>(A, j, l);
		b = interpolate_line<N2, NS2>(A, j, l+NS1);
		return 0.5*a + 0.5*b;
	}
}

static void get_othnomal_cmpo_magnt(
		float *BX, float *BY, float *BZ,
		const double *BP, const double *BQ, const double *BR)
{
	for (int k = 0; k < NR; k++)
	for (int j = 0; j < NQ; j++)
	for (int i = 0; i < NP; i++) {
		const int l0 = j-1 + (NQ-1)*(k-1 + (NR-1)*i);
		const int l1 = k-1 + (NR-1)*(i-1 + (NP-1)*j);
		const int l2 = i-1 + (NP-1)*(j-1 + (NQ-1)*k);

		const double p = (double)i;
		const double q = (double)j;
		const double r = (double)k;

		// physical components vertical with grids
		double b0 = interpolate_plane<NQ, NR, 1, NQ-1>(BP, j, k, l0);
		double b1 = interpolate_plane<NR, NP, 1, NR-1>(BQ, k, i, l1);
		double b2 = interpolate_plane<NP, NQ, 1, NP-1>(BR, i, j, l2);

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

int output(const double *Ep, const double *Eq, const double *Er,
           const double *BP, const double *BQ, const double *BR,
           int nt)
{
	// The order of values is as (0,0,0)(1,0,0)(0,1,0)(1,1,0)(0,0,1)...
	float *EBX = (float *)malloc(NP*NQ*NR*sizeof(float));
	float *EBY = (float *)malloc(NP*NQ*NR*sizeof(float));
	float *EBZ = (float *)malloc(NP*NQ*NR*sizeof(float));
	if (EBX == NULL || EBY == NULL || EBZ == NULL ) {
		free(EBX); free(EBY); free(EBZ);
		return ENOMEM;
	}

	char *outpath = getenv("MAXOL_OUT_PATH");
	if (outpath == NULL) {
		outpath = ".";
	}

	char fpath[1024];
	int fd;
	int _errno;
	char header[] = "### Maxol ###\n";
	int asize[] = {NP, NQ, NR};
	float t;

	// Output electric field

	if (sprintf(fpath, "%s/%08dE", outpath, nt) == EOF) {
		free(EBX); free(EBY); free(EBZ);
		return ENAMETOOLONG;
	}

	if ((fd = open(fpath, O_WRONLY | O_CREAT,
			S_IFREG | S_IRUSR | S_IRGRP | S_IROTH)) == -1)
		goto err;

	if (write(fd, header, sizeof(header)) == -1)
		goto err;

	if (write(fd, asize, 3*sizeof(int)) == -1)
		goto err;

	if (write(fd, &nt, sizeof(int)) == -1)
		goto err;

	t = ((float)nt*DT);
	if (write(fd, &t, sizeof(float)) == -1)
		goto err;

	get_othnomal_cmpo_elect(EBX, EBY, EBZ, Ep, Eq, Er);

	if (write(fd, EBX, sizeof(float)*NP*NQ*NR) == -1)
		goto err;

	if (write(fd, EBY, sizeof(float)*NP*NQ*NR) == -1)
		goto err;

	if (write(fd, EBZ, sizeof(float)*NP*NQ*NR) == -1)
		goto err;

	close(fd);

	// Output magnetic flux densities

	if (sprintf(fpath, "%s/%08dB", outpath, nt) == EOF) {
		free(EBX); free(EBY); free(EBZ);
		return ENAMETOOLONG;
	}

	if ((fd = open(fpath, O_WRONLY | O_CREAT,
			S_IFREG | S_IRUSR | S_IRGRP | S_IROTH)) == -1)
		goto err;

	if (write(fd, header, sizeof(header)) == -1)
		goto err;

	if (write(fd, asize, 3*sizeof(int)) == -1)
		goto err;

	t = ((float)nt-0.5)*DT;
	if (write(fd, &t, sizeof(float)) == -1)
		goto err;

	get_othnomal_cmpo_magnt(EBX, EBY, EBZ, BP, BQ, BR);

	if (write(fd, EBX, sizeof(float)*NP*NQ*NR) == -1)
		goto err;

	if (write(fd, EBY, sizeof(float)*NP*NQ*NR) == -1)
		goto err;

	if (write(fd, EBZ, sizeof(float)*NP*NQ*NR) == -1)
		goto err;

	close(fd);
	free(EBX); free(EBY); free(EBZ);
	return 0;

err:
	_errno = errno;
	free(EBX); free(EBY); free(EBZ);
	close(fd);
	return _errno;
}
