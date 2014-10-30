#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>

#include "comp_param.h"
#include "coordinate.h"

// Linear-extrapolate if on boundary else linear-interpolate.
template <int N, int NS>
double interpolate_line(const double *A, int i, int l)
{
	switch (i) {
	case 0:   return   1.5*A[l+NS] - 0.5*A[l+2*NS];
	case N-1: return - 0.5*A[l-NS] + 1.5*A[l];
	default:  return 0.5*(A[l] + A[l+NS]);
	}
}

template <int c>
void get_othnomal_component_electric(
		float *E, const double *E0, const double *E1, const double *E2)
{
	for (int k = 0; k < NR; k++)
	for (int j = 0; j < NQ; j++)
	for (int i = 0; i < NP; i++) {
		const int l0 = j + NQ*(k + NR*(i-1));
		const int l1 = k + NR*(i + NP*(j-1));
		const int l2 = i + NP*(j + NQ*(k-1));

		const double e0 = interpolate_line<NP, NQ*NR>(E0, i, l0);
		const double e1 = interpolate_line<NQ, NR*NQ>(E1, j, l1);
		const double e2 = interpolate_line<NR, NP*NP>(E2, k, l2);

		const int l = i + NP*(j + NQ*k);

		const double di = (double)i;
		const double dj = (double)j;
		const double dk = (double)k;

		E[l] = (float)(
				unit_covariant_basic_vector<c, 0, 0>(di, dj, dk) * e0 +
				unit_covariant_basic_vector<c, 1, 0>(di, dj, dk) * e1 +
				unit_covariant_basic_vector<c, 2, 0>(di, dj, dk) * e2);
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
		return 0.5*(a + b);
	}
}

template <int c>
void get_othnomal_component_magnetic(
		float *B, const double *B0, const double *B1, const double *B2)
{
	for (int k = 0; k < NR; k++)
	for (int j = 0; j < NQ; j++)
	for (int i = 0; i < NP; i++) {
		const int l0 = j-1 + (NQ-1)*(k-1 + (NR-1)*i);
		const int l1 = k-1 + (NR-1)*(i-1 + (NP-1)*j);
		const int l2 = i-1 + (NP-1)*(j-1 + (NQ-1)*k);

		const double b0 = interpolate_plane<NQ, NR, 1, NQ-1>(B0, j, k, l0);
		const double b1 = interpolate_plane<NR, NP, 1, NR-1>(B1, k, i, l1);
		const double b2 = interpolate_plane<NP, NQ, 1, NP-1>(B2, i, j, l2);

		const int l = i + NP*(j + NQ*k);

		const double di = (double)i;
		const double dj = (double)j;
		const double dk = (double)k;

		B[l] = (float)(
				unit_contravariant_basic_vector<c, 0, 0>(di, dj, dk) * b0 +
				unit_contravariant_basic_vector<c, 1, 0>(di, dj, dk) * b1 +
				unit_contravariant_basic_vector<c, 2, 0>(di, dj, dk) * b2);
	}
}

int output(const double *Ep, const double *Eq, const double *Er,
           const double *BP, const double *BQ, const double *BR,
           int nt)
{
	// The order of values is as (0,0,0)(1,0,0)(0,1,0)(1,1,0)(0,0,1)...
	float *EB = (float *)malloc(NP*NQ*NR*sizeof(float));
	if (EB == NULL) {
		fprintf(stderr, "Failed to allocate memory. %s:%d\n",
				__FILE__, __LINE__);
		return ENOMEM;
	}

	char fpath[1024];
	int fd;
	int _errno;
	int asize[] = {NP, NQ, NR};
	float t;

	// Output electric field

	if (sprintf(fpath, "%s/E%08d.bin", DIR_OUT, nt) == EOF) {
		fprintf(stderr, "Path %s/E%08d.bin is too long. %s:%d\n",
				DIR_OUT, nt, __FILE__, __LINE__);
		free(EB);
		return ENAMETOOLONG;
	}

	if ((fd = open(fpath, O_WRONLY | O_CREAT,
			S_IFREG | S_IRUSR | S_IRGRP | S_IROTH)) == -1)
		goto open_err;

	if (write(fd, asize, 3*sizeof(int)) == -1)
		goto write_err;

	t = ((float)nt*DT);
	if (write(fd, &t, sizeof(float)) == -1)
		goto write_err;

	// for x
	get_othnomal_component_electric<0>(EB, Ep, Eq, Er);
	if (write(fd, EB, sizeof(float)*NP*NQ*NR) == -1)
		goto write_err;

	// for y
	get_othnomal_component_electric<1>(EB, Ep, Eq, Er);
	if (write(fd, EB, sizeof(float)*NP*NQ*NR) == -1)
		goto write_err;

	// for z
	get_othnomal_component_electric<2>(EB, Ep, Eq, Er);
	if (write(fd, EB, sizeof(float)*NP*NQ*NR) == -1)
		goto write_err;

	close(fd);

	// Output magnetic flux densities

	if (sprintf(fpath, "%s/B%08d.bin", DIR_OUT, nt) == EOF) {
		fprintf(stderr, "Path %s/B%08d.bin is too long. %s:%d\n",
				DIR_OUT, nt, __FILE__, __LINE__);
		free(EB);
		return ENAMETOOLONG;
	}

	if ((fd = open(fpath, O_WRONLY | O_CREAT,
			S_IFREG | S_IRUSR | S_IRGRP | S_IROTH)) == -1)
		goto open_err;

	if (write(fd, asize, 3*sizeof(int)) == -1)
		goto write_err;

	t = ((float)nt-0.5)*DT;
	if (write(fd, &t, sizeof(float)) == -1)
		goto write_err;

	// for x
	get_othnomal_component_magnetic<0>(EB, BP, BQ, BR);
	if (write(fd, EB, sizeof(float)*NP*NQ*NR) == -1)
		goto write_err;

	// for y
	get_othnomal_component_magnetic<1>(EB, BP, BQ, BR);
	if (write(fd, EB, sizeof(float)*NP*NQ*NR) == -1)
		goto write_err;

	// for z
	get_othnomal_component_magnetic<2>(EB, BP, BQ, BR);
	if (write(fd, EB, sizeof(float)*NP*NQ*NR) == -1)
		goto write_err;

	close(fd);
	free(EB);
	return 0;

open_err:
	_errno = errno;
	fprintf(stderr, "Failed to open a file: %s. %s:%d\n",
					fpath, __FILE__, __LINE__);
	free(EB);
	return _errno;

write_err:
	_errno = errno;
	fprintf(stderr, "Failed to write to a file: %s. %s:%d\n",
			fpath, __FILE__, __LINE__);
	close(fd);
	free(EB);
	return _errno;
}
