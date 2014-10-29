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

template <int d, int N0, int N1, int N2>
void get_othnomal_component_electric(
		float *E, const double *E0, const double *E1, const double *E2)
{
	for (int k = 0; k < N2; k++)
	for (int j = 0; j < N1; j++)
	for (int i = 0; i < N0; i++) {
		const int l0 = j + N1*(k + N2*(i-1));
		const int l1 = k + N2*(i + N0*(j-1));
		const int l2 = i + N0*(j + N1*(k-1));

		const double e0 = interpolate_line<N0, N1*N2>(E0, i, l0);
		const double e1 = interpolate_line<N1, N2*N0>(E1, j, l1);
		const double e2 = interpolate_line<N2, N0*N1>(E2, k, l2);

		int _i = i, _j = j, _k = k;
		swap<d>(&_i, &_j, &_k);
		const int _l = _i + NP*(_j + NQ*_k); // Here, NP,NQ,NR are OK.

		E[_l] = (float)(
				unit_covariant_basic_vector<d, 0, 0>(
					(double)_i, (double)_j, (double)_k) * e0 +
				unit_covariant_basic_vector<d, 1, 0>(
					(double)_i, (double)_j, (double)_k) * e1 +
				unit_covariant_basic_vector<d, 2, 0>(
					(double)_i, (double)_j, (double)_k) * e2);
	}
}

// Linear-extrapolate if on boundary else linear-interpolate.
template <int N1, int N2>
double interpolate_plane(const double *E, int i, int j, int l)
{
	// TODO
}

template <int d, int N0, int N1, int N2>
void get_othnomal_component_magnetic(
		float *B, const double *B0, const double *B1, const double *B2)
{
	for (int k = 0; k < N2; k++)
	for (int j = 0; j < N1; j++)
	for (int i = 0; i < N0; i++) {
		const int l0 = j-1 + (N1-1)*(k-1 + (N2-1)*i);
		const int l1 = k-1 + (N2-1)*(i-1 + (N0-1)*j);
		const int l2 = i-1 + (N0-1)*(j-1 + (N1-1)*k);

		const double b0 = interpolate_plane<N1, N2>(B0, i, l0);
		const double b1 = interpolate_plane<N1, N2>(B1, j, l1);
		const double b2 = interpolate_plane<N1, N2>(B2, k, l2);

		int _i = i, _j = j, _k = k;
		swap<d>(&_i, &_j, &_k);
		const int _l = _i + NP*(_j + NQ*_k); // Here, NP,NQ,NR are OK.

		B[_l] = (float)(
				unit_contravariant_basic_vector<d, 0, 0>(
					(double)_i, (double)_j, (double)_k) * b0 +
				unit_contravariant_basic_vector<d, 1, 0>(
					(double)_i, (double)_j, (double)_k) * b1 +
				unit_contravariant_basic_vector<d, 2, 0>(
					(double)_i, (double)_j, (double)_k) * b2);
	}
}

int output(const double *Ep, const double *Eq, const double *Er,
           const double *BP, const double *BQ, const double *BR,
           int nt)
{
	// The order of values is as (0,0,0)(1,0,0)(0,1,0)(1,1,0)(0,0,1)...
	float *EB = (float *)calloc(NP*NQ*NR, sizeof(float));
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
	get_othnomal_component_electric<0, NP, NQ, NR>(EB, Ep, Eq, Er);
	if (write(fd, EB, sizeof(float)*NP*NQ*NR) == -1)
		goto write_err;

	// for y
	get_othnomal_component_electric<1, NQ, NR, NP>(EB, Eq, Er, Ep);
	if (write(fd, EB, sizeof(float)*NP*NQ*NR) == -1)
		goto write_err;

	// for z
	get_othnomal_component_electric<2, NR, NP, NQ>(EB, Er, Ep, Eq);
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
	get_othnomal_component_magnetic<0, NP, NQ, NR>(EB, BP, BQ, BR);
	if (write(fd, EB, sizeof(float)*NP*NQ*NR) == -1)
		goto write_err;

	// for y
	get_othnomal_component_magnetic<1, NQ, NR, NP>(EB, BQ, BR, BP);
	if (write(fd, EB, sizeof(float)*NP*NQ*NR) == -1)
		goto write_err;

	// for z
	get_othnomal_component_magnetic<2, NR, NP, NQ>(EB, BR, BP, BQ);
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
