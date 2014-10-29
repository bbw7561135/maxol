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

template <int d, int N0, int N1, int N2>
void get_othnomal_component_electric(
		float *E, const double *E0, const double *E1, const double *E2)
{
	for (int k = 0; k < N2; k++)
	for (int j = 0; j < N1; j++)
	for (int i = 0; i < N0; i++) {
		const int l0 = j   + N1*(k + N2*(i-1));
		const int l1 = k   + N2*(i + N0*(j-1));
		const int l2 = i   + N0*(j + N1*(k-1));

		// Linear-extrapolate if on boundary else linear-interpolate.
		double e0, e1, e2;
		switch (i) {
		case 0:    e0 =   1.5*E0[l0+N1*N2] - 0.5*E0[l0+2*N1*N2]; break;
		case N0-1: e0 = - 0.5*E0[l0-N1*N2] + 1.5*E0[l0];         break;
		default:   e0 = 0.5 * (E0[l0] + E0[l0+N1*N2]);
		}
		switch (j) {
		case 0:    e1 =   1.5*E1[l1+N2*N0] - 0.5*E1[l1+2*N2*N0]; break;
		case N1-1: e1 = - 0.5*E1[l1-N2*N0] + 1.5*E1[l1];         break;
		default:   e1 = 0.5 * (E1[l1] + E1[l1+N2*N0]);
		}
		switch (k) {
		case 0:    e2 =   1.5*E2[l2+N0*N1] - 0.5*E1[l2+2*N0*N1]; break;
		case N2-1: e2 = - 0.5*E2[l2-N0*N1] + 1.5*E2[l2];         break;
		default:   e2 = 0.5 * (E2[l2] + E2[l2+N0*N1]);
		}

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

		// Linear-extrapolate if on boundary else linear-interpolate.
		double b0, b1, b2;
		switch (i) {
		case 0:
			b0 = 1.5*B0[l0+(N1-1)*(N2-1)] - 0.5*B0[l0+2*(N1-1)*(N2-1)];
			break;
		case N0-1:
			b0 = - 0.5*B0[l0-(N1-1)*(N2-1)] + 1.5*B0[l0];
			break;
		default:
			b0 = 0.5*(B0[l0] + B0[l0+(N1-1)*(N2-1)]);
		}
		switch (j) {
		case 0:
			b1 = 1.5*B1[l1+(N2-1)*(N0-1)] - 0.5*B1[l1+2*(N2-1)*(N0-1)];
			break;
		case N1-1:
			b1 = - 0.5*B1[l1-(N2-1)*(N0-1)] + 1.5*B1[l1];
			break;
		default:
			b1 = 0.5 * (B1[l1] + B1[l1+(N2-1)*(N0-1)]);
		}
		switch (k) {
		case 0:
			b2 = 1.5*B2[l2+(N0-1)*(N1-1)] - 0.5*B1[l2+2*(N0-1)*(N1-1)];
			break;
		case N2-1:
			b2 = - 0.5*B2[l2-(N0-1)*(N1-1)] + 1.5*B2[l2];
			break;
		default:
			b2 = 0.5 * (B2[l2] + B2[l2+(N0-1)*(N1-1)]);
		}

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
	double t;

	// Output electric field

	if (sprintf(fpath, "%s/E%08d.bin", DIR_OUT, nt) == EOF) {
		fprintf(stderr, "Path %s/E%08d.bin is too long. %s:%d\n",
				DIR_OUT, nt, __FILE__, __LINE__);
		free(EB);
		return ENAMETOOLONG;
	}

	fd = open(fpath, O_WRONLY | O_CREAT,
			S_IFREG | S_IRUSR | S_IRGRP | S_IROTH);
	if (fd == -1)
		goto open_err;

	t = ((double)nt*DT);
	if (write(fd, &t, sizeof(double)) == -1)
		goto write_err;

	// for x
	get_othnomal_component_electric<0, NP, NQ, NR>(EB, Ep, Eq, Er);
	if (write(fd, &EB, sizeof(float)*NP*NQ*NR) == -1)
		goto write_err;

	// for y
	get_othnomal_component_electric<1, NQ, NR, NP>(EB, Eq, Er, Ep);
	if (write(fd, &EB, sizeof(float)*NP*NQ*NR) == -1)
		goto write_err;

	// for z
	get_othnomal_component_electric<2, NR, NP, NQ>(EB, Er, Ep, Eq);
	if (write(fd, &EB, sizeof(float)*NP*NQ*NR) == -1)
		goto write_err;

	close(fd);

	// Output magnetic flux densities

	if (sprintf(fpath, "%s/B%08d.bin", DIR_OUT, nt) == EOF) {
		fprintf(stderr, "Path %s/B%08d.bin is too long. %s:%d\n",
				DIR_OUT, nt, __FILE__, __LINE__);
		free(EB);
		return ENAMETOOLONG;
	}

	fd = open(fpath, O_WRONLY | O_CREAT,
			S_IFREG | S_IRUSR | S_IRGRP | S_IROTH);
	if (fd == -1)
		goto open_err;

	t = ((double)nt-0.5)*DT;
	if (write(fd, &t, sizeof(double)) == -1)
		goto write_err;

	// for x
	get_othnomal_component_magnetic<0, NP, NQ, NR>(EB, BP, BQ, BR);
	if (write(fd, &EB, sizeof(float)*NP*NQ*NR) == -1)
		goto write_err;

	// for y
	get_othnomal_component_magnetic<1, NQ, NR, NP>(EB, BQ, BR, BP);
	if (write(fd, &EB, sizeof(float)*NP*NQ*NR) == -1)
		goto write_err;

	// for z
	get_othnomal_component_magnetic<2, NR, NP, NQ>(EB, BR, BP, BQ);
	if (write(fd, &EB, sizeof(float)*NP*NQ*NR) == -1)
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
