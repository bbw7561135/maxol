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

extern double dt;

void get_othnomal_cmpo_elect(
		float *EX, float *EY, float *EZ,
		const double *EP, const double *EQ, const double *ER);

void get_othnomal_cmpo_magnt(
		float *BX, float *BY, float *BZ,
		const double *BP, const double *BQ, const double *BR);

int output_grid(int nt);

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
		char outpath_def[] = ".";
		outpath = outpath_def;
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

	t = ((float)nt*dt);

	if (write(fd, header, sizeof(header)) == -1) goto err;
	if (write(fd, asize,  3*sizeof(int))  == -1) goto err;
	if (write(fd, &nt,    sizeof(int))    == -1) goto err;
	if (write(fd, &t,     sizeof(float))  == -1) goto err;

	get_othnomal_cmpo_elect(EBX, EBY, EBZ, Ep, Eq, Er);

	if (write(fd, EBX, NP*NQ*NR*sizeof(float)) == -1) goto err;
	if (write(fd, EBY, NP*NQ*NR*sizeof(float)) == -1) goto err;
	if (write(fd, EBZ, NP*NQ*NR*sizeof(float)) == -1) goto err;

	close(fd);

	// Output magnetic flux densities

	if (sprintf(fpath, "%s/%08dB", outpath, nt) == EOF) {
		free(EBX); free(EBY); free(EBZ);
		return ENAMETOOLONG;
	}

	if ((fd = open(fpath, O_WRONLY | O_CREAT,
			S_IFREG | S_IRUSR | S_IRGRP | S_IROTH)) == -1)
		goto err;

	t = ((float)nt-0.5)*dt;

	if (write(fd, header, sizeof(header)) == -1) goto err;
	if (write(fd, asize,  3*sizeof(int))  == -1) goto err;
	if (write(fd, &nt,    sizeof(int))    == -1) goto err;
	if (write(fd, &t,     sizeof(float))  == -1) goto err;

	get_othnomal_cmpo_magnt(EBX, EBY, EBZ, BP, BQ, BR);

	if (write(fd, EBX, NP*NQ*NR*sizeof(float)) == -1) goto err;
	if (write(fd, EBY, NP*NQ*NR*sizeof(float)) == -1) goto err;
	if (write(fd, EBZ, NP*NQ*NR*sizeof(float)) == -1) goto err;

	close(fd);
	free(EBX); free(EBY); free(EBZ);

	// Output grid data only at initialization.
	if (nt == 0) return output_grid(nt);

	return 0;

err:
	_errno = errno;
	free(EBX); free(EBY); free(EBZ);
	close(fd);
	return _errno;
}

int output_grid(int nt) {
	// The order of values is as (0,0,0)(1,0,0)(0,1,0)(1,1,0)(0,0,1)...
	float *X = (float *)malloc(NP*NQ*NR*sizeof(float));
	float *Y = (float *)malloc(NP*NQ*NR*sizeof(float));
	float *Z = (float *)malloc(NP*NQ*NR*sizeof(float));

	if (X == NULL || Y == NULL || Z == NULL ) {
		free(X); free(Y); free(Z);
		return ENOMEM;
	}

	char *outpath = getenv("MAXOL_OUT_PATH");
	if (outpath == NULL) {
		char outpath_def[] = ".";
		outpath = outpath_def;
	}

	char fpath[1024];
	int fd;
	int _errno;
	char header[] = "### Maxol ###\n";
	int asize[] = {NP, NQ, NR};
	float t = (float)nt*dt;

	if (sprintf(fpath, "%s/%08dG", outpath, nt) == EOF) {
		free(X); free(Y); free(Z);
		return ENAMETOOLONG;
	}

	if ((fd = open(fpath, O_WRONLY | O_CREAT,
			S_IFREG | S_IRUSR | S_IRGRP | S_IROTH)) == -1)
		goto err;

	if (write(fd, header, sizeof(header)) == -1) goto err;
	if (write(fd, asize,  3*sizeof(int))  == -1) goto err;
	if (write(fd, &nt,    sizeof(int))    == -1) goto err;
	if (write(fd, &t,     sizeof(float))  == -1) goto err;

	for (int r = 0; r < NR; r++)
	for (int q = 0; q < NQ; q++)
	for (int p = 0; p < NP; p++) {
		const int l = p + NP*(q+NQ*r);
		X[l] = (float)othnormal_position<0, 0>(p, q, r);
		Y[l] = (float)othnormal_position<1, 0>(p, q, r);
		Z[l] = (float)othnormal_position<2, 0>(p, q, r);
	}

	if (write(fd, X, NP*NQ*NR*sizeof(float)) == -1) goto err;
	if (write(fd, Y, NP*NQ*NR*sizeof(float)) == -1) goto err;
	if (write(fd, Z, NP*NQ*NR*sizeof(float)) == -1) goto err;

	close(fd);
	free(X); free(Y); free(Z);

	return 0;

err:
	_errno = errno;
	free(X); free(Y); free(Z);
	close(fd);
	return _errno;
}
