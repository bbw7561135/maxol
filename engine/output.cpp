#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>

#include "coordinate.h"

extern double dt;

struct vector othnomal_electric_field(int i, int j, int k,
		const double *Ep, const double *Eq, const double *Er);
struct vector othnomal_magnetic_flux(int i, int j, int k,
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

	if (EBX == NULL || EBY == NULL || EBZ == NULL) {
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

	for (int k; k < NR; k++)
	for (int j; j < NQ; j++)
	for (int i; i < NP; i++) {
		int l = i + NP*(j + NQ*k);
		struct vector e = othnomal_electric_field(i, j, k, Ep, Eq, Er);
		EBX[l] = (float)e._0;
		EBY[l] = (float)e._1;
		EBZ[l] = (float)e._2;
	}

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

	for (int k; k < NR; k++)
	for (int j; j < NQ; j++)
	for (int i; i < NP; i++) {
		int l = i + NP*(j + NQ*k);
		struct vector b = othnomal_magnetic_flux(i, j, k, BP, BQ, BR);
		EBX[l] = (float)b._0;
		EBY[l] = (float)b._1;
		EBZ[l] = (float)b._2;
	}

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

	if (X == NULL || Y == NULL || Z == NULL) {
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
		const int l = p + NP*(q + NQ*r);
		struct vector position =
				othnormal_position<0>({(double)p, (double)q, (double)r});
		X[l] = (float)position._0;
		Y[l] = (float)position._1;
		Z[l] = (float)position._2;
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
