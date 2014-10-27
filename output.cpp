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

int output(const double *Ep, const double *Eq, const double *Er,
           const double *BP, const double *BQ, const double *BR,
           int nt)
{
	double *EB = (double *) malloc(NP*NQ*NR*sizeof(double));
	if (EB == NULL) {
		fprintf(stderr, "Failed to allocate memory. %s:%d\n",
				__FILE__, __LINE__);
		return ENOMEM;
	}

	// Output magnetic flux densities

	char fpath[1024];
	if (sprintf(fpath, "%s/B%08d.bin", DIR_OUT, nt) == EOF) {
		fprintf(stderr, "Path %s/B%08d.bin is too long. %s:%d\n",
				DIR_OUT, nt, __FILE__, __LINE__);
		free(EB);
		return ENAMETOOLONG;
	}

	const int fdB = open(fpath, O_WRONLY | O_CREAT);
	if (fdB == -1) {
		fprintf(stderr, "Failed to open a file: %s. %s:%d\n",
				fpath, __FILE__, __LINE__);
		free(EB);
		return errno;
	}

	const double tB = ((double)nt-0.5)*DT;
	if (write(fdB, &tB, sizeof(double)) == -1) {
		fprintf(stderr, "Failed to write to a file: %s. %s:%d\n",
				fpath, __FILE__, __LINE__);
		close(fdB);
		free(EB);
		return errno;
	}

	// TODO

	close(fdB);

	// Output electric field

	if (sprintf(fpath, "%s/E%08d.bin", DIR_OUT, nt) == EOF) {
		fprintf(stderr, "Path %s/E%08d.bin is too long. %s:%d\n",
				DIR_OUT, nt, __FILE__, __LINE__);
		free(EB);
		return ENAMETOOLONG;
	}

	const int fdE = open(fpath, O_WRONLY | O_CREAT);
	if (fdE == -1) {
		fprintf(stderr, "Failed to open a file: %s. %s:%d\n",
						fpath, __FILE__, __LINE__);
		free(EB);
		return errno;
	}

	const double tE = ((double)nt*DT);
	if (write(fdE, &tE, sizeof(double)) == -1) {
		fprintf(stderr, "Failed to write to a file: %s. %s:%d\n",
				fpath, __FILE__, __LINE__);
		close(fdE);
		free(EB);
		return errno;
	}

	// TODO

	close(fdE);
	free(EB);

	return 0;
}
