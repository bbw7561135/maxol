/*
 * Convert a binary output file to text.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>

#include "../engine/coordinate.h"

int main(int argc, char **argv)
{
	if (argc < 2) {
		fprintf(stderr, "Arguments are few.\n");
		return 1;
	}

	char *path = argv[1];

	int fd = open(path, O_RDONLY);
	if (fd == -1) {
		fprintf(stderr, "Failed to open %s (%s)\n", path, strerror(errno));
		return 1;
	}

	int size[3];
	int np, nq, nr;
	float time;
	float *ax = NULL, *ay = NULL, *az = NULL;
	int _errno = 0;

	if (read(fd, size, 3*sizeof(int)) != (ssize_t)(3*sizeof(int))) {
		_errno = errno;
		goto read_err;
	}

	np = size[0];
	nq = size[1];
	nr = size[2];

	if (read(fd, &time, sizeof(float)) != (ssize_t)sizeof(float)) {
		_errno = errno;
		goto read_err;
	}

	fprintf(stdout, "#time= %E\n", time);

	ax = (float *)malloc(np*nq*nr*sizeof(float));
	ay = (float *)malloc(np*nq*nr*sizeof(float));
	az = (float *)malloc(np*nq*nr*sizeof(float));

	if (ax == NULL || ay == NULL || az == NULL) {
		free(ax); free(ay); free(az);
		close(fd);
		fprintf(stderr, "Failed to allocate memory\n");
		return 1;
	}

	if (read(fd, ax, np*nq*nr*sizeof(float)) !=
			(ssize_t)(np*nq*nr*sizeof(float))) {
		_errno = errno;
		goto read_err;
	}

	if (read(fd, ay, np*nq*nr*sizeof(float)) !=
			(ssize_t)(np*nq*nr*sizeof(float))) {
		_errno = errno;
		goto read_err;
	}

	if (read(fd, az, np*nq*nr*sizeof(float)) !=
			(ssize_t)(np*nq*nr*sizeof(float))) {
		_errno = errno;
		goto read_err;
	}

	fprintf(stdout,
			"p	q	r	x	y	z	x_cmpo	y_cmpo	z_cmpo\n");

	for (int k = 0; k < nr; k++)
	for (int j = 0; j < nq; j++)
	for (int i = 0; i < np; i++) {
		fprintf(stdout, "%d	%d	%d	%E	%E	%E	%E	%E	%E\n",
				i, j, k,
				othnormal_position<0, 0>(i, j, k),
				othnormal_position<1, 0>(i, j, k),
				othnormal_position<2, 0>(i, j, k),
				ax[i+np*(j+nq)],
				ay[i+np*(j+nq)],
				az[i+np*(j+nq)]);
	}

	close(fd);
	free(ax); free(ay), free(az);
	return 0;

read_err:
	close(fd);
	free(ax); free(ay), free(az);
	fprintf(stderr, "Failed to read from %s (%s)\n", path, strerror(_errno));
	return 1;
}
