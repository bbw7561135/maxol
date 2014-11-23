/*
 * Maxol - A Maxwell's equations solver with FDTD method
 *         in a curvilinear coordinates system
 *
 * Author: Kazuho FUJII (kazuho.fujii@gmail.com)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <errno.h>

#include "maxol.h"

static int cpu_sec(clock_t start_clk) {
	return (int)((clock()-start_clk)/(CLOCKS_PER_SEC));
}

double dt;

int main(int argc, char **argv)
{
	if (argc < 3) {
		fprintf(stderr, "Usage: maxol NT_END NT_OUT [COURANT]\n");
		return 1;
	}

	int nt_end = strtol(argv[1], NULL, 10);
	int nt_out = strtol(argv[2], NULL, 10);

	if (nt_end == 0 || nt_out == 0) {
		fprintf(stderr, "The 1st & 2nd arguments must be non-zero integers\n");
		return 1;
	}

	// determine and set time skip
	double courant;
	const double courant_def = 0.8;
	if (argc >= 4) {
		courant = strtod(argv[3], NULL);
		if (courant == 0.0 || courant > 1.0) {
			fprintf(stderr,
					"Failed to convert the 3rd argument to double precision. Use default courant number %f\n",
					courant_def);
			courant = courant_def;
		}
	} else {
		courant = courant_def;
	}

	set_dt(courant);
	fprintf(stderr, "Set dt = %E (Courant factor = %f)\n", dt, courant);

	fprintf(stdout, "### Maxol ###\n");

	/*
	 * In this code, I use modified Yee's lattice system.
	 * I modified Yee's lattice for curvilinear coordinates system.
	 *
	 * Electric field vectors (Ep, Eq, Er) are contravariant physical
	 * components; Magnetic flux density vectors (BP, BQ, BR) are covariant
	 * physical components in curvilinear coordinates.
	 *
	 *              ^ Er     ^ BR
	 *          +---|-------+|
	 *         /    |      /||
	 *       //     *     / ||
	 *      //           /  |
	 *     /+-----------+  ======>
	 * BP L |           |   +   Eq
	 *      |           |  /
	 *      |    /      | /
	 *      |   /       |/
	 *      +--/--=======>
	 *        L Ep      BQ
	 *
	 * The covariant basic vectors are along to the computational grids.
	 * The contravariant basic vectors are vertical to the computational grids.
	 *
	 * Grid position are like this.
	 *
	 * (p) 0   1  ... NP-1
	 *                      (q)
	 *   +---+---+---+---+
	 *   | +-----------+ |   0
	 *   +-|-+---+---+-|-+
	 *   | | |   |   | | |   1
	 *   +-|-+---+---+-|-+
	 *   | | |   |   | | |   ...
	 *   +-|-+---+---+-|-+
	 *   | +-----------+ |   NQ-1
	 *   +---+-|-+---+---+
	 *         |
	 *         Boundary
	 *
	 * Position 0 and NP-1 (or NQ-1, NR-1) mean boundary.
	 *
	 * Ep starts from p = 1/2.
	 * Eq starts from q = 1/2.
	 * Er starts from r = 1/2.
	 * BP starts from q = 1/2, r = 1/2.
	 * BQ starts from r = 1/2, q = 1/2.
	 * BR starts from p = 1/2, q = 1/2.
	 *
	 * In arrays, the order of values is different with component.
	 * In Ep and BP, ordered as(0,0,0) (0,1,0)...(0,0,1)(0,1,1)...(1,1,1)...
	 * In Eq and BQ, ordered as(0,0,0) (0,0,1)...(1,0,0)(1,0,1)...(1,1,1)...
	 * In Er and BR, ordered as(0,0,0) (1,0,0)...(0,1,0)(1,1,0)...(1,1,1)...
	 *
	 * This is for code simplification and caching optimization.
	 */

	// Contravariant physical components of electric field.
	double *Ep = (double *)malloc(NQ*NR*(NP-1)*sizeof(double));
	double *Eq = (double *)malloc(NR*NP*(NQ-1)*sizeof(double));
	double *Er = (double *)malloc(NP*NQ*(NR-1)*sizeof(double));

	// Covariant physical components of magnetic flux density.
	double *BP = (double *)malloc((NQ-1)*(NR-1)*NP*sizeof(double));
	double *BQ = (double *)malloc((NR-1)*(NP-1)*NQ*sizeof(double));
	double *BR = (double *)malloc((NP-1)*(NQ-1)*NR*sizeof(double));

	if (Ep == NULL || Eq == NULL || Er == NULL ||
		BP == NULL || BQ == NULL || BR == NULL) {
		fprintf(stderr, "Failed to allocate memory. %s:%d\n",
				__FILE__, __LINE__);
		return 1;
	}

	set_init_cond(Ep, Eq, Er, BP, BQ, BR);

	int err = output(Ep, Eq, Er, BP, BQ, BR, /*nt=*/0);
	if (err) {
		fprintf(stderr, "Failed to save a data. %s:%d (%s)\n",
				__FILE__, __LINE__, strerror(err));
		return 2;
	}

	clock_t start_clk = clock(); // ignore err;

	fprintf(stdout, "cpu_time	nt	time	electric_energy	magnetic_energy");
	fprintf(stdout, "	momemtum_x	momemtum_y	momentum_z");
	fprintf(stdout, "	div_electric	div_magnetic");
	fprintf(stdout, "\n");
	fprintf(stdout, "%d	%d	%E	%E	%E",
			cpu_sec(start_clk), 0, 0.0, 0.0, 0.0);
	struct vector mom = total_momuntum(Ep, Eq, Er, BP, BQ, BR);
	fprintf(stdout, "	%E	%E	%E", mom._0, mom._1, mom._2);
	double div_elect = div_rms_elect(Ep, Eq, Er);
	double div_magnt = div_rms_magnt(BP, BQ, BR);
	fprintf(stdout, "	%E	%E", div_elect, div_magnt);
	fprintf(stdout, "\n");

	// Main loop
	for (int nt = 1; nt <= nt_end; nt++) {
		/*
		 * For time evolution, I use leap-frog method.
		 * Thus, B is at step of nt - 0.5
		 */
		double eng_magnt = evolute_magnt(BP, BQ, BR, Ep, Eq, Er);

		// Boundary condition of magnet flux is not required.
		bound_cond_magnt(BP, BQ, BR);

		double eng_elect = evolute_elect(Ep, Eq, Er, BP, BQ, BR);

		bound_cond_elect(Ep, Eq, Er);

		// Save data
		if (nt % nt_out == 0) {
			int err = output(Ep, Eq, Er, BP, BQ, BR, nt);
			if (err) {
				fprintf(stderr, "Failed to save a data. %s:%d (%s)\n",
						__FILE__, __LINE__, strerror(err));
				return 2;
			}
			fprintf(stderr, "Saved to file. nt=%d\n", nt);
		}

		fprintf(stdout, "%d	%d	%E	%E	%E",
				cpu_sec(start_clk), nt, (double)nt*dt, eng_elect, eng_magnt);
		mom = total_momuntum(Ep, Eq, Er, BP, BQ, BR);
		fprintf(stdout, "	%E	%E	%E", mom._0, mom._1, mom._2);

		div_elect = div_rms_elect(Ep, Eq, Er);
		div_magnt = div_rms_magnt(BP, BQ, BR);
		fprintf(stdout, "	%E	%E", div_elect, div_magnt);
		fprintf(stdout, "\n");
	}

	return 0;
}
