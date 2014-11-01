/*
 * Maxol: A Maxwell's equations solver with FDTD method
 *        in a curvilinear coordinates system
 *
 * Author: Kazuho FUJII (kazuho.fujii@gmail.com)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "maxol.h"

static int cpu_sec(clock_t start_clk) {
	return (int)((clock()-start_clk)/(CLOCKS_PER_SEC));
}

int main(int argc, char **argv)
{
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
	 *       //     *     / |*
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
	 */

	// Contravariant physical components of electric field.
	double *Ep, *Eq, *Er;
	// Covariant physical components of magnetic flux density.
	double *BP, *BQ, *BR;

	/*
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
	Ep = (double *)malloc(NQ*NR*(NP-1)*sizeof(double));
	Eq = (double *)malloc(NR*NP*(NQ-1)*sizeof(double));
	Er = (double *)malloc(NP*NQ*(NR-1)*sizeof(double));

	BP = (double *)malloc((NQ-1)*(NR-1)*NP*sizeof(double));
	BQ = (double *)malloc((NR-1)*(NP-1)*NQ*sizeof(double));
	BR = (double *)malloc((NP-1)*(NQ-1)*NR*sizeof(double));

	if (Ep == NULL || Eq == NULL || Er == NULL ||
		BP == NULL || BQ == NULL || BR == NULL) {
		fprintf(stderr, "Failed to allocate memory. %s:%d\n",
				__FILE__, __LINE__);
		free(Ep); free(Eq); free(Er);
		free(BP); free(BQ); free(BR);
		return 1;
	}

	set_init_cond(Ep, Eq, Er, BP, BQ, BR);

	int err = output(Ep, Eq, Er, BP, BQ, BR, /*nt=*/0);
	if (err) {
		fprintf(stderr, "Failed to save a data. %s:%d (%s)\n",
				__FILE__, __LINE__, strerror(err));
		return 2;
	}
	fprintf(stderr, "Saved to file. nt=%d\n", 0);

	clock_t start_clk = clock(); // ignore err;

	fprintf(stdout, "cpu_time	nt	time\n");
	fprintf(stdout, "%d	%d	%f\n", cpu_sec(start_clk), 0, 0.0);

	// Main loop
	for (int nt = 1; nt <= NT_END; nt++) {
		/*
		 * For time evolution, I use leap-frog method.
		 * Thus, B is at step of nt - 0.5
		 */
		evolute_magnt(BP, BQ, BR, Ep, Eq, Er);

		// Boundary condition of magnet flux is not required.
		// bound_cond_magnt(BP, BQ, BR, ((double)nt-0.5)*DT);

		evolute_elect(Ep, Eq, Er, BP, BQ, BR);

		bound_cond_elect(Ep, Eq, Er, (double)nt*DT);

		// Save data
		if (nt%NT_OUT == 0) {
			int err = output(Ep, Eq, Er, BP, BQ, BR, nt);
			if (err) {
				fprintf(stderr, "Failed to save a data. %s:%d (%s)\n",
						__FILE__, __LINE__, strerror(err));
				return 2;
			}
			fprintf(stderr, "Saved to file. nt=%d\n", nt);
		}

		fprintf(stdout, "%d	%d	%f\n",
				cpu_sec(start_clk), nt, (double)nt*DT);
	}

	free(Ep); free(Eq); free(Er);
	free(BP); free(BQ); free(BR);

	return 0;
}
