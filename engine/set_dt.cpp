#include <float.h>

#include "coordinate.h"
#include "../config/comp_param.h"
#include "../config/phys_param.h"

extern double dt;

void set_dt(double courant) {
	double dx_min = DBL_MAX;

	// Get the smallest grid size
	for (int r = 0; r < NR; r++)
	for (int q = 0; q < NQ; q++)
	for (int p = 0; p < NP; p++) {
		// TODO: Is this right?
		const struct vector pos = {(double)p, (double)q, (double)r};
		const double dd0 = covariant_metric_tensor<0, 0, 0>(pos);
		const double dd1 = covariant_metric_tensor<1, 1, 0>(pos);
		const double dd2 = covariant_metric_tensor<2, 2, 0>(pos);
		const double dx = 1.0/sqrt(1.0/dd0 + 1.0/dd1 + 1.0/dd2);
		dx_min = dx < dx_min ? dx : dx_min;
	}

	// celeritas of light
	const double c = 1.0/sqrt(magnetic_permeability*electric_permittivity);
	dt = courant*dx_min/c;
}
