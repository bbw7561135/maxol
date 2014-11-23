#include "coordinate.h"
#include "vector.h"
#include "../config/comp_param.h"
#include "../config/phys_param.h"

struct vector othnomal_electric_field(int i, int j, int k,
		const double *Ep, const double *Eq, const double *Er);
struct vector othnomal_magnetic_flux(int i, int j, int k,
		const double *BP, const double *BQ, const double *BR);

struct vector poynting_vector(
		int i, int j, int k,
		// electric field vectors
		const double *Ep, const double *Eq, const double *Er,
		// magnetic flux density vectors
		const double *BP, const double *BQ, const double *BR
		)
{
	struct vector e = othnomal_electric_field(i, j, k, Ep, Eq, Er);
	struct vector b = othnomal_magnetic_flux(i, j, k, BP, BQ, BR);
	return exterior_product(e, b);
}

struct vector total_momuntum(
		// electric field vectors
		const double *Ep, const double *Eq, const double *Er,
		// magnetic flux density vectors
		const double *BP, const double *BQ, const double *BR)
{
	struct vector mom = {0.0, 0.0, 0.0};

	for (int k = 0; k < NR; k++)
	for (int j = 0; j < NQ; j++)
	for (int i = 0; i < NP; i++) {
		struct vector poy = poynting_vector(i, j, k, Ep, Eq, Er, BP, BQ, BR);
		double jac = jacobian<0>({(double)i, (double)j, (double)k});
		mom = vector_add(mom, vector_multiple(poy, jac));
	}

	return vector_multiple(mom, magnetic_permeability);
}
