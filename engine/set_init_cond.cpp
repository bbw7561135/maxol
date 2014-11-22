#include <math.h>
#include <stdio.h>

#include "coordinate.h"
#include "variable.h"
#include "../config/comp_param.h"
#include "../config/my_init_cond.h"

template <int c, int N0, int N1, int N2>
static void set_init_cond_elect(double *E)
{
	for (int i = 0; i < N0-1; i++)
	for (int k = 0; k < N2;   k++)
	for (int j = 0; j < N1;   j++) {
		const struct vector pos = {(double)i+0.5, (double)j, (double)k};

		// othnormal component
		const struct vector p_xyz = othnormal_position<c>(pos);
		const struct vector e_xyz = init_cond_elect(p_xyz);

		// unphysical components along to the grid.
		// Contravariant basic vector is a member of coordinate transformation
		// matrix here.
		const double e = inner_product(
				contravariant_basic_vector<c, c>(pos), e_xyz);

		// Convert to physical length.
		*elect_field<N1, N2>(i, j, k, E) =
				e * vector_length(covariant_basic_vector<c, c>(pos));
	}
}

template <int c, int N0, int N1, int N2>
static void set_init_cond_magnt(double *B)
{
	for (int i = 0; i < N0;   i++)
	for (int k = 0; k < N2-1; k++)
	for (int j = 0; j < N1-1; j++) {
		const struct vector pos = {(double)i, (double)j+0.5, (double)k+0.5};

		// othnormal component
		const struct vector p_xyz = othnormal_position<c>(pos);
		const struct vector b_xyz = init_cond_magnt(p_xyz);

		// physical component vertical with the grid
		// Covariant basic vector is a member of coordinate transformation
		// matrix here.
		const double b = inner_product(
				covariant_basic_vector<c, c>(pos), b_xyz);

		// Convert to physical length.
		*magnt_flux<N1, N2>(i, j, k, B) =
				b * vector_length(contravariant_basic_vector<c, c>(pos));
	}
}

void set_init_cond(double *Ep, double *Eq, double *Er,
                   double *BP, double *BQ, double *BR)
{
	set_init_cond_elect<0, NP, NQ, NR>(Ep);
	set_init_cond_elect<1, NQ, NR, NP>(Eq);
	set_init_cond_elect<2, NR, NP, NQ>(Er);

	set_init_cond_magnt<0, NP, NQ, NR>(BP);
	set_init_cond_magnt<1, NQ, NR, NP>(BQ);
	set_init_cond_magnt<2, NR, NP, NQ>(BR);
}
