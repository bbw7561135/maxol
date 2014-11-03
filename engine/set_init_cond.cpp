#include <math.h>
#include <stdio.h>

#include "coordinate.h"
#include "../config/comp_param.h"
#include "../config/my_init_cond.h"

template <int c, int N0, int N1, int N2>
static void set_init_cond_elect(double *E)
{
	for (int i = 0; i < N0; i++)
	for (int k = 0; k < N2; k++)
	for (int j = 0; j < N1; j++) {
		const double p0 = (double)i + 0.5;
		const double p1 = (double)j;
		const double p2 = (double)k;

		const double x = othnormal_position<0, c>(p0, p1, p2);
		const double y = othnormal_position<1, c>(p0, p1, p2);
		const double z = othnormal_position<2, c>(p0, p1, p2);

		// othnormal component
		const double ex = init_cond_elect_x(x, y, z);
		const double ey = init_cond_elect_y(x, y, z);
		const double ez = init_cond_elect_z(x, y, z);
		const int l = j + N1*(k + N2*i);

		// unphysical components along to the grid.
		// Contravariant basic vector is a member of coordinate transformation
		// matrix here.
		const double e = contravariant_basic_vector<0, c, c>(p0, p1, p2) * ex +
						 contravariant_basic_vector<1, c, c>(p0, p1, p2) * ey +
						 contravariant_basic_vector<2, c, c>(p0, p1, p2) * ez;

		// Convert to physical length.
		E[l] = e * len_of_covariant_basic_vector<c, c>(p0, p1, p2);
	}
}

template <int c, int N0, int N1, int N2>
static void set_init_cond_magnt(double *B)
{
	for (int i = 0; i < N0; i++)
	for (int k = 0; k < N2; k++)
	for (int j = 0; j < N1; j++) {
		const double p0 = (double)i;
		const double p1 = (double)j + 0.5;
		const double p2 = (double)k + 0.5;

		const double x = othnormal_position<0, c>(p0, p1, p2);
		const double y = othnormal_position<1, c>(p0, p1, p2);
		const double z = othnormal_position<2, c>(p0, p1, p2);

		// othnormal component
		const double bx = init_cond_magnt_x(x, y, z);
		const double by = init_cond_magnt_y(x, y, z);
		const double bz = init_cond_magnt_z(x, y, z);

		int l = i + N0*(j + N1*k);

		// physical component vertical with the grid
		// Covariant basic vector is a member of coordinate transformation
		// matrix here.
		const double b = covariant_basic_vector<0, c, c>(p0, p1, p2) * bx +
						 covariant_basic_vector<1, c, c>(p0, p1, p2) * by +
						 covariant_basic_vector<2, c, c>(p0, p1, p2) * bz;

		// Convert to physical length.
		B[l] = b * len_of_contravariant_basic_vector<c, c>(p0, p1, p2);
	}
}

void set_init_cond(double *Ep, double *Eq, double *Er,
                   double *BP, double *BQ, double *BR)
{
	set_init_cond_elect<0, NP-1, NQ, NR>(Ep);
	set_init_cond_elect<1, NQ-1, NR, NP>(Eq);
	set_init_cond_elect<2, NR-1, NP, NQ>(Er);

	set_init_cond_magnt<0, NP, NQ-1, NR-1>(BP);
	set_init_cond_magnt<1, NQ, NR-1, NP-1>(BQ);
	set_init_cond_magnt<2, NR, NP-1, NQ-1>(BR);
}
