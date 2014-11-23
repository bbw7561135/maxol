#include <math.h>

#include "../config/comp_param.h"
#include "coordinate.h"
#include "variable.h"
#include "vector.h"

double poisson_elect(double *Ep, double *Eq, double *Er) {
	return 0.0;
}

double div_rms_elect(const double *Ep, const double *Eq, const double *Er)
{
	double sq_div = 0.0;

	for (int k = 1; k < NR-1; k++)
	for (int j = 1; j < NQ-1; j++)
	for (int i = 1; i < NP-1; i++) {
		double p = (double)i;
		double q = (double)j;
		double r = (double)k;

		double div_pu = (*elect_field<0, NP, NQ, NR>(i-1, j, k, Ep)) *
					jacobian<0>({p-0.5, q, r}) /
					vector_length(covariant_basic_vector<0, 0>({p-0.5, q, r}));

		double div_pd = (*elect_field<0, NP, NQ, NR>(i, j, k, Ep)) *
					jacobian<0>({p+0.5, q, r}) /
					vector_length(covariant_basic_vector<0, 0>({p+0.5, q, r}));

		double div_qu = (*elect_field<1, NP, NQ, NR>(i, j-1, k, Eq)) *
					jacobian<0>({p, q-0.5, r}) /
					vector_length(covariant_basic_vector<1, 0>({p, q-0.5, r}));

		double div_qd = (*elect_field<1, NP, NQ, NR>(i, j, k, Eq)) *
					jacobian<0>({p, q+0.5, r}) /
					vector_length(covariant_basic_vector<1, 0>({p, q+0.5, r}));

		double div_ru = (*elect_field<2, NP, NQ, NR>(i-1, j, k, Er)) *
					jacobian<0>({p, q, r-0.5}) /
					vector_length(covariant_basic_vector<2, 0>({p, q, r-0.5}));

		double div_rd = (*elect_field<2, NP, NQ, NR>(i, j, k, Er)) *
					jacobian<0>({p, q, r+0.5}) /
					vector_length(covariant_basic_vector<2, 0>({p, q, r+0.5}));

		double div = (- div_pu + div_pd
					  - div_qu + div_qd
					  - div_ru + div_rd) / jacobian<0>({p, q, r});
		sq_div += div*div;
	}

	return sqrt(sq_div/((NP-2)*(NQ-2)*(NR-2)));
}

double div_rms_magnt(const double *BP, const double *BQ, const double *BR)
{
	double sq_div = 0.0;

	for (int k = 0; k < NR-1; k++)
	for (int j = 0; j < NQ-1; j++)
	for (int i = 0; i < NP-1; i++) {
		double p = (double)i + 0.5;
		double q = (double)j + 0.5;
		double r = (double)k + 0.5;

		double div_pu = (*magnt_flux<0, NP, NQ, NR>(i  , j, k, BP)) *
				jacobian<0>({p-0.5, q, r}) *
				vector_length(contravariant_basic_vector<0, 0>({p-0.5, q, r}));

		double div_pd = (*magnt_flux<0, NP, NQ, NR>(i+1, j, k, BP)) *
				jacobian<0>({p+0.5, q, r}) *
				vector_length(contravariant_basic_vector<0, 0>({p+0.5, q, r}));

		double div_qu = (*magnt_flux<1, NP, NQ, NR>(i, j  , k, BQ)) *
				jacobian<0>({p, q-0.5, r}) *
				vector_length(contravariant_basic_vector<1, 0>({p, q-0.5, r}));

		double div_qd = (*magnt_flux<1, NP, NQ, NR>(i, j+1, k, BQ)) *
				jacobian<0>({p, q+0.5, r}) *
				vector_length(contravariant_basic_vector<1, 0>({p, q+0.5, r}));

		double div_ru = (*magnt_flux<2, NP, NQ, NR>(i, j, k  , BR)) *
				jacobian<0>({p, q, r-0.5}) *
				vector_length(contravariant_basic_vector<2, 0>({p, q, r-0.5}));

		double div_rd = (*magnt_flux<2, NP, NQ, NR>(i, j, k+1, BR)) *
				jacobian<0>({p, q, r+0.5}) *
				vector_length(contravariant_basic_vector<2, 0>({p, q, r+0.5}));

		double div = (- div_pu + div_pd
					  - div_qu + div_qd
					  - div_ru + div_rd) / jacobian<0>({p, q, r});
		sq_div += div*div;
	}

	return sqrt(sq_div/((NP-1)*(NQ-1)*(NR-1)));
}
