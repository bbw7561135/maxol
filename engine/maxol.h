#ifndef MAXOL_H_
#define MAXOL_H_

#include "vector.h"
#include "../config/comp_param.h"
#include "../config/phys_param.h"

void set_dt(double courant);

void set_init_cond(double *Ep, double *Eq, double *Er,
                   double *BP, double *BQ, double *BR);

double evolute_elect(
		// electric fields at nt
		double *Ep, double *Eq, double *Er,
		 // magnetic flux densities at nt + 0.5
		const double *BP, const double *BQ, const double *BR);
double evolute_magnt(
		// magnetic flux densities at nt
		double *BP, double *BQ, double *BR,
		// electric fields at nt + 0.5
		const double *Ep, const double *Eq, const double *Er);

void bound_cond_elect(double *Ep, double *Eq, double *Er);
void bound_cond_magnt(double *BP, double *BQ, double *BR);

int output(const double *Ep, const double *Eq, const double *Er,
           const double *BP, const double *BQ, const double *BR,
           int nt);

struct vector total_momuntum(
		// electric field vectors
		const double *Ep, const double *Eq, const double *Er,
		// magnetic flux density vectors
		const double *BP, const double *BQ, const double *BR);

double div_rms_elect(const double *Ep, const double *Eq, const double *Er);
double div_rms_magnt(const double *BP, const double *BQ, const double *BR);

#endif /* MAXOL_H_ */
