#ifndef ASSIGN5FUNCTIONS_H
#define ASSIGN5FUNCTIONS_H
#include <gsl/gsl_complex.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_multiroots.h>

//struct chemical_params {double x0; double y0; double z0; double k1; double k2;};
int chemical(gsl_vector *q, void *p, gsl_vector * f);
void print_state (size_t iter, gsl_multiroot_fsolver * s);
double solver(double x0, double y0, double z0, double k1, double k2);
void polynomial_solver(double x_0, double y_0, double z_0, double k1, double k2);

#endif
