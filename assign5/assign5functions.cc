#include <iostream>
#include <fstream>
#include <cmath>
#include <iostream>
#include "assign5functions.h"
#include <gsl/gsl_complex.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>

struct chemical_params{double x0; double y0; double z0; double k1; double k2;};

void polynomial_solver(double x_0, double y_0, double z_0, double k1, double k2)
{
  std::cout << "Starting the retrieval through the polynomial method" << std::endl;
  double c1 = 2*x_0 - y_0; //constants through the problem
  double c2 = x_0 + z_0; //constants through the problem			 
  double a = -(c1*k1 + k2)/(2*k1);
  double b = c2*k2/(k1);
  double c = -(k2*c2*c2)/(2*k1);
  
  double r0, r1, r2;

  gsl_poly_solve_cubic (a, b, c, &r0, &r1, &r2);
  // z3 +az2 +bz+c=0

  double x0answer = r0;
  double z0answer = c2 - r0;
  double y0answer = 2*r0 - c1;

  //double ftest = k1*x0answer*x0answer*y0answer - k2*z0answer*z0answer;

  std::cout << "The value of x is " << x0answer << std::endl;
  std::cout << "The value of y is " << y0answer << std::endl;
  std::cout << "The value of z is " << z0answer << std::endl;
}
  

void print_state (size_t iter, gsl_multiroot_fsolver * s)
{
  printf ("iter = %3u (x,y,z) = % .3f % .3f % .3f"
          " f(x,y,z) = % .3e % .3e % .3e\n",
          iter,
          gsl_vector_get (s->x, 0), 
          gsl_vector_get (s->x, 1),
	  gsl_vector_get (s->x, 2), 
          gsl_vector_get (s->f, 0),
          gsl_vector_get (s->f, 1), 
          gsl_vector_get (s->f, 2));
}

int chemical(const gsl_vector *q, void *p, gsl_vector *f) {
   struct chemical_params * params = (struct chemical_params *)p;
   const double x0 = (params-> x0);
   const double y0 = (params-> y0);
   const double z0 = (params-> z0);
   const double k1 = (params-> k1);
   const double k2 = (params-> k2);
   
   const double x = gsl_vector_get(q,0);
   const double y = gsl_vector_get(q,1);
   const double z = gsl_vector_get(q,2);

   
   const double c1 = 2*x0 - y0; //constants through the problem
   const double c2 = x0 + z0; //constants through the problem
   //equations to be solved
  //k1*(x^2)*y - k2*(z^2) = 0
  //2x - y - c1 = 0
  //x + z = c2
   
   gsl_vector_set(f, 0, k1*x*x*y - k2*z*z);
   gsl_vector_set (f, 1, 2*x - y - c1);
   gsl_vector_set (f, 2, x + z - c2);
   
   return GSL_SUCCESS;
}

double solver(double x0, double y0, double z0, double k1, double k2)
{
  std::cout << "Starting the retrieval through the multidimensional method" << std::endl;
  int status;
  size_t i, iter = 0;

  const size_t dim = 3;
  const int dimaux = 3;
  struct chemical_params p = {x0, y0, z0, k1, k2};
  gsl_multiroot_function f = {&chemical, dim, &p};

  gsl_vector *x = gsl_vector_alloc (dim);
  gsl_vector_set (x, 0, x0);
  gsl_vector_set (x, 1, y0);
  gsl_vector_set (x, 2, z0);

  const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrid;
  gsl_multiroot_fsolver *s  = gsl_multiroot_fsolver_alloc (T, dimaux);
  gsl_multiroot_fsolver_set (s, &f, x);

  do
    {
      iter++;
      status = gsl_multiroot_fsolver_iterate (s);

      print_state (iter, s);

      if (status)   /* check if solver is stuck */
        break;

      status = 
        gsl_multiroot_test_residual (s->f, 1e-7);
    }
  while (status == GSL_CONTINUE && iter < 1000);

  printf ("status = %s\n", gsl_strerror (status));

  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);
  return 0;
}


  




  
  
  
