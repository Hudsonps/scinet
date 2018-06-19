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


int main()
{

  //equations to be solved
  //k1*(x^2)*y - k2*(z^2) = 0
  //2x - y - c1 = 0
  //x + z = c2

  //cubic equation to be solved:
  // [2k1, -c1 k1 + k2, 2c2 k2, -k2] * [x3, x2, x1, 1]^T = 0
  //The gsl routine assumes that the coefficient of x3 is one. This can be arranged by dividing the whole equation by k1

  double x_0 = 0.5;
  double y_0 = 1;
  double z_0 = 0;
  double k1 = 1;
  double k2 = 0.7;

  polynomial_solver(x_0, y_0, z_0, k1, k2);

  solver(x_0, y_0, z_0, k1, k2);

  

}
