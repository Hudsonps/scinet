This code solves the equations below:

  //k1*(x^2)*y - k2*(z^2) = 0
  //2x - y - c1 = 0
  //x + z = c2

  // [2k1, -c1 k1 + k2, 2c2 k2, -k2] * [x3, x2, x1, 1]^T = 0


It does so by using the gsl routines. Two major functions were implemented. The function solver solves the first equation. The function polynomial_solver solves the second equation. These functions use gsl under the hood.

