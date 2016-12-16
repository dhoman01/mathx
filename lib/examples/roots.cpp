#include "roots.hpp"
#include <cstdlib>
#include <iostream>

int main(){
  /* Function to use is (x-1)^2 - 3 */
  auto f = [](double x){ return std::pow(x-1,2) - 3; };
  auto df = [](double x){ return 2 * (x - 1); };
  auto g = [](double x){  return -1 * (2 - std::pow(x,2)) / 2; };

  /* Set maximum iterations */
  int maxiter = 10000;

  /* Set tolerance */
  double tol = std::pow(10,-8);

  /* Find root using bisect */
  double bisect = mathx::roots::bisect(f, 0, 4, f(0), f(4), tol, maxiter);

  /* Find root using fixed point iteration */
  double fixed_point_iter = mathx::roots::fixed_point_iter(g, f, 2, tol, maxiter);

  /* Find root using Newton's method */
  double newtons_method = mathx::roots::newtons_method(f, df, -1, tol, maxiter);

  /* Find root using the Secant method */
  double secant_method = mathx::roots::secant_method(f, -1, 0, tol, maxiter);

  /* Find root using the globalization of the Secant method */
  double hybrid_method = mathx::roots::hybrid_method(f, -4, 4, tol, maxiter);

  /* Print results */
  std::cout << "bisect           = " << bisect << std::endl;
  std::cout << "fixed_point_iter = " << fixed_point_iter << std::endl;
  std::cout << "newtons_method   = " << newtons_method << std::endl;
  std::cout << "secant_method    = " << secant_method << std::endl;
  std::cout << "hybrid_method    = " << hybrid_method << std::endl;

  /***** OUTPUT *****/
  // bisect           = 2.73205
  // fixed_point_iter = -0.732051
  // newtons_method   = -0.732051
  // secant_method    = -0.732051
  // hybrid_method    = -0.732051

  return EXIT_SUCCESS;
}
