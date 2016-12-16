#include "utils.hpp"
#include <cstdlib>
#include <iostream>

int main(){
  /* Calculate double machine precision */
  double eps_double = mathx::utils::precision::machine_epsilon<double>();

  /* Calculate single machine precision */
  float eps_float = mathx::utils::precision::machine_epsilon<float>();

  /* Print results */
  std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << "eps_double = " << eps_double << std::endl;
  std::cout << std::setprecision(std::numeric_limits<float>::digits10 + 1) << "eps_float = " << eps_float << std::endl;

  /***** OUTPUTS *****/
  // eps_double = 2.220446049250313e-16
  // eps_float = 1.192093e-07
  return EXIT_SUCCESS;
}
