#include "interpolation.hpp"
#include <cstdlib>
#include <iostream>

int main(){
  /* Generate some data */
  int n = 4;
  auto f = [](double t){ return 3 * std::pow(t, 2) - 2 * std::pow(t, 3); };
  mathx::array<double> t;
  mathx::array<double> ft;
  double spacing = 3.0 / (n - 1);
  for(double i = -1.0; i <= 2.0; i += spacing){
    t.push(i);
    ft.push(f(i));
  }

  /* Get divided difference
  *  table for generated data
  */
  mathx::matrix<double> diff_table = mathx::interpolation::divided_differences(t, ft);

  /* Print table */
  std::cout << "diff_table" << std::endl;
  std::cout << diff_table.to_string() << std::endl;

  /* Extract Newton's form coefficients */
  mathx::array<double> coeff = mathx::interpolation::newtons_coeff(diff_table);

  /* Print the extracted coefficients */
  std::cout << n << ": coeff: " << coeff.to_string() << std::endl;

  /* Evaluate Newton's form poly at x */
  double fx = mathx::interpolation::eval_newtons(10, t, coeff);

  /* Print results */
  std::cout << "f(t)  = " << f(10) << std::endl;
  std::cout << "f(x)  = " << fx << std::endl;
  std::cout << "e_abs = " << std::abs(f(10) - fx) << std::endl;

  /***** OUTPUT *****/
  // diff_table
  // 5          0          0          0
  // 0          -5         0          0
  // 1          1          3          0
  // -4         -5         -3         -2

  // 4: coeff: [ 5 -5 3 -2  ]^T
  // f(t)  = -1700
  // f(x)  = -1700
  // e_abs = 0

  return EXIT_SUCCESS;
}
