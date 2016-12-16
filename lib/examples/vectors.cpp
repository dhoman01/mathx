#include "vectors.hpp"
#include <cstdlib>
#include <iostream>

int main(){
  /* Define vectors v and w */
  mathx::array<double> v(3, 2);
  mathx::array<double> w(3, 6);

  /* Calculate dot product of two vectors, v and w */
  double dot = mathx::vectors::dot_product(v, w);

  /* Calculate cross product of two vectors, v and w */
  mathx::array<double> cross = mathx::vectors::cross_product(v, w);

  /* Calculate l-2 norm of vector v */
  double norm = mathx::vectors::norm(v);

  /* Calculate l-1 norm of vector v */
  double one_norm = mathx::vectors::one_norm(v);

  /* Calculate l-max norm of vector v */
  double infinity_norm = mathx::vectors::infinity_norm(v);

  /* Normalize vector v */
  mathx::array<double> v_norm = mathx::vectors::normalize(v);

  /* Print results */
  std::cout << "<v,w>     = " << dot << std::endl;
  std::cout << "v X w     = " << cross.to_string() << std::endl;
  std::cout << "||v||_2   = " << norm << std::endl;
  std::cout << "||v||_1   = " << one_norm << std::endl;
  std::cout << "||v||_max = " << infinity_norm << std::endl;
  std::cout << "normalized: " << v_norm.to_string() << std::endl;

  /***** OUTPUTS *****/
  // <v,w>     = 36
  // v X w     = [ 0 0 0  ]^T
  // ||v||_2   = 3.4641
  // ||v||_1   = 6
  // ||v||_max = 2
  // normalized: [ 0.57735 0.57735 0.57735  ]^T

  return EXIT_SUCCESS;
}
