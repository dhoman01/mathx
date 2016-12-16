#include "goodrand.hpp"
#include <cstdlib>
#include <iostream>

int main(){
  /* Get a random double in the range
   * [0,1] with a uniform distribution
   */
  double random_double = goodrand::get_rand(0.0, 1.0);

  /* Get a random integer in the range
   * [-10,10] with a uniform distribution
   */
  int random_int = goodrand::get_rand(-10,10);

  /* Print random double */
  std::cout << "random double is " << random_double << std::endl;

  /* Print random integer */
  std::cout << "random integer is " << random_int << std::endl;

  /***** OUTPUT *****/
  // random double is 0.224257
  // random integer is -4
  return EXIT_SUCCESS;
}
