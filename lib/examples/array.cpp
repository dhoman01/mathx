#include "array.hpp"
#include <cstdlib>
#include <iostream>

int main(){
  /* Declare an array */
  mathx::array<double> v;

  /* Use intializer list */
  mathx::array<double> w = {11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};

  /* Push i in [0,10] with
     step of 1 into array */
  for(double i = 0; i <= 10; i += 1)
    v.push(i);

  /* Multiply every element
     by 2, * operand */
  v = v * 2;

  /* Print a string representation
     of the array */
  std::cout << "v: " << v.to_string() << std::endl;


  /* Dot to arrays together using * operand */
  double dot = v * w;

  /* Print the result of <v,w> */
  std::cout << "<v,w> = " << dot << std::endl;

  /***** OUTPUT *****/
  // v: [ 0 2 4 6 8 10 12 14 16 18 20  ]^T
  // <v,w> = 1980
  return EXIT_SUCCESS;
}
