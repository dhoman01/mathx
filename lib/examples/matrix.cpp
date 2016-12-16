#include "matrix.hpp"
#include <cstdlib>
#include <iostream>

int main(){
  /* Declare a matrix */
  mathx::matrix<double> A(5,5);

  /* Use intializer list */
  mathx::matrix<double> I = {{1,0,0,0,0},{0,1,0,0,0},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,0,1}};

  /* Fill array with numbers
     i in [0,4] and j in [0,4]
     with a step of 1          */
  for(int i = 0; i <= 4; i++)
    for(int j = 0; j <= 4; j++)
      A[i][j] = j;

  /* Get number of rows */
  int m = A.rows();

  /* Get number of cols */
  int n = A.cols();

  /* Print a string representation
     of the matrix */
  std::cout << "A" << std::endl;
  std::cout << A.to_string() << std::endl;


  /* Print the One-Norm of A */
  std::cout << "||A||_1 = " << A.one_norm() << std::endl;

  /* Print the Infinity-Norm of A */
  std::cout << "||A||_max = " << A.infinity_norm() << std::endl;

  /* Print the One-Norm of I */
  std::cout << "||I||_1 = " << I.one_norm() << std::endl;

  /* Print the Infinity-Norm of I */
  std::cout << "||I||_max = " << I.infinity_norm() << std::endl;

  /***** OUTPUT *****/
  // A
  // 0          1          2          3          4
  // 0          1          2          3          4
  // 0          1          2          3          4
  // 0          1          2          3          4
  // 0          1          2          3          4

  // ||A||_1 = 20
  // ||A||_max = 10
  // ||I||_1 = 1
  // ||I||_max = 1
  return EXIT_SUCCESS;
}
