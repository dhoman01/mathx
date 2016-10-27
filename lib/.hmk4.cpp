#include "mathx.hpp"
#include <iostream>
#include <iomanip>
#include <cstdlib>

void problem5(int, bool debug = false);

int main(){
  for(int n = 10; n <= 160; n *= 2)
    problem5(n);

  return EXIT_SUCCESS;
}

void problem5(int n, bool debug){
    std::cout << "n = " << n << std::endl;
    mathx::matrix<double> A(n, n,true);

    if(debug){
      std::cout << "A: " << std::endl;
      for(int i = 0; i < A.rows(); i++) {
        for(int j = 0; j < A.cols(); j++)
          std::cout << std::setw(10) << std::left << A.get(i,j) << " ";
        std::cout << std::endl;
      }
      std::cout << std::endl;
    }


    mathx::array<double> x(n, 1);
    mathx::array<double> b = mathx::linsolv::product(A, x);
    if(debug){
      std::cout << "b: " << std::endl;
      for(int i = 0; i < b.size(); i++)
        std::cout << b[i] << std::endl;

      std::cout << std::endl;
    }

    mathx::array<double> x_1 = mathx::linsolv::solve(A,b);
    if(debug){
      std::cout << "x_1: " << std::endl;
      for(int i = 0; i < x_1.size(); i++)
        std::cout << x_1[i] << std::endl;
    }

    std::cout << "Error (one-norm): " << mathx::vectors::oneNorm(x - x_1) << std::endl;
    std::cout << "Error (two-norm): " << mathx::vectors::euclideanLength(x - x_1) << std::endl;
    std::cout << "Error (max-norm): " << mathx::vectors::maxNorm(x - x_1) << std::endl;
}
