#include "mathx.hpp"
#include "goodrand.hpp"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <chrono>
#include <tuple>

mathx::array<double> PR(mathx::matrix<double>, double);

int main(){

  mathx::matrix<double> T = {{0,1,0,0},
                             {1,0,1,0},
                             {1,1,0,1},
                             {0,0,1,0}};
  double d = 0.5;
  auto rank = PR(T, d);
  std::cout << rank.to_string() << std::endl;  

  return EXIT_SUCCESS;
}

mathx::array<double> PR(mathx::matrix<double> A, double d){
  mathx::matrix<double> A_hat = A;
  for(int i = 0; i < A_hat.rows();i++){
    for(int j = 0; j < A_hat.cols(); j++){
      A_hat[i][j] = (1-d) + d*A[i][j];
    }
  }

  mathx::array<double> v0(A_hat.rows(), 1.0);

  std::pair<double, mathx::array<double>> pageRank = mathx::linsolv::power_method(A_hat, v0, std::pow(10, -16), 10000, true);
  std::cout << "A_hat: " << std::endl;
  std::cout << A_hat.to_string() << std::endl;
  std::cout << "Page Rank (a,b,c,d,e,f)^-1" << std::endl;
  std::cout << "Page Rank eigenvalue " << pageRank.first << std::endl;
  return pageRank.second;
}