#include <cstdlib>  // Gets machines EXIT_SUCCESS
#include <iostream>
#include <iomanip>
#include "mathx.hpp"

void problem7();

int main(){
  auto f = [](double x){ return std::sqrt(x) - 1.1; };
  double root = mathx::roots::bisect(f, 0.0, 2.0, f(0), f(2), std::pow(10,-8), 50);
  std::cout << "root " << root << std::endl;
  std::cout << "e_abs " << mathx::utils::error::e_abs(1.21, root) << std::endl;
  std::cout << std::endl;
  problem7();
  return EXIT_SUCCESS;
}

void problem7(){
  mathx::array<double> coeff = {-512, 2304, -4608, 5376, -4032, 2016, -672, 144, -18, 1};
  std::cout << std::setw(10) << "i," << std::setw(13) << "Nested Eval," << std::setw(13) << " Direct Eval" << std::endl;
  for(double i = 1.92; i <= 2.08; i+= (2.08 - 1.92)/161)
    std::cout << std::setw(10) << i << "," << std::setw(13) << mathx::utils::poly::nested_eval(coeff, i) << ", "  << std::setw(13) <<  std::pow((i - 2), 9) << std::endl;
}
