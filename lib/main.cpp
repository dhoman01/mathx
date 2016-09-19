/** @file main.cpp */
/** @mainpage
* Use Mathx to solve your computantional mathematic's problems.
* @author Dustin E. Homan
* @date 09/05/2016
* @version 0.1
*/
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include "mathx.hpp"

template <typename T>
void printElement(char align, T t, const int& width, char fill) {
  std::cout << (align == 'L' ? std::left : std::right) << std::setw(width)
            << std::setfill(fill) << t;
}

/**
* @brief This is the main driver for Mathx.
* @details Use this to solve your problems. It includes all the relevant
* headers. Referer to the documentation to determine what fundtion you need.
*/
int main() {
  mathx::complex<double> w(3.0001, 4.0001);
  mathx::complex<double> z(2.9999, 3.9999);
  double u[5] = {1, 1, -1.5, 100, 100};
  double v[5] = {0.99, 1.01, -1.2, 99.99, 99};
  printElement('R', "u", 6, ' ');
  printElement('R', "v", 6, ' ');
  printElement('R', "Absolute Error", 15, ' ');
  printElement('R', "Relative Error", 15, ' ');
  std::cout << std::endl;
  printElement('R', "-", 45, '-');
  std::cout << std::endl;
  for (int i = 0; i < 5; i++) {
    printElement('R', u[i], 6, ' ');
    printElement('R', v[i], 6, ' ');
    printElement('R', mathx::utils::error::e_abs(u[i], v[i]), 15, ' ');
    printElement('R', mathx::utils::error::e_rel(u[i], v[i]), 15, ' ');
    std::cout << std::endl;
  }

  std::cout << std::endl;

  std::cout << "w=" << w << std::endl;
  std::cout << "z=" << z << std::endl;
  std::cout << "|w|=" << mathx::utils::absolute_value(w) << std::endl;
  std::cout << "w+z=" << w + z << std::endl;
  std::cout << "w-z=" << w - z << std::endl;
  std::cout << "e_abs>" << mathx::utils::error::e_abs(w, z) << std::endl;
  std::cout << "e_rel>" << mathx::utils::error::e_rel(w, z) << std::endl;

  return EXIT_SUCCESS;
}
