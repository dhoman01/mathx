#include "mathx.hpp"
#include "goodrand.hpp"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <chrono>
#include <tuple>

void problem1(bool debug=false);
void problem2(bool debug=false);
void problem3(bool debug=false);
void problem4(bool debug=false);

void doPrintRunningTime(std::chrono::duration<double> child_runningtime);

int main(){
  // problem1();
  // problem2();
  // problem3();
  problem4();
  return EXIT_SUCCESS;
}

void problem1(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem One                   " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  for(int n = 10; n <= 160; n *= 2){
    mathx::matrix<double> A(n, n, true);
    mathx::array<double> v0(n, goodrand::getRand(-1.0, 1.0));

    auto start = std::chrono::steady_clock::now();
    std::pair<double, mathx::array<double>> maxEigen = mathx::linsolv::power_method(A, v0, std::pow(10, -16), 10000);
    auto end = std::chrono::steady_clock::now();

    mathx::array<double> vstar = mathx::linsolv::product(A, maxEigen.second);

    if(debug) std::cout << "eigenvalue: " << maxEigen.first << std::endl;
    std::cout << "\nn = " << n << std::endl;
    // Av= ƛv
    std::cout << "Av=\u03BBv" << std::endl;
    std::cout << "error = " << mathx::vectors::euclideanLength(vstar - maxEigen.second * maxEigen.first) << std::endl;
    doPrintRunningTime(end - start);
  }
}

void problem2(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem Two                   " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  for(int n = 100; n <= 1600; n *= 2){
    mathx::matrix<double> A(n, n, true);
    mathx::array<double> v0(n, goodrand::getRand(-1.0, 1.0));

    auto start = std::chrono::steady_clock::now();
    std::pair<double, mathx::array<double>> minEigen = mathx::linsolv::inverse_power_method(A, v0, goodrand::getRand(-1.0, 1.0), std::pow(10, -16), 10000);
    auto end = std::chrono::steady_clock::now();

    mathx::array<double> vstar = mathx::linsolv::product(A, minEigen.second);

    if(debug) std::cout << "eigenvalue: " << minEigen.first << std::endl;
    std::cout << "\nn = " << n << std::endl;
    // Av=lambda * v
    std::cout << "Av=\u03BBv" << std::endl;
    std::cout << "error = " << mathx::vectors::euclideanLength(vstar - minEigen.second * minEigen.first) << std::endl;
    doPrintRunningTime(end - start);
  }
}

void problem3(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem Three                 " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  for(int n = 10; n <= 160; n *= 2){
    mathx::matrix<double> A(n, n, true);
    mathx::array<double> v0(n, goodrand::getRand(-1.0, 1.0));

    std::pair<double, mathx::array<double>> maxEigen = mathx::linsolv::power_method(A, v0, std::pow(10, -16), 10000, true);
  }
}

void problem4(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem Four                  " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  for(int n = 100; n <= 1600; n *= 2){
    mathx::matrix<double> A(n, n, true);
    mathx::array<double> v0(n, goodrand::getRand(-1.0, 1.0));

    std::pair<double, mathx::array<double>> minEigen = mathx::linsolv::inverse_power_method(A, v0, goodrand::getRand(-1.0, 1.0), std::pow(10, -16), 10000, true);
  }
}

void doPrintRunningTime(std::chrono::duration<double> child_runningtime)
{
  auto sec = std::chrono::duration_cast<std::chrono::seconds>(child_runningtime);
  child_runningtime -= sec;
  auto mill = std::chrono::duration_cast<std::chrono::milliseconds>(child_runningtime);
  child_runningtime -= mill;
  auto micro = std::chrono::duration_cast<std::chrono::microseconds>(child_runningtime);
  std::cout << "\n***     Time spent executing: " << sec.count() << " seconds ";
  std::cout << mill.count() << " milliseconds and " << micro.count() << " microseconds     ***" << std::endl;
}
