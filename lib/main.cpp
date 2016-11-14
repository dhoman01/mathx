#include "mathx.hpp"
#include "goodrand.hpp"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <chrono>
#include <tuple>

void problem1(bool debug = false);

mathx::matrix<double> generate_nonsquare(int n, int m);
void doPrintRunningTime(std::chrono::duration<double> child_runningtime);

int main(){
  // problem1();
  mathx::matrix<double> A = {{1.,0.},{1.,1.},{1.,2.}};
  // mathx::matrix<double> A = {{1.,2.,1.},{0.,1.,2.},{1.,2.,0.}};
  std::cout << "A" << std::endl;
  std::cout << A.to_string() << std::endl;
  std::pair<mathx::matrix<double>, mathx::matrix<double>> p = mathx::linsolv::qr_factorization(A);
  std::cout << "Q" << std::endl;
  std::cout << p.first.to_string() << std::endl;
  std::cout << "R" << std::endl;
  std::cout << p.second.to_string() << std::endl;
  return EXIT_SUCCESS;
}

void problem1(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem One                   " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  for(int i = 10; i <= 160; i *= 2){
    int n = i;
    int m = i - goodrand::getRand(2, i - 1);
    std::cout << "(" << n << " x " << m << ")" << std::endl;
    mathx::matrix<double> A = generate_nonsquare(n,m);
    mathx::array<double> x(m, 1);
    mathx::array<double> b = mathx::linsolv::product(A,x);

    auto start = std::chrono::steady_clock::now();
    mathx::array<double> x_1 = mathx::linsolv::least_squares(A,b);
    mathx::array<double> r = b - mathx::linsolv::product(A, x);
    auto end = std::chrono::steady_clock::now();

    if(debug){
      std::cout << "A" << std::endl;
      std::cout << A.to_string() << std::endl;
      std::cout << "b: " << b.to_string() << std::endl;
      std::cout << "x: " << x.to_string() << std::endl;
      std::cout << "x_1: " << x_1.to_string() << std::endl;
    }

    std::cout << "Error (one-norm): " << mathx::vectors::oneNorm(x - x_1) << std::endl;
    std::cout << "Error (two-norm): " << mathx::vectors::euclideanLength(x - x_1) << std::endl;
    std::cout << "Error (max-norm): " << mathx::vectors::maxNorm(x - x_1) << std::endl;
    std::cout << std::endl;
    if(debug) std::cout << "r = " << r.to_string() << std::endl;
    std::cout << "Residual (two-norm): " << mathx::vectors::euclideanLength(r) << std::endl;
    doPrintRunningTime(end - start);
  }
}

mathx::matrix<double> generate_nonsquare(int n, int m){
  mathx::matrix<double> A(n, m, 0.);
  for(int i = 0; i < n; i++)
    for(int j = 0; j < m; j++)
      A[i][j] = goodrand::getRand(-15., 15.);

  return A;
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
