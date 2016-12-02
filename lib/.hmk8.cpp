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
void problem5(bool debug=false);

double cast2ms(std::chrono::duration<double> duration);
mathx::matrix<double> generateRandom(int n, bool debug = false);
mathx::matrix<double> generateSPD(int n, bool debug = false);

void doPrintRunningTime(std::chrono::duration<double> child_runningtime);

int main(){
  problem1();
  problem2();
  problem3();
  problem4();
  problem5();
  return EXIT_SUCCESS;
}

void problem1(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem One                   " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  for(int n = 1000; n <= 2000; n+=100){
    std::cout << "\nn = " << n << std::endl;
    mathx::matrix<double> A(n, n, true);
    mathx::array<double> x(n, 1);
    mathx::array<double> b = mathx::linsolv::product(A, x);
    mathx::array<double> x0(n, goodrand::getRand(-5, 5));

    auto start = std::chrono::steady_clock::now();
    mathx::array<double> xstar = mathx::linsolv::jacobi(A, b, x0, std::pow(10, -8), 1000);
    auto end = std::chrono::steady_clock::now();

    std::cout << "error " << mathx::vectors::euclideanLength(x - xstar) << std::endl;
    if(debug){
      std::cout << "x*: ";
      xstar.to_string();
      std::cout << std::endl;
    }
    doPrintRunningTime(end - start);
  }
}

void problem2(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem Two                   " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  for(int n = 1000; n <= 2000; n+=100){
    std::cout << "\nn = " << n << std::endl;
    mathx::matrix<double> A(n, n, true);
    mathx::array<double> x(n, 1);
    mathx::array<double> b = mathx::linsolv::product(A, x);
    mathx::array<double> x0(n, goodrand::getRand(-5, 5));

    auto start = std::chrono::steady_clock::now();
    mathx::array<double> xstar = mathx::linsolv::gauss_seidel(A, b, x0, std::pow(10, -8), 1000);
    auto end = std::chrono::steady_clock::now();

    std::cout << "error " << mathx::vectors::euclideanLength(x - xstar) << std::endl;
    if(debug){
      std::cout << "x*: ";
      xstar.to_string();
      std::cout << std::endl;
    }
    doPrintRunningTime(end - start);
  }
}

void problem3(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem Three                 " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  for(int n = 1000; n <= 10000; n+=1000){
    mathx::matrix<double> A(n, n, true);
    mathx::array<double> x(n, 1);
    mathx::array<double> b = mathx::linsolv::product(A, x);
    mathx::array<double> x0(n, 5);

    std::cout << "Jacobi (n, iter)" << std::endl;
    mathx::array<double> Jxstar = mathx::linsolv::jacobi(A, b, x0, std::pow(10, -8), 1000, true);
    std::cout << "Gauss-Seidel (n, iter)" << std::endl;
    mathx::array<double> Gxstar = mathx::linsolv::gauss_seidel(A, b, x0, std::pow(10, -8), 1000, true);

    if(debug){
      std::cout << "Jacobi x*: ";
      Jxstar.to_string();
      std::cout << std::endl;
      std::cout << "G-S x*: ";
      Gxstar.to_string();
      std::cout << std::endl;
    }
  }
}

void problem4(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem Four                  " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "\nn,LU,Jacobi,Gauss-Seidel" << std::endl;
  for(int n = 1000; n <= 5000; n+=500){
    if(debug) std::cout << "Generating Matrix with " << n << " rows" << std::endl;
    mathx::matrix<double> A = generateRandom(n, debug);
    // mathx::matrix<double> A(n,n,true);
    mathx::array<double> x(n, 1);
    mathx::array<double> b = mathx::linsolv::product(A, x);
    mathx::array<double> x0(n, 2);
    mathx::matrix<double> LU;

    if(debug) std::cout << "Solving system using LU Factorization" << std::endl;
    auto luStart = std::chrono::steady_clock::now();
    mathx::array<double> luXstar = mathx::linsolv::solve(A, b, LU, 0);
    auto luEnd = std::chrono::steady_clock::now();
    if(debug) std::cout << "Solved system with error " << mathx::vectors::euclideanLength(luXstar - x) << std::endl;

    if(debug) std::cout << "Solving system using Jacobi Iteration" << std::endl;
    auto jacStart = std::chrono::steady_clock::now();
    mathx::array<double> jacXstar = mathx::linsolv::jacobi(A, b, x0, std::pow(10, -8), 1000);
    auto jacEnd = std::chrono::steady_clock::now();
    if(debug) std::cout << "Solved system with error " << mathx::vectors::euclideanLength(jacXstar - x) << std::endl;

    if(debug) std::cout << "Solving system using Gauss-Seidel Iteration" << std::endl;
    auto gsStart = std::chrono::steady_clock::now();
    mathx::array<double> gsXstar = mathx::linsolv::gauss_seidel(A, b, x0, std::pow(10, -8), 1000);
    auto gsEnd = std::chrono::steady_clock::now();
    if(debug) std::cout << "Solved system with error " << mathx::vectors::euclideanLength(gsXstar - x) << std::endl;

    std::cout << n << ", " << cast2ms(luEnd - luStart) << ", " << cast2ms(jacEnd - jacStart) << ", " << cast2ms(gsEnd - gsStart) << std::endl;

    if(debug){
      std::cout << "LU x*: ";
      luXstar.to_string();
      std::cout << std::endl;
      std::cout << "Jacobi x*: ";
      jacXstar.to_string();
      std::cout << std::endl;
      std::cout << "G-S x*: ";
      gsXstar.to_string();
      std::cout << std::endl;
    }
  }
}

void problem5(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem Five                   " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  for(int n = 10; n <= 160; n *= 2){
    std::cout << "\nn = " << n << std::endl;
    mathx::matrix<double> A = generateSPD(n);
    mathx::array<double> x(n, 1);
    mathx::array<double> b = mathx::linsolv::product(A, x);
    mathx::array<double> x0(n, 2);

    auto start = std::chrono::steady_clock::now();
    mathx::array<double> xstar = mathx::linsolv::cgm(A, b, x0, std::pow(10, -8), 1000);
    auto end = std::chrono::steady_clock::now();

    std::cout << "error " << mathx::vectors::euclideanLength(x - xstar) << std::endl;
    if(debug){
      std::cout << "x*: ";
      xstar.to_string();
      std::cout << std::endl;
    }
    doPrintRunningTime(end - start);
  }
}

double cast2ms(std::chrono::duration<double> duration){
  return std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
}

mathx::matrix<double> generateRandom(int n, bool debug){
  if(goodrand::getRand(0.0, 1.0) < .75){
    if(debug) std::cout << "Generating a random diagonally dominant matrix" << std::endl;
    return mathx::matrix<double>(n,n,true);
  } else {
    mathx::matrix<double> A = mathx::matrix<double>(n,n);
    if(debug) std::cout << "Generating a random matrix of dim(" << A.rows() << ", " << A.cols() << ")" << std::endl;
    for(int i = 0; i < n; i++)
      for(int j = 0; j < n; j++)
        A[i][j] = goodrand::getRand(-5.0, 5.0);

    return A;
  }
}

mathx::matrix<double> generateSPD(int n, bool debug){
  mathx::matrix<double> A = mathx::matrix<double>(n,n, true);

  if(debug)
    std::cout << "A: \n" + A.to_string() << std::endl;

  return mathx::linsolv::mult_transpose(A);
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
