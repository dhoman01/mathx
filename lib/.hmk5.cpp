


#include "mathx.hpp"
#include "goodrand.hpp"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <chrono>

void problem1(bool debug = false);
void problem2(bool debug = false);
void problem3(bool debug = false);
void problem4(bool debug = false);

void doPrintRunningTime(std::chrono::duration<double> child_runningtime);
void generateDiagonallyDom(int n, mathx::matrix<double>& A, mathx::array<double>& b, bool debug = false);
void generateSym(int n, mathx::matrix<double>& A, mathx::array<double>& b, bool debug = false);

int main(){
  problem1();
  std::cout << "\n\n\n";
  problem2();
  std::cout << "\n\n\n";
  problem3();
  std::cout << "\n\n\n";
  problem4(true);

  return EXIT_SUCCESS;
}

void problem1(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem One                   " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  for(int n = 10; n <= 160; n *= 2){
    std::cout << "\nn = " << n << std::endl;
    mathx::matrix<double> A;
    mathx::array<double> b;

    // Generate a diagonaly dominant matrix and multiply by x = 1
    generateDiagonallyDom(n, A, b, debug);

    // Perform LU factorization with partial pivoting then
    // use Foward Substitution and Back Substitution to find solution
    mathx::matrix<double> LU;
    auto start = std::chrono::steady_clock::now();
    mathx::array<double> x_1 = mathx::linsolv::solve(A, b, LU, 2);
    auto end = std::chrono::steady_clock::now();

    if(debug){
      std::cout << "x_1: " << std::endl;
      for(int i = 0; i < x_1.size(); i++)
        std::cout << x_1[i] << std::endl;
    }

    doPrintRunningTime(end - start);
    mathx::array<double> x(n, 1);
    std::cout << "Error (one-norm): " << mathx::vectors::oneNorm(x - x_1) << std::endl;
    std::cout << "Error (two-norm): " << mathx::vectors::euclideanLength(x - x_1) << std::endl;
    std::cout << "Error (max-norm): " << mathx::vectors::maxNorm(x - x_1) << std::endl;
    std::cout << "LU: " << (LU.hasPivoted ? "has pivoted" : "has not pivoted") << std::endl;
    if(debug){
      for(int i = 0; i < LU.rows(); i++) {
        for(int j = 0; j < LU.cols(); j++)
          std::cout << std::setw(10) << std::left << LU.get(i,j) << " ";
        std::cout << std::endl;
      }
      std::cout << std::endl;
    }
  }
}

void problem2(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem Two                   " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  for(int n = 10; n <= 160; n *= 2){
    std::cout << "\nn = " << n << std::endl;
    mathx::matrix<double> A;
    mathx::array<double> b;

    // Generate a diagonaly dominant matrix and multiply by x = 1
    generateDiagonallyDom(n, A, b, debug);

    // Perform GE with partial pivoting then
    // use Back Substitution to find solution
    auto start = std::chrono::steady_clock::now();
    mathx::array<double> x_1 = mathx::linsolv::solve(A, b, 2);
    auto end = std::chrono::steady_clock::now();

    if(debug){
      std::cout << "x_1: " << std::endl;
      for(int i = 0; i < x_1.size(); i++)
        std::cout << x_1[i] << std::endl;
    }
    doPrintRunningTime(end - start);
    mathx::array<double> x(n, 1);
    std::cout << "Error (one-norm): " << mathx::vectors::oneNorm(x - x_1) << std::endl;
    std::cout << "Error (two-norm): " << mathx::vectors::euclideanLength(x - x_1) << std::endl;
    std::cout << "Error (max-norm): " << mathx::vectors::maxNorm(x - x_1) << std::endl;
    std::cout << "A: " << (A.hasPivoted ? "has pivoted" : "has not pivoted") << std::endl;
    if(debug){
      for(int i = 0; i < A.rows(); i++) {
        for(int j = 0; j < A.cols(); j++)
          std::cout << std::setw(10) << std::left << A.get(i,j) << " ";
        std::cout << std::endl;
      }
      std::cout << std::endl;
    }
  }
}

void problem3(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem Three                   " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  for(int n = 10; n <= 160; n *= 2){
    std::cout << "\nn = " << n << std::endl;
    mathx::matrix<double> A;
    mathx::array<double> b;
    generateSym(n, A, b, debug);

    mathx::matrix<double> LU;
    auto start = std::chrono::steady_clock::now();
    mathx::array<double> x_1 = mathx::linsolv::solve(A, b, LU, 2);
    auto end = std::chrono::steady_clock::now();

    if(debug){
      std::cout << "x_1: " << std::endl;
      for(int i = 0; i < x_1.size(); i++)
        std::cout << x_1[i] << std::endl;
    }
    doPrintRunningTime(end - start);
    mathx::array<double> x(n, 1);
    std::cout << "Error (one-norm): " << mathx::vectors::oneNorm(x - x_1) << std::endl;
    std::cout << "Error (two-norm): " << mathx::vectors::euclideanLength(x - x_1) << std::endl;
    std::cout << "Error (max-norm): " << mathx::vectors::maxNorm(x - x_1) << std::endl;
    std::cout << "LU: " << (LU.hasPivoted ? "has pivoted" : "has not pivoted") << std::endl;
    if(debug){
      for(int i = 0; i < LU.rows(); i++) {
        for(int j = 0; j < LU.cols(); j++)
          std::cout << std::setw(10) << std::left << LU.get(i,j) << " ";
        std::cout << std::endl;
      }
      std::cout << std::endl;
    }
  }
}

void problem4(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem Four                   " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  mathx::matrix<int> A(4,4,0);
  A.set(0,0,5);
  A.set(0,1,6);
  A.set(0,2,7);
  A.set(0,3,8);
  A.set(1,1,4);
  A.set(1,2,3);
  A.set(1,3,2);
  A.set(2,3,1);
  A.set(3,2,-1);
  A.set(3,3,-2);
  mathx::array<int> b = {26,9,1,-3};

  auto start = std::chrono::steady_clock::now();
  mathx::linsolv::lu(A,b,1);
  auto end = std::chrono::steady_clock::now();
  if(debug){
    std::cout << "LU: " << std::endl;
    for(int i = 0; i < A.rows(); i++) {
      for(int j = 0; j < A.cols(); j++)
        std::cout << std::setw(10) << std::left << A.get(i,j) << " ";
      std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "b: " << std::endl;
    for(int i = 0; i < b.size(); i++)
      std::cout << b[i] << std::endl;
  }
  auto start2 = std::chrono::steady_clock::now();
  mathx::array<int> y = mathx::linsolv::forward_substitution(A, b, true);
  mathx::array<int> x_1 = mathx::linsolv::back_substitution(A, y);
  auto end2 = std::chrono::steady_clock::now();

  if(debug){
    std::cout << "\nx_1: " << std::endl;
    for(int i = 0; i < x_1.size(); i++)
      std::cout << x_1[i] << std::endl;
  }
  doPrintRunningTime((end - start) + (end2 -start2));
  mathx::array<int> x(4, 1);
  std::cout << "Error (one-norm): " << mathx::vectors::oneNorm(x - x_1) << std::endl;
  std::cout << "Error (two-norm): " << mathx::vectors::euclideanLength(x - x_1) << std::endl;
  std::cout << "Error (max-norm): " << mathx::vectors::maxNorm(x - x_1) << std::endl;
}

void generateDiagonallyDom(int n, mathx::matrix<double>& A, mathx::array<double>& b, bool debug){
  A = mathx::matrix<double>(n, n,true);
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
  b = mathx::linsolv::product(A, x);
  if(debug){
    std::cout << "b: " << std::endl;
    for(int i = 0; i < b.size(); i++)
      std::cout << b[i] << std::endl;

    std::cout << std::endl;
  }
}

void generateSym(int n, mathx::matrix<double>& A, mathx::array<double>& b, bool debug){
  A = mathx::matrix<double>(n,n);
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      A.set(i, j, ( i > j ? A.get(j, i) : goodrand::getRand(1.0, 2.0)));
    }
  }

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
  b = mathx::linsolv::product(A, x);
  if(debug){
    std::cout << "b: " << std::endl;
    for(int i = 0; i < b.size(); i++)
      std::cout << b[i] << std::endl;

    std::cout << std::endl;
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
