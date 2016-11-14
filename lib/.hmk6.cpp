#include "mathx.hpp"
#include "goodrand.hpp"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <chrono>
#include <tuple>

void problem1(bool debug = false);
void problem2(bool debug = false);
void problem3(bool debug = false);
void problem4(bool debug = false);
void problem5(bool debug = false);
void problem6(bool debug = false);
void problem7(bool debug = false);

void generateSym(int n, mathx::matrix<double>& A, bool debug = false);
void generateSPD(int n, mathx::matrix<double>& A, bool debug = false);
void doPrintRunningTime(std::chrono::duration<double> child_runningtime);
std::tuple<mathx::array<double>,mathx::array<double>,mathx::array<double>> generateTriDiag(int n);

int main(){
  problem1();
  std::cout << std::endl;
  problem2();
  std::cout << std::endl;
  problem3();
  std::cout << std::endl;
  problem4();
  std::cout << std::endl;
  problem5();
  std::cout << std::endl;
  problem6();
  std::cout << std::endl;
  problem7();
  return EXIT_SUCCESS;
}

void problem1(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem One                   " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  for(int n = 10; n <= 160; n *= 2){
    std::cout << "\nn=" << n << std::endl;
    mathx::matrix<double> A;
    generateSPD(n, A, debug);

    mathx::array<double> x(n, 1);
    mathx::array<double> b = mathx::linsolv::product(A, x);
    if(debug){
      std::cout << "A"<< std::endl;
      std::cout << A.to_string() << std::endl;
      std::cout << "b" << std::endl;
      for(int i = 0; i < b.size(); i++)
        std::cout << b[i] << std::endl;
    }

    auto start = std::chrono::steady_clock::now();
    mathx::array<double> x_1 = mathx::linsolv::solve(A, b);
    auto end = std::chrono::steady_clock::now();

    std::cout << "Error (one-norm): " << mathx::vectors::oneNorm(x - x_1) << std::endl;
    std::cout << "Error (two-norm): " << mathx::vectors::euclideanLength(x - x_1) << std::endl;
    std::cout << "Error (max-norm): " << mathx::vectors::maxNorm(x - x_1) << std::endl;
    doPrintRunningTime(end - start);
  }
}

void problem2(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem Two                   " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  mathx::matrix<double> A;
  generateSym(4, A, debug);
  std::cout << A.to_string() << std::endl;
  std::cout << "Matrix is " << (mathx::linsolv::is_spd(A) ? "" : " not ") << "s.p.d." << std::endl;
  std::cout << std::endl;
  generateSPD(4, A, debug);
  std::cout << A.to_string() << std::endl;
  std::cout << "Matrix is " << (mathx::linsolv::is_spd(A) ? "" : " not ") << "s.p.d." << std::endl;
  if(debug)
    std::cout << "A\n" + A.to_string() << std::endl;
}

void problem3(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem Three                 " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  mathx::matrix<double> A(3,3,true);
  if(debug){
    for(int j = 0; j < A.cols(); j++){
      std::cout << "Column " << j << ": ";
      double sum = 0;
      for(int i = 0; i < A.rows(); i++){
        std::cout << A[i][j] << (i < A.rows() - 1 ? " + " : " = ");
        sum += A[i][j];
      }
      std::cout << sum << std::endl;
    }
  }
  std::cout << "A:\n" + A.to_string() << std::endl;
  std::cout << "\n||A||_1 = " << A.one_norm() << std::endl;
}

void problem4(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem Four                  " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  mathx::matrix<double> A(3,3,true);
  if(debug){
    for(int i = 0; i < A.rows(); i++){
      std::cout << "Row " << i << ": ";
      double sum = 0;
      for(int j = 0; j < A.cols(); j++){
        std::cout << A[i][j] << (j < A.cols() - 1 ? " + " : " = ");
        sum += A[i][j];
      }
      std::cout << sum << std::endl;
    }
  }
  std::cout << "A:\n" + A.to_string() << std::endl;
  std::cout << "\n||A||_infty = " << A.infinity_norm() << std::endl;
}

void problem5(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem Five                  " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  mathx::matrix<double> A;
  generateSym(10, A);
  std::cout << "A:\n" + A.to_string() << std::endl;
  std::cout << "k(A) = " << mathx::linsolv::kappa(A) << std::endl;
}

void problem6(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem Six                   " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "n,count" << std::endl;
  for(int n = 10; n <= 2560; n *= 2){
      mathx::matrix<double> A;
      generateSym(n, A);
      mathx::linsolv::kappa(A);
  }
}

void problem7(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem Seven                 " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  for(int n = 10; n <= 160; n *= 2){
    std::tuple<mathx::array<double>,mathx::array<double>,mathx::array<double>> A = generateTriDiag(n);
    mathx::array<double> x (n,1);
    mathx::array<double> b = mathx::linsolv::product(std::get<0>(A), std::get<1>(A), std::get<2>(A), x);

    if(debug){
      std::cout << "al: ";
      for(int j = 0; j < b.size(); j++){
        std::cout << std::get<0>(A)[j] << " ";
      }
      std::cout << std::endl;
      std::cout << "am: ";
      for(int i = 0; i < b.size(); i++){
        std::cout << std::get<1>(A)[i] << " ";
      }
      std::cout << std::endl;
      std::cout << "au: ";
      for(int i = 0; i < b.size(); i++){
        std::cout << std::get<2>(A)[i] << " ";
      }
      std::cout << std::endl;
      std::cout << "b: ";
      for(int i = 0; i < b.size(); i++){
        std::cout << b[i] << " ";
      }
      std::cout << std::endl;
    }

    auto start = std::chrono::steady_clock::now();
    mathx::array<double> x_1 = mathx::linsolv::solve(std::get<0>(A), std::get<1>(A), std::get<2>(A), b);
    auto end = std::chrono::steady_clock::now();

    if(debug){
      std::cout << "x_1: ";
      for(int i = 0; i < x_1.size(); i++){
        std::cout << x_1[i] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << "\nn=" << n << std::endl;
    std::cout << "Error (one-norm): " << mathx::vectors::oneNorm(x - x_1) << std::endl;
    std::cout << "Error (two-norm): " << mathx::vectors::euclideanLength(x - x_1) << std::endl;
    std::cout << "Error (max-norm): " << mathx::vectors::maxNorm(x - x_1) << std::endl;
    doPrintRunningTime(end - start);
  }
}

std::tuple<mathx::array<double>,mathx::array<double>,mathx::array<double>> generateTriDiag(int n){
  // Initialize containers and set values
  mathx::array<double> al = mathx::array<double>(n, 1);
  mathx::array<double> au = mathx::array<double>(n, 1);
  mathx::array<double> am = mathx::array<double>(n, -2);

  // Set padding
  al[0] = 0;
  au[n - 1] = 0;

  return std::make_tuple(al,am,au);
}

void generateSym(int n, mathx::matrix<double>& A, bool debug){
  A = mathx::matrix<double>(n,n);
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      A[i][j] = ( i > j ? A[j][i] : goodrand::getRand(1.0, 2.0));
    }
  }

  if(debug)
    std::cout << "A: \n" + A.to_string() << std::endl;
}

void generateSPD(int n, mathx::matrix<double>& A, bool debug){
  A = mathx::matrix<double>(n,n, true);
  for(int i = 0; i < n; i++){
    for(int j = i + 1; j < n; j++){
      A[i][j] = A[j][i];
    }
  }

  if(debug)
    std::cout << "A: \n" + A.to_string() << std::endl;
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
