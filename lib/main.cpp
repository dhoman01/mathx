#include "mathx.hpp"
#include "goodrand.hpp"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <chrono>
#include <tuple>

void problem1(bool debug = false);
void problem2(bool debug = false);
void problem4(bool debug = false);
void problem5(bool debug = false);
void problem6(bool debug = false);
void problem7(bool debug = false);
void problem8(bool debug = false);
void problem9(bool debug = false);

mathx::matrix<double> generate_nonsquare(int n, int m);
void doPrintRunningTime(std::chrono::duration<double> child_runningtime);

int main(){
  problem1();
  problem2();
  problem4();
  problem5();
  problem6();
  problem7();
  problem8();
  problem9();
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

void problem2(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem Two                   " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  for(int n = 2; n <= 10; n++){
    std::cout << "n = " << n << std::endl;
    mathx::matrix<double> A(n,n,true);
    mathx::matrix<double> Q = mathx::linsolv::qr_factorization(A);
    mathx::matrix<double> I = mathx::linsolv::mult_transpose(Q);
    if(debug){
      std::cout << "A" << std::endl;
      std::cout << A.to_string() << std::endl;
      std::cout << "Q" << std::endl;
      std::cout << Q.to_string() << std::endl;
      std::cout << "I" << std::endl;
      std::cout << I.to_string() << std::endl;
    }

    std::cout << "One-Norm: " << I.one_norm() << std::endl;
    std::cout << "Infty-Norm: " << I.infinity_norm() << std::endl;
  }
}

void problem4(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem Four                  " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  for(int n = 2; n <= 5; n++){
    std::cout << "n = " << n << std::endl;
    mathx::matrix<double> A(n,n,true);
    mathx::array<double> x(n,1);
    mathx::array<double> b = mathx::linsolv::product(A,x);

    auto start = std::chrono::steady_clock::now();
    mathx::array<double> x_1 = mathx::linsolv::least_squares_QR(A,b);
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

void problem5(bool debug){
    std::cout << "----------------------------------------------" << std::endl;
    std::cout << "                Problem Five                  " << std::endl;
    std::cout << "----------------------------------------------" << std::endl;
    std::cout << "\nn,l_2-error" << std::endl;
    for(int n = 10; n <= 160; n *= 2){
      mathx::matrix<double> A(n,n,true);
      mathx::array<double> x(n,1);
      mathx::array<double> b = mathx::linsolv::product(A,x);

      mathx::array<double> x_1 = mathx::linsolv::least_squares_QR(A,b);
      if(debug){
        std::cout << "A" << std::endl;
        std::cout << A.to_string() << std::endl;
        std::cout << "b: " << b.to_string() << std::endl;
        std::cout << "x: " << x.to_string() << std::endl;
        std::cout << "x_1: " << x_1.to_string() << std::endl;
      }

      std::cout << n << "," << mathx::vectors::euclideanLength(x - x_1) << std::endl;
    }
}

void problem6(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem Six                   " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  int m = 21;
  mathx::array<double> tt(21);
  mathx::array<double> bb(21);
  tt[0] = 0;
  double pi = 3.1415926535897;
  for(int i = 0; i < 21; i++){
    bb[i] = std::cos(2 * pi * tt[i]);
    tt[i + 1] = tt[i] + 1./ (m - 1.);
  }

  mathx::matrix<double> x_i(21,4,1.);
  for(int i = 0; i < x_i.rows(); i++){
    for(int j = 0; j < x_i.cols(); j++){
      x_i[i][j] = std::pow(tt[i], j + 1);
    }
  }

  mathx::array<double> x = mathx::linsolv::least_squares(x_i,bb);

  std::cout << "x: ";
  x.to_string();
  std::cout << std::endl;
}

void problem7(bool debug){
    std::cout << "----------------------------------------------" << std::endl;
    std::cout << "                Problem Seven                 " << std::endl;
    std::cout << "----------------------------------------------" << std::endl;
    mathx::matrix<double> x_i(11, 1, 1.);

    if(debug) std::cout << "x_i: " << x_i.to_string() << std::endl;
    mathx::array<double> y_i = {0.9,1.01,1.05,0.97,0.98,0.95,0.01,-0.1,0.02,-0.1,0.0};
    if(debug) std::cout << "y_i: " << y_i.to_string() << std::endl;
    mathx::array<double> x = mathx::linsolv::least_squares(x_i,y_i);
    std::cout << "x: ";
    x.to_string();
    std::cout << std::endl;
}

void problem8(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem Eight                 " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  double pi = 3.1415926535897;
  auto f = [pi](double t){ return .05 * std::sin(1000 * t) + .5 * std::cos(pi * t) - .4 * std::sin(10 * t);};
  mathx::array<double> y(101);
  mathx::array<double> t(101, 0);
  std::cout << "t,f(t)" << std::endl;
  for(int i = 0; i < 101; i++){
    y[i] = f(t[i]);
    std::cout << t[i] << "," << y[i] << std::endl;
    t[i + 1] = t[i] + .01;
  }

  mathx::matrix<double> x_i(101,5,1.);
  for(int i = 1; i < x_i.rows(); i++){
    for(int j = 0; j < x_i.cols(); j++){
      x_i[i][j] = std::pow(t[i], j + 1);
    }
  }

  mathx::array<double> x = mathx::linsolv::least_squares(x_i,y);

  std::cout << "x: ";
  x.to_string();
  std::cout << std::endl;
}

void problem9(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem Nine                  " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  auto f = [](double t){ return -11. + (55./3.)*t-(17./2.)*(t*t)+(7./6.)*(t*t*t); };
  mathx::array<double> y(33);
  mathx::array<double> t(33, 0);
  t[0] = .9;
  std::cout << "t,f(t)" << std::endl;
  for(int i = 0; i < 33; i++){
    y[i] = f(t[i]);
    std::cout << t[i] << "," << y[i] << std::endl;
    t[i + 1] = t[i] + .1;
  }

  mathx::matrix<double> x_i(33,4,1.);
  for(int i = 1; i < x_i.rows(); i++){
    for(int j = 0; j< x_i.cols(); j++){
      x_i[i][j] = std::pow(t[i], j + 1);
    }
  }

  mathx::array<double> x = mathx::linsolv::least_squares(x_i,y);
  std::cout << "x: ";
  x.to_string();
  std::cout << std::endl;
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
