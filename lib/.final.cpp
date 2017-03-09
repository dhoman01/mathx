#include "mathx.hpp"
#include "goodrand.hpp"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <chrono>
#include <tuple>

void problem2(bool debug=false);
void problem3(bool debug=false);
void problem4(bool debug=false);
void problem5(bool debug=false);

mathx::matrix<double> generateRandom(int n, bool debug = false);
mathx::matrix<double> generateSPD(int n);
void doPrintRunningTime(std::chrono::duration<double> child_runningtime);

int main(){
  problem2();
  problem3();
  problem4();
  problem5();


  return EXIT_SUCCESS;
}

void problem2(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem Two                   " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;

  // Test on matrices that produce
  // small (k(I) = 1) and large
  // condition numbers
  {
    mathx::matrix<double> I(3,3);
    for(int i = 0; i < 3; i++)
      for(int j = 0; j < 3; j++)
        I[i][j] = (i == j ? 1 : 0);

    double id_kappa = mathx::linsolv::kappa(I);
    std::cout << "Small condition number: " << id_kappa << std::endl; // k(I) = 1

    mathx::matrix<double> A = {{10,20},{10000.0001, 20}};
    double kappa = mathx::linsolv::kappa(A);
    std::cout << "Large condition number: " << kappa << std::endl; // k(A) ~ 5 x 10^6
  }

  std::cout << "\n\n\n" << std::endl;

  for(int n = 16; n <= 256; n *= 2){
    // Random matrix of size n x n
    mathx::matrix<double> A(n, n, true);

    // Find condition numbers of A
    double one_kappa = mathx::linsolv::kappa(A);
    double inf_kappa = mathx::linsolv::kappa(A, 1);

    std::cout << "n = " << n << std::endl;
    std::cout << "One-Condition Number: " << one_kappa << std::endl;
    std::cout << "Infinity-Condition Number: " << inf_kappa << std::endl;
  }
}

void problem3(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem Three                 " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  // Function to use is (x-1)^2 - 3
  auto f = [](double x){ return std::pow(x-1,2) - 3; };

  // Time the hybrid method on large interval
  auto startHybrid = std::chrono::steady_clock::now();
  double root = mathx::roots::hybrid_method(f, 0, 4, std::pow(10,-16), 10000);
  auto endHybrid = std::chrono::steady_clock::now();

  std::cout << "hybrid result = " << root << std::endl;
  doPrintRunningTime(endHybrid - startHybrid);

  // Time the secant method alone on same interval
  // expecting that it should fail as x0 and x1
  // are not "sufficiently" close to a root
  auto startSec = std::chrono::steady_clock::now();
  double root2 = mathx::roots::secant_method(f, 0, 4, std::pow(10,-16), 10000);
  auto endSec = std::chrono::steady_clock::now();
  std::cout << "secant result = " << root2 << std::endl;
  doPrintRunningTime(endSec - startSec);
}

void problem4(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem Four                  " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  int better = 0;
  for(int n = 16; n <= 256; n *= 2){
    mathx::matrix<double> A = generateRandom(n);
    mathx::array<double> x(n,1);
    mathx::array<double> b = mathx::linsolv::matmul(A,x);

    mathx::array<double> no_pivot = mathx::linsolv::solve(A, b, 0);
    mathx::array<double> pivot = mathx::linsolv::solve(A, b, 2);

    if(debug){
      std::cout << "A" << std::endl;
      std::cout << A.to_string() << std::endl;
    }

    double no_pivot_error = mathx::vectors::norm(x - no_pivot);
    double pivot_error = mathx::vectors::norm(x - pivot);
    std::cout << "Error with no pivot:                " << no_pivot_error << std::endl;
    std::cout << "Error with scaled-partial pivoting: " << pivot_error << std::endl;
    if(pivot_error > no_pivot_error) better++;
  }

  std::cout << "Pivoting was better " << better << " times." << std::endl;
}

void problem5(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem Five                  " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  for(int n = 16; n <= 256; n *= 2){
    mathx::matrix<double> A = generateSPD(n);
    mathx::array<double> x(n,1);
    mathx::array<double> b = mathx::linsolv::matmul(A, x);

    mathx::array<double> xstar = mathx::linsolv::least_squares_QR(A, b);

    if(debug){
      std::cout << "A" << std::endl;
      std::cout << A.to_string() << std::endl;
    }

    double error = mathx::vectors::norm(x - xstar);
    std::cout << "Error: " << error << std::endl;
  }
}

mathx::matrix<double> generateRandom(int n, bool debug){
  mathx::matrix<double> A = mathx::matrix<double>(n,n);
  if(debug) std::cout << "Generating a random matrix of dim(" << A.rows() << ", " << A.cols() << ")" << std::endl;
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      A[i][j] = goodrand::get_rand(-5.0, 5.0);

  return A;
}

/* Method to generate a matrix
 * that is s.p.d.
 */
mathx::matrix<double> generateSPD(int n){
  mathx::matrix<double> A = mathx::matrix<double>(n, n, true);

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
