#include "mathx.hpp"
#include "goodrand.hpp"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <chrono>
#include <tuple>


mathx::matrix<double> generateSPD(int n, bool debug = false);

void doPrintRunningTime(std::chrono::duration<double> child_runningtime);

int main(){
  // mathx::matrix<double> A = {{7,3,1},{-3,10,2},{1,7,-15}};
  // mathx::array<double> b = {3,4,2};
  mathx::matrix<double> A = generateSPD(10);
  mathx::array<double> x(10,1);
  mathx::array<double> b = mathx::linsolv::product(A,x);
  mathx::array<double> x0(10, 5);
  mathx::array<double> xstar = mathx::linsolv::cgm(A,b,x0,std::pow(10,-8),1000);
  std::cout << "error " << mathx::vectors::euclideanLength(x - xstar) << std::endl;
  std::cout << "xstar ";
  xstar.to_string();
  std::cout << std::endl;
  return EXIT_SUCCESS;
}


mathx::matrix<double> generateSPD(int n, bool debug){
  mathx::matrix<double>A = mathx::matrix<double>(n,n, true);
  for(int i = 0; i < n; i++){
    for(int j = i + 1; j < n; j++){
      A[i][j] = A[j][i];
    }
  }

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
