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

void doPrintRunningTime(std::chrono::duration<double> child_runningtime);

int main(){
  problem1();
  problem2();
  problem3();

  return EXIT_SUCCESS;
}

void problem1(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem One                   " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  const double PI = 3.14159265358979323846;
  auto f = [](double t){ return 3 * std::pow(t, 2) - 2 * std::pow(t, 3); };
  auto g = [PI](double x){ return std::sin(PI * x); };
  for(int n = 4; n <= 32; n *= 2){
    mathx::array<double> t;
    mathx::array<double> ft;
    mathx::array<double> x;
    mathx::array<double> fx;
    double spacing = 3.0 / (n - 1);
    for(double i = -1.0; i <= 2.0; i += spacing){
      t.push(i);
      ft.push(f(i));
    }

    for(double i = 0.0; i <= 3.0; i += spacing){
      x.push(i);
      fx.push(g(i));
    }

    mathx::matrix<double> diff_table_t = mathx::interpolation::divided_differences(t, ft);
    mathx::matrix<double> diff_table_x = mathx::interpolation::divided_differences(x, fx);

    std::cout << "diff_table_t" << std::endl;
    std::cout << diff_table_t.to_string() << std::endl;
    std::cout << "diff_table_x" << std::endl;
    std::cout << diff_table_x.to_string() << std::endl;

    if(debug){
      std::cout << "t: " << t.to_string() << std::endl;
      std::cout << "x: " << x.to_string() << std::endl;
    }
  }
}

void problem2(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem Two                   " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  const double PI = 3.14159265358979323846;
  auto f = [](double t){ return 3 * std::pow(t, 2) - 2 * std::pow(t, 3); };
  auto g = [PI](double x){ return std::sin(PI * x); };
  for(int n = 4; n <= 32; n *= 2){
    mathx::array<double> t;
    mathx::array<double> ft;
    mathx::array<double> x;
    mathx::array<double> fx;
    double spacing = 3.0 / (n - 1);
    for(double i = -1.0; i <= 2.0; i += spacing){
      t.push(i);
      ft.push(f(i));
    }

    for(double i = 0.0; i <= 3.0; i += spacing){
      x.push(i);
      fx.push(g(i));
    }

    mathx::matrix<double> diff_table_t = mathx::interpolation::divided_differences(t, ft);
    mathx::matrix<double> diff_table_x = mathx::interpolation::divided_differences(x, fx);

    mathx::array<double> coeff_t = mathx::interpolation::newtons_coeff(diff_table_t);
    mathx::array<double> coeff_x = mathx::interpolation::newtons_coeff(diff_table_x);

    std::cout << n << ": coeff t: " << coeff_t.to_string() << std::endl;
    std::cout << n << ": coeff x: " << coeff_x.to_string() << std::endl;

    if(debug){
      std::cout << "t: " << t.to_string() << std::endl;
      std::cout << "x: " << x.to_string() << std::endl;
      std::cout << "diff_table_t" << std::endl;
      std::cout << diff_table_t.to_string() << std::endl;
      std::cout << "diff_table_x" << std::endl;
      std::cout << diff_table_x.to_string() << std::endl;
    }
  }
}

void problem3(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem Three                 " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  const double PI = 3.14159265358979323846;
  auto f = [](double t){ return 3 * std::pow(t, 2) - 2 * std::pow(t, 3); };
  auto g = [PI](double x){ return std::sin(PI * x); };
  mathx::array<double> t;
  mathx::array<double> ft;
  mathx::array<double> x;
  mathx::array<double> fx;

  t.push(goodrand::get_rand(-1.0, 0.0));
  x.push(goodrand::get_rand(0.0, 1.0));
  ft.push(t[0]);
  fx.push(x[0]);

  for(double i = 1; i < 16; i++){
    double ti = goodrand::get_rand(t[i - 1], 2.0);
    double xi = goodrand::get_rand(t[i - 1], 3.0);
    t.push(ti);
    ft.push(f(ti));
    x.push(xi);
    fx.push(g(xi));
  }

  mathx::matrix<double> diff_table_t = mathx::interpolation::divided_differences(t, ft);
  mathx::matrix<double> diff_table_x = mathx::interpolation::divided_differences(x, fx);

  mathx::array<double> coeff_t = mathx::interpolation::newtons_coeff(diff_table_t);
  mathx::array<double> coeff_x = mathx::interpolation::newtons_coeff(diff_table_x);

  std::cout << "coeff t: " << coeff_t.to_string() << std::endl;
  std::cout << "coeff x: " << coeff_x.to_string() << std::endl;

  if(debug){
    std::cout << "t: " << t.to_string() << std::endl;
    std::cout << "x: " << x.to_string() << std::endl;
    std::cout << "diff_table_t" << std::endl;
    std::cout << diff_table_t.to_string() << std::endl;
    std::cout << "diff_table_x" << std::endl;
    std::cout << diff_table_x.to_string() << std::endl;
  }
}
