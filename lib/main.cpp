/** @file main.cpp */
/** @mainpage
* Use Mathx to solve your computantional mathematic's problems.
* @author Dustin E. Homan
* @date 09/05/2016
* @version 0.1
*/
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include "mathx.hpp"

/**
* @brief This is the main driver for Mathx.
* @details Use this to solve your problems. It includes all the relevant
* headers. Referer to the documentation to determine what fundtion you need.
*/

std::string problem_one();
std::string problem_three();
double problem_five(int, double);
double rec_evaluation(double, int, int, double);

int main() {
  auto one = problem_one();
  auto three = problem_three();
  auto five = problem_five(20, std::pow(10, -5));

  // std::cout << "Problem One\n-----------------" << std::endl;
  // std::cout << one << std::endl;
  // std::cout << "\n\nProblem Three\n----------------" << std::endl;
  // std::cout << three << std::endl;
  std::cout << "\n\nProblem Five\n----------------" << std::endl;
  std::cout << five << std::endl;

  return EXIT_SUCCESS;
}

std::string problem_one() {
  std::stringstream ss;
  auto f = [](double x) { return std::exp(-2 * x); };
  auto df = [](double x) { return -2 * std::exp(-2 * x); };
  auto x = 0.5;

  ss << "h, Absolute Error(one-sided diff.)" << std::endl;
  for (int h_exp = 0; h_exp >= -16; h_exp--) {
    auto h = std::pow(10, h_exp);
    ss << std::scientific << h << ", "
       << mathx::utils::error::one_sided_difference(f, df, x, h) << std::endl;
  }

  return ss.str();
}

std::string problem_three() {
  std::stringstream ss;
  auto f = [](double x) { return std::sin(x); };
  auto df = [](double x) { return std::cos(x); };
  auto x = 1.2;

  ss << "h, Absolute Error(centeral-diff.)" << std::endl;
  for (int h_exp = 0; h_exp >= -16; h_exp--) {
    auto h = std::pow(10, h_exp);
    ss << std::scientific << h << ", "
       << mathx::utils::error::central_difference(f, df, x, h) << std::endl;
  }

  return ss.str();
}

double problem_five(int n, double des) {
  int n1 = n + std::abs(std::log10(des)) + 1;
  int e = 0;
  double eval = rec_evaluation(0.0, n1, e, des);
  return eval;
}

double rec_evaluation(double y_n, int n, int e, double des) {
  if (std::pow(.1, e) <= des) return y_n;
  return rec_evaluation(.1 * ((1.0 / n - (y_n))), --n, ++e, des);
}
