/** @file main.cpp */
/** @mainpage
* Use Mathx to solve your computantional mathematic's problems.
* @author Dustin E. Homan
* @date 09/05/2016
* @version 0.1
*/
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include "mathx.hpp"

template <typename T>
void printElement(char align, T t, const int& width, char fill) {
  std::cout << (align == 'L' ? std::left : std::right) << std::setw(width)
            << std::setfill(fill) << t;
}

void problemOne();
void problemTwo();
void problemThree();
void problemFive();
void problemSix();
double round(double, int);
std::vector<double> round(std::vector<double>, int);
std::string vectorToString(std::vector<double>);

/**
* @brief This is the main driver for Mathx.
* @details Use this to solve your problems. It includes all the relevant
* headers. Referer to the documentation to determine what fundtion you need.
*/
int main() {
  problemOne();
  // problemTwo();
  // problemThree();
  // problemFive();
  // problemSix();
  return EXIT_SUCCESS;
}

void problemOne() {
  mathx::complex<double> w(3.0001, 4.0001);
  mathx::complex<double> z(2.9999, 3.9999);
  double u[5] = {1, 1, -1.5, 100, 100};
  double v[5] = {0.99, 1.01, -1.2, 99.99, 99};
  printElement('R', "u", 6, ' ');
  printElement('R', "v", 6, ' ');
  printElement('R', "Absolute Error", 15, ' ');
  printElement('R', "Relative Error", 15, ' ');
  std::cout << std::endl;
  printElement('R', "-", 45, '-');
  std::cout << std::endl;
  for (int i = 0; i < 5; i++) {
    printElement('R', u[i], 6, ' ');
    printElement('R', v[i], 6, ' ');
    printElement('R', mathx::utils::error::e_abs(u[i], v[i]), 15, ' ');
    printElement('R', mathx::utils::error::e_rel(u[i], v[i]), 15, ' ');
    std::cout << std::endl;
  }

  std::cout << std::endl;

  std::cout << "w=" << w << std::endl;
  std::cout << "z=" << z << std::endl;
  std::cout << "|w|=" << mathx::utils::absolute_value(w) << std::endl;
  std::cout << "w+z=" << w + z << std::endl;
  std::cout << "w-z=" << w - z << std::endl;
  std::cout << "e_abs>" << mathx::utils::error::e_abs(w, z) << std::endl;
  std::cout << "e_rel>" << mathx::utils::error::e_rel(w, z) << std::endl;
}

void problemTwo() {
  std::vector<double> one = {4, 5, 6};
  std::vector<double> two = {3.33, 7, 8};
  std::vector<double> three = {.5, .33, .25};
  std::cout << "\n\nVector one: " << vectorToString(one) << std::endl;
  std::cout << "Vector two: " << vectorToString(two) << std::endl;
  std::cout << "Vector three: " << vectorToString(three) << std::endl;
  std::cout << std::endl;
  std::cout << "Euclidean Lengths:" << std::endl;
  std::cout << "One: " << mathx::vectors::euclideanLength(one) << std::endl;
  std::cout << "Two: " << mathx::vectors::euclideanLength(two) << std::endl;
  std::cout << "Three: " << mathx::vectors::euclideanLength(three) << std::endl;
  std::cout << std::endl;
  std::cout << "One norms:" << std::endl;
  std::cout << "One: " << mathx::vectors::oneNorm(one) << std::endl;
  std::cout << "Two: " << mathx::vectors::oneNorm(two) << std::endl;
  std::cout << "Three: " << mathx::vectors::oneNorm(three) << std::endl;
  std::cout << std::endl;
  std::cout << "Max norms:" << std::endl;
  std::cout << "One: " << mathx::vectors::maxNorm(one) << std::endl;
  std::cout << "Two: " << mathx::vectors::maxNorm(two) << std::endl;
  std::cout << "Three: " << mathx::vectors::maxNorm(three) << std::endl;
}

void problemThree() {
  std::vector<double> one = {4, 5, 6};
  std::vector<double> two = {3.33, 7, 8};
  std::vector<double> three = {.5, .33, .25};
  std::cout << "\n\nVector one: " << vectorToString(one) << std::endl;
  std::cout << "Vector two: " << vectorToString(two) << std::endl;
  std::cout << "Vector three: " << vectorToString(three) << std::endl;
  std::cout << std::endl;
  std::cout << "Dot Products:" << std::endl;
  std::cout << "One dot Two: " << mathx::vectors::dotProduct(one, two)
            << std::endl;
  std::cout << "One dot Three: " << mathx::vectors::dotProduct(one, three)
            << std::endl;
  std::cout << "Two dot Three: " << mathx::vectors::dotProduct(two, three)
            << std::endl;
  std::cout << std::endl;
  std::cout << "Cross Products:" << std::endl;
  std::cout << "One cross Two: "
            << vectorToString(mathx::vectors::crossProduct(one, two))
            << std::endl;
  std::cout << "One cross Three: "
            << vectorToString(mathx::vectors::crossProduct(one, three))
            << std::endl;
  std::cout << "Two cross One: "
            << vectorToString(mathx::vectors::crossProduct(two, one))
            << std::endl;
  std::cout << "Two cross Three: "
            << vectorToString(mathx::vectors::crossProduct(two, three))
            << std::endl;
  std::cout << "Three cross One: "
            << vectorToString(mathx::vectors::crossProduct(three, one))
            << std::endl;
  std::cout << "Three cross two: "
            << vectorToString(mathx::vectors::crossProduct(three, two))
            << std::endl;
}

std::string vectorToString(std::vector<double> v) {
  std::stringstream ss;
  ss << "( ";
  for (int i = 0; i < v.size(); i++) {
    ss << v[i] << (i == v.size() - 1 ? " )" : ", ");
  }
  return ss.str();
}

void problemFive() {
  // Approximation of f'(x) where f(x)=sin(x);
  auto df = [](double x, double h) {
    return (2.0 * std::cos((2.0 * x + h) / 2.0) * std::sin(h / 2.0)) / h;
  };

  // Output column headers
  std::cout << std::right << std::setw(30) << "h";
  std::cout << std::right << std::setw(30) << "approximation";
  std::cout << std::right << std::setw(30) << "absolute error" << std::endl;

  double actual = std::cos(1.2);
  // Compute the approximation of f'(1.2) for h=1e-20,1e-19,...,1
  for (int h_exp = -20; h_exp < 1; h_exp++) {
    double h = std::pow(10, h_exp);
    double approx = df(1.2, h);
    std::cout << std::setprecision(20) << std::right << std::setw(30) << h;
    std::cout << std::setprecision(20) << std::right << std::setw(30) << approx;
    std::cout << std::setprecision(20) << std::right << std::setw(30)
              << mathx::utils::error::e_abs(approx, actual) << std::endl;
  }
}

void problemSix() {
  double value = 6.55638;
  double digits = 3;
  std::vector<double> values = {3.352346, 2.2346234, 4.3426462, 42362.2436,
                                4.2436234};
  std::cout << "Value:   " << value << std::endl;
  std::cout << "Rounded: " << round(value, digits) << std::endl;
  std::cout << "Values:  " << vectorToString(values) << std::endl;
  std::cout << "Rounded: " << vectorToString(round(values, digits))
            << std::endl;
}

double round(double value, int digits) {
  double shift = std::pow(10, digits);
  return std::floor(value * shift + 0.5) / shift;
}

std::vector<double> round(std::vector<double> values, int digits) {
  std::for_each(values.begin(), values.end(),
                [&](auto& value) { value = round(value, digits); });

  return values;
}
