#ifndef UTILS_HPP
#define UTILS_HPP

#include <cmath>
#include <cstdlib>
#include <exception>
#include <typeinfo>
#include "mathx.hpp"

namespace mathx {
/*! The utils namespace contains useful utility methods that are seperated by
 * namespaces. */
namespace utils {

// Return the absolute value of a real number
template <typename T>
T absolute_value(T x) {
  return x >= 0 ? x : -1 * x;
}

// Return the absolute value of a complex number
template <typename T>
double absolute_value(mathx::complex<T> x) {
  return std::sqrt((x.real * x.real) + (x.imaginary * x.imaginary));
}

/*! The precision namespace contains useful methods dealing with precision */
namespace precision {

/**
* @brief This method calculates the machine's epsilon based on the data type T.
* @details The method calculates the smallest value such that when added
* to 1 it is not equal to 1.
* Example of use:
* \code{.cpp}
* #include "mathx.hpp"
*
* double dbl_eps = mathx::utils::precision::machine_epsilon<double>();
* float flt_eps = mathx::utils::precision::machine_epsilon<float>();
*
* std::cout << "Double precision=" << dbl_eps << std::endl;
* std::cout << "Single precision=" << flt_eps << std::endl;
* // Output:
* // Double precision=2.220446e-16
* // Single precision=1.192093e-07
* \endcode
*/
template <typename T>
T machine_epsilon() {
  int pow = 0;
  T eps = 1;
  while (eps + 1 != 1) {
    eps /= 2;
    --pow;
  }

  return std::pow(2, pow + 1);
};
}

/*! The error namespace contains useful methods for calulating error. */
namespace error {

// clang-format off
/**
* @brief Calculates the absolute error of an approximation of \f$f'(x)\f$ using
* a one-sided difference formula.
* @details A one-sided difference formula approximates the derivative in the
* following manner:
* \f[f'(x)=\frac{f(x+h)-f(x)}{h}+\frac{f''(\xi)h^2}{2}+\cdots\f]
* This leads to the following absolute error for the approximation:
* \f[e_{abs}=\bigg|f'(x)-\frac{f(x+h)-f(x)}{h}\bigg|\leq\frac{f''(\xi)h^2}{2}\Rightarrow
* e_{abs}\leq Ch\f]
* Where \f$\xi\f$ is some unknown constant. If \f$\xi\f$ was known then the
* absolute error would be known exactly.
* This approximation is \f$O(h)\f$.
* Example of use:
* \code{.cpp}
* #include "mathx.hpp"
*
* auto f = [](double x){ return std::exp(-2*x); };
* auto df = [](double x){ return -2*std::exp(-2*x); };
* auto e_abs = mathx::utils::error::one_sided_difference(f, df, 0.5, std::pow(10,-8));
*
* std::cout << "Absolute error: " << std::endl;
* // Output: Absolute error: 7.447261e-10
* \endcode
* @param f A function of the type double(double)
* @param df The derivative of \f$f\f$. Also of the type double(double). If the
* derivative is not known (or easily calculable) pass a function returning the
* value considered as correct for \f$f'(x)\f$.
* @param x The value to evaluate \f$f\f$ and \f$f'\f$ at
* @param h The amount to add to x. h is typically expressed as 1.e-8 and
* \f$0<h<1\f$
*/
// clang-format off
double one_sided_difference(function *f, function *df, double x, double h) {
  return absolute_value(df(x) - (f(x + h) - f(x)) / h);
}
// clang-format off
/**
* @brief Calculates the absolute error of an approximation of \f$f'(x)\f$ using a
* central-difference formula.
* @details A central-difference formula approximates the derivative in the
* following manner:
* \f[f'(x)=\frac{f(x+h)-f(x-h)}{2h}+\frac{f^{(3)}(\xi)h^3}{12}+\cdots\f]
* This leads to the following absolute error for the approximation:
* \f[e_{abs}=\bigg|f'(x)-\frac{f(x+h)-f(x-h)}{2h}\bigg|\leq\frac{f^{(3)}(\xi)h^3}{12}\Rightarrow
* e_{abs}\leq Ch^2\f]
* Where \f$\xi\f$ is some unknown constant. If \f$\xi\f$ was known then the
* absolute error would be known exactly.
* This approximation is \f$O(h^2)\f$.
* Example of use:
* \code{.cpp}
* #include "mathx.hpp"
*
* auto f = [](double x){ return std::exp(-2*x); };
* auto df = [](double x){ return -2*std::exp(-2*x); };
* auto e_abs = mathx::utils::error::central_difference(f, df, 0.5, std::pow(10,-8));
*
* std::cout << "Absolute error: " << std::endl;
* // Output: Absolute error: 2.030832e-09
* \endcode
* @param f A function of the type double(double)
* @param df The derivative of \f$f\f$. Also of the type double(double). If the
* derivative is not known (or easily calculable) pass a function returning the
* value considered as correct for \f$f'(x)\f$.
* @param x The value to evaluate \f$f\f$ and \f$f'\f$ at
* @param h The amount to add to x. h is typically expressed as 1.e-8 and
* \f$0<h<1\f$
*/
// clang-format on
double central_difference(function *f, function *df, double x, double h) {
  if (h == 0) throw divide_by_zero;
  return absolute_value(df(x) - (f(x + h) - f(x - h)) / (2 * h));
}

// clang-format off
/**
* @brief Calculates the absolute error in the approximation of one number  of type T by another.
* @details Absolute error of two numbers of type T is defined as:
* \f[|x-x_0|<\epsilon_{abs}, \textrm{ where $x$ and $y$ are type T}
* @param x A number approximation of type T and is greater than 0
* @param x0 An approximation of x
*/
// clang-format on
template <typename T>
T e_abs(T x, T x0) {
  return absolute_value(x - x0);
}

template <typename T>
T e_abs(mathx::complex<T> x, mathx::complex<T> x0) {
  return absolute_value(x - x0);
}

template <typename T>
T e_rel(T x, T x0) {
  if (x == 0) throw divide_by_zero;
  return absolute_value(x - x0) / absolute_value(x);
}

template <typename T>
T e_rel(mathx::complex<T> x, mathx::complex<T> x0) {
  if (x.real == 0 && x.imaginary == 0) throw divide_by_zero;
  return absolute_value(x - x0) / absolute_value(x);
}
}
}
}

#endif
