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
* @returns \f$\varepsilon\f$ - the machine precision
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

/** @example precision.cpp
* This example shows how to find machine precision for both single and double.
*/

/*! The error namespace contains useful methods for calulating error. */
namespace error {


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
* @param f A function of the type double(double)
* @param df The derivative of \f$f\f$. Also of the type double(double). If the
* derivative is not known (or easily calculable) pass a function returning the
* value considered as correct for \f$f'(x)\f$.
* @param x The value to evaluate \f$f\f$ and \f$f'\f$ at
* @param h The amount to add to x. h is typically expressed as 1.e-8 and
* \f$0<h<1\f$
* @returns e - the error of a one sided difference approximation
*/

double one_sided_difference(function *f, function *df, double x, double h) {
  return absolute_value(df(x) - (f(x + h) - f(x)) / h);
}

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
* @param f A function of the type double(double)
* @param df The derivative of \f$f\f$. Also of the type double(double). If the
* derivative is not known (or easily calculable) pass a function returning the
* value considered as correct for \f$f'(x)\f$.
* @param x The value to evaluate \f$f\f$ and \f$f'\f$ at
* @param h The amount to add to x. h is typically expressed as 1.e-8 and
* \f$0<h<1\f$
* @returns e - the error of a central difference approximation
*/

double central_difference(function *f, function *df, double x, double h) {
  if (h == 0) throw divide_by_zero;
  return absolute_value(df(x) - (f(x + h) - f(x - h)) / (2 * h));
}


/**
* @brief Calculates the absolute error in the approximation of one number  of type T by another.
* @details Absolute error of two numbers of type T is defined as:
* \f[|x-x_0|<\epsilon_{abs}, \textrm{ where $x$ and $y$ are type T}
* @param x A number approximation of type T and is greater than 0
* @param x0 An approximation of x
* @returns e - the absolute error of the two inputs
*/
template <typename T>
T e_abs(T x, T x0) {
  return absolute_value(x - x0);
}

/**
* @brief Calculates the absolute error in the approximation of one number  of type mathx::complex<T> by another.
* @details Absolute error of two numbers of type mathx::complex<T> is defined as:
* \f[|x-x_0|<\epsilon_{abs}, \textrm{ where $x$ and $y$ are type mathx::complex<T>}
* @param x A number approximation of type mathx::complex<T> and is greater than 0
* @param x0 An approximation of x
* @returns e - the absolute error of the two inputs
*/
template <typename T>
T e_abs(mathx::complex<T> x, mathx::complex<T> x0) {
  return absolute_value(x - x0);
}

/**
* @brief Calculates the relative error in the approximation of one number  of type T by another.
* @details Relative error of two numbers of type T is defined as:
* \f[\frac{|x-x_0|}{|x|}<\epsilon_{rel}, \textrm{ where $x$ and $y$ are type T}
* @param x A number approximation of type T and is greater than 0
* @param x0 An approximation of x
* @returns e - the relative error of the two inputs
*/
template <typename T>
T e_rel(T x, T x0) {
  if (x == 0) throw divide_by_zero;
  return absolute_value(x - x0) / absolute_value(x);
}

/**
* @brief Calculates the relative error in the approximation of one number  of type mathx::complex<T> by another.
* @details Relative error of two numbers of type mathx::complex<T> is defined as:
* \f[\frac{|x-x_0|}{|x|}<\epsilon_{rel}, \textrm{ where $x$ and $y$ are type mathx::complex<T>}
* @param x A number approximation of type mathx::complex<T> and is greater than 0
* @param x0 An approximation of x
* @returns e - the relative error of the two inputs
*/
template <typename T>
T e_rel(mathx::complex<T> x, mathx::complex<T> x0) {
  if (x.real == 0 && x.imaginary == 0) throw divide_by_zero;
  return absolute_value(x - x0) / absolute_value(x);
}

}

namespace poly {

/**
* @brief This method evaluates a polynomial using the nested method.
* @details The nested method of evaluating a polynomial is as follows:
* \f[p_n(x)=(\cdots((c_nx+c_{n-1})x + c_{n-2})x \cdots)x + c_0 \f]
* @param coeff - an array of the coefficients of the polynomial. \f$p(x)=c_0+c_1x+c_2x^2+\cdots+c_nx^2 \f$
* @param x - the value to evaluate \f$p\f$ at
* @returns f(x) - the evaluation of \f$f\f$ at \f$x\f$
*/
template<typename T>
T nested_eval(array<T> coeff, T x){
  // Initialize p
  T p = coeff[coeff.size() - 1];

  // Loop through the coefficients "unwrapping"
  // the evaluation
  for(int i = coeff.size() - 2; i >= 0; i--)
    p = p*x + coeff[i];

  // Return answer
  return p;
}

}

}

}

#endif
