#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include <sstream>
#include <stdexcept>
#include <iostream>

namespace mathx {

/*! Polynomial interpolation is very useful in real-world applications of mathematics. This is because there are often times when the desired function is unknown, but data can be gathered to map inputs to outputs. Interpolation allows us to take a set of points \f$(x_i,y_i)\f$ and interpolate the generating polynomial. The method used in this namespace is finding the Newton's form coefficients from a divided differences table. */
namespace interpolation {

/**
* @brief Compute the divided difference table for points \f$(x_i,y_i)\f$
* @details The divided difference table entries are defined as \f{eqnarray*}{f[x_i]&=&f(x_i)\\f[x_i,\cdots,x_j]&=&\frac{f[x_{i+1},\cdots,x_j]-f[x_i,\cdots,x_{j-1}]}{x_j-x_i}.\f}
* @param x - An array of x points in a cartesian pair
* @param f - An array of y points in a cartesian pair, where y=f(x)
* @returns diff_table - a matrix<T> storing the divied difference table of the input data
*/
template<typename T>
matrix<T> divided_differences(array<T>& x, array<T>& f){
  // Matrix representing the difference
  // table. NOTE: diff[i][j] j > i will
  // be junk data as it is never assigned
  matrix<T> diff(x.size(), f.size());

  // Assign the first column of the
  // table to be equal to f(x_i)
  for(int j = 0; j < f.size(); j++)
    diff[j][0] = f[j];

  // Compute the rest of the columns
  // of the table using the previous
  // entries
  for(int i = 1; i < x.size(); i++){
    for(int j = 1; j <= i; j++){
      diff[i][j] = (diff[i][j-1] - diff[i-1][j-1]) / (x[i] - x[i - j]);
    }
  }

  // Return the computed table
  return diff;
}

/**
* @brief Return the Newton's Form coefficients from a divied differences table
* @details The diagonal of a divided differences table are the coefficients of the Newton's Form polynomial of the function that produced the inputs for the table.
* @param diff_table - A divided differences table
* @returns coeff - an array<T> that contains the Newton's form coefficients of the input table
*/
template<typename T>
array<T> newtons_coeff(matrix<T>& diff_table){
  // Extract the diagonal
  array<T> coeff;
  for(int i = 0; i < diff_table.rows(); i++)
    coeff.push(diff_table[i][i]);

  // Return the coefficients
  return coeff;
}

/**
* @brief evaluates a Newton's form polynomial at \f$x\f$ for the given coefficients
*/
template<typename T>
T eval_newtons(double x, array<T>& xi, array<T>& coeff){
  int np1 = xi.size();
  T p = coeff[np1 - 1];
  for(int j = np1 - 2; j >= 0; j--){
    p = p * (x - xi[j]) + coeff[j];
  }

  return p;
}

}

/** @example interpolation.cpp
* This example shows how to use the methods of the interpolation namespace.
*/

}

#endif
