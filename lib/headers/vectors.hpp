#ifndef VECTORS_HPP
#define VECTORS_HPP

#include <cmath>
#include <exception>
#include "array.hpp"

namespace mathx {

/*! The vectors namespace contains useful vector utilities */
namespace vectors {

/**
* @brief Calculates the dot procuct of two vectors
* @param v - input vector
* @param w - input vector
* @returns s - the result of \f$<\textbf{v},\textbf{w}>\f$
*/
template <typename T>
T dot_product(array<T> v, array<T> w) {
  if (v.size() != w.size()) throw std::runtime_error("Vector dot products are only defined for vectors of the same length");
  T product = 0;
  for (int i = 0; i < v.size(); i++) {
    product += v[i] * w[i];
  }

  return product;
}

/**
* @brief Calculates the cross procuct of two vectors
* @param v - input vector
* @param w - input vector
* @returns u - an array<T> that is the result of \f$\textbf{v}\times\textbf{w}\f$
* @throws A std::runtime_error if length of either vector is not 3
*/
template <typename T>
array<T> cross_product(array<T> v, array<T> w) {
  if (v.size() != w.size() && v.size() != 3) throw std::runtime_error("Vector cross product is only defined for vectors of length 3");
  array<T> vxw = {v[1] * w[2] - v[2] * w[1], v[2] * w[0] - v[0] * w[2],
                        v[0] * w[1] - v[1] * w[0]};
  return vxw;
}

/**
* @brief Calculates the \f$l_2\f$-norm of a vector
* @param v - input vector
* @returns \f$||v||_2\f$ - the resulting norm
*/
template <typename T>
T norm(array<T> v) {
  return std::sqrt(dot_product(v, v));
}

/**
* @brief Calculates the \f$l_1\f$-norm of a vector
* @param v - input vector
* @returns \f$||v||_1\f$ - the resulting norm
*/
template <typename T>
T one_norm(array<T> v) {
  T norm = 0;

  for (T e : v)
    norm += std::abs(e);

  return norm;
}

/**
* @brief Calculates the \f$l_\infty\f$-norm of a vector
* @param v - input vector
* @returns \f$||v||_\infty\f$ - the resulting norm
*/
template <typename T>
T infinity_norm(array<T> v) {
  T max = 0;
  for (int i = 0; i < v.size(); i++) {
    T x = std::abs(v[i]);
    max = x > max ? x : max;
  }

  return max;
}

/**
* @brief Normalizes a vector
* @param v - input vector
* @returns v - an array<T> that is the result of the normalization
*/
template<typename T>
array<T> normalize(array<T> v){
  array<T> normal = v;
  T n = norm(v);
  for(int i = 0; i < normal.size(); i++){
    normal[i] /= n;
  }

  return normal;
}

}

/** @example vectors.cpp
* This example shows how to use the utility functions found in the ::vectors namespace.
*/

}

#endif
