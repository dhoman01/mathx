#ifndef VECTORS_HPP
#define VECTORS_HPP

#include <cmath>
#include <exception>
#include "array.hpp"

namespace mathx {
namespace vectors {
class not_equal_exception : public std::exception {
  virtual const char *what() const throw() {
    return "Your vectors need to be the same dimensions to perform this "
           "operation!";
  }
} not_equal_exception;

template <typename T>
T dotProduct(array<T> v, array<T> w) {
  if (v.size() != w.size()) throw not_equal_exception;
  T product = 0;
  for (int i = 0; i < v.size(); i++) {
    product += v[i] * w[i];
  }

  return product;
}

template<typename T>
T dot_product(T* v, T* w, int size){
  T product = 0;
  for(int i = 0; i < size; i++){
    product += v[i] * w[i];
  }

  return product;
}

template <typename T>
T euclideanLength(array<T> v) {
  return std::sqrt(dotProduct(v, v));
}

template <typename T>
T euclideanLength(T* v, int size){
  return std::sqrt(dot_product(v,v,size));
}

template <typename T>
array<T> crossProduct(array<T> v, array<T> w) {
  if (v.size() != w.size() && v.size() != 3) throw not_equal_exception;
  array<T> vxw = {v[1] * w[2] - v[2] * w[1], v[2] * w[0] - v[0] * w[2],
                        v[0] * w[1] - v[1] * w[0]};
  return vxw;
}

template <typename T>
T oneNorm(array<T> v) {
  T norm = 0;

  for (T e : v)
    norm += std::abs(e);

  return norm;
}

template <typename T>
T maxNorm(array<T> v) {
  T max = 0;
  for (int i = 0; i < v.size(); i++) {
    T x = std::abs(v[i]);
    max = x > max ? x : max;
  }

  return max;
}
}
}

#endif
