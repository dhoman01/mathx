#ifndef VECTORS_HPP
#define VECTORS_HPP

#include <cmath>
#include <exception>
#include "array.hpp"

namespace mathx {
namespace vectors {

template <typename T>
T dot_product(array<T> v, array<T> w) {
  if (v.size() != w.size()) throw std::runtime_error("Vector dot products are only defined for vectors of the same length");
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
T norm(array<T> v) {
  return std::sqrt(dot_product(v, v));
}

template <typename T>
T norm(T* v, int size){
  return std::sqrt(dot_product(v,v,size));
}

template<typename T>
array<T> normalize(array<T> v){
  array<T> normal = v;
  T norm = norm(v);
  for(int i = 0; i < normal.size(); i++){
    normal[i] /= norm;
  }

  return normal;
}

template <typename T>
array<T> cross_product(array<T> v, array<T> w) {
  if (v.size() != w.size() && v.size() != 3) throw std::runtime_error("Vector cross product is only defined for vectors of length 3");
  array<T> vxw = {v[1] * w[2] - v[2] * w[1], v[2] * w[0] - v[0] * w[2],
                        v[0] * w[1] - v[1] * w[0]};
  return vxw;
}

template <typename T>
T one_norm(array<T> v) {
  T norm = 0;

  for (T e : v)
    norm += std::abs(e);

  return norm;
}

template <typename T>
T max_norm(array<T> v) {
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
