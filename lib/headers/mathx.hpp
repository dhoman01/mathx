#ifndef MATHX_HPP
#define MATHX_HPP

#include <cmath>
#include <iostream>
#include "vectors.hpp"
#include "array.hpp"
#include "matrix.hpp"
#include "linsolv.hpp"
#include "interpolation.hpp"

/*! @mainpage Introduction
* \section Introduction
* Mathx is a software package developed to solve several problems in the realm of computational linear algebra. This package was developed as part of the MATH 5610 course at Utah State University, Fall 2016. The contents of this manual include the methods used to solve problems such as finding roots of polynomials, solving linear systems, using linear systems to solve the least-squares problem, and function approximation.\n\n
* Mathx is organized based on namespace. Each namespace includes multiple functions related to that namespace. Descriptions of the namespaces can be found at the beginning of each detailed documentation of that namespace. The source code is included inline after the signature, brief, description, and parameter listing of each function.
*/
namespace mathx {

/**
* @brief Defines the type of a function that has one double parameter and
* returns a double.
*/
typedef double function(double);

class divide_by_zero : public std::exception {
  virtual const char* what() const throw() { return "Cannot divide by zero!"; }
} divide_by_zero;

template <typename T>
struct complex {
  complex<T>(T r, T i) : real(r), imaginary(i){};
  T real;
  T imaginary;

  complex operator+(const complex& rhs) {
    return complex(real + rhs.real, imaginary + rhs.imaginary);
  }
  complex operator-(const complex& rhs) {
    return complex(real - rhs.real, imaginary - std::abs(rhs.imaginary));
  }

  friend std::ostream& operator<<(std::ostream& os, const complex& obj) {
    os << obj.real << (obj.imaginary > 0 ? "+" : "") << obj.imaginary << "i";
    return os;
  }

  friend std::istream& operator>>(std::istream& is, const complex<T>& obj) {
    is >> obj.real >> obj.imaginary;
    return is;
  }

  friend bool operator==(const complex<T>& lhs, const complex<T>& rhs) {
    return lhs.real == rhs.real && lhs.imaginary == rhs.imaginary;
  }
};
}

#include "utils.hpp"
#include "roots.hpp"

#endif
