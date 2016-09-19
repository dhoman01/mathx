#ifndef MATHX_HPP
#define MATHX_HPP

#include <cmath>
#include <iostream>

namespace mathx {
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
};
}
#include "integration.hpp"
#include "utils.hpp"

#endif
