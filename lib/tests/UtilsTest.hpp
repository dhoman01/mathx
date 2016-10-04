#include <cmath>
#include <iostream>
#include <limits>
#include "gtest/gtest.h"
#include "mathx.hpp"
#include "utils.hpp"

TEST(MachineEpsilonTest, CorrectDoublePrecision) {
  EXPECT_EQ(std::numeric_limits<double>::epsilon(),
            mathx::utils::precision::machine_epsilon<double>());
}

TEST(MachineEpsilonTest, CorrectSinglePrecision) {
  EXPECT_EQ(std::numeric_limits<float>::epsilon(),
            mathx::utils::precision::machine_epsilon<float>());
}

TEST(DerivativeAbsoluteErrorTest, OneSidedDifference) {
  auto f = [](double x) { return std::sin(x); };
  auto df = [](double x) { return std::cos(x); };
  auto x = 1.2;
  auto h = std::pow(10, -8);
  auto expected = std::abs(df(x) - (f(x + h) - f(x)) / h);
  EXPECT_EQ(expected, mathx::utils::error::one_sided_difference(f, df, x, h));
}

TEST(DerivativeAbsoluteErrorTest, CentralDifference) {
  auto f = [](double x) { return std::sin(x); };
  auto df = [](double x) { return std::cos(x); };
  auto x = 1.2;
  auto h = std::pow(10, -8);
  auto expected = std::abs(df(x) - (f(x + h) - f(x - h)) / (2 * h));
  EXPECT_EQ(expected, mathx::utils::error::central_difference(f, df, x, h));
}

TEST(AbsoluteValueTest, PositiveIn) {
  EXPECT_EQ(2, mathx::utils::absolute_value(2));
}

TEST(AbsoluteValueTest, NegativeIn) {
  EXPECT_EQ(2, mathx::utils::absolute_value(-2));
}

TEST(AbsoluteValueTest, ComplexIn) {
  mathx::complex<double> w(3.0001, 4.0001);
  EXPECT_NEAR(5.00014, mathx::utils::absolute_value(w), 0.0000001);
}

TEST(AbsoluteErrorTest, RealIn) {
  EXPECT_NEAR(0.01, mathx::utils::error::e_abs(1.0, 0.99), 0.00001);
}

TEST(AbsoluteErrorTest, ComplexIn) {
  mathx::complex<double> w(3.0001, 4.0001);
  mathx::complex<double> z(2.9999, 3.9999);
  EXPECT_NEAR(0.000282843, mathx::utils::error::e_abs(w, z), 0.00001);
}

TEST(RelateiveErrorTest, RealIn) {
  EXPECT_NEAR(0.01, mathx::utils::error::e_rel(1.0, 0.99), 0.00001);
}

TEST(RelateiveErrorTest, ComplexIn) {
  mathx::complex<double> w(3.0001, 4.0001);
  mathx::complex<double> z(2.9999, 3.9999);
  EXPECT_NEAR(0.000056567, mathx::utils::error::e_rel(w, z), 0.00001);
}

TEST(PolyTest, NestedEvalTest){
  mathx::array<double> coeff = {-512, 2304, -4608, 5376, -4032, 2016, -672, 144, -18, 1};
  double x = 5;
  double p = mathx::utils::poly::nested_eval(coeff, x);
  EXPECT_NEAR(19683, p, 0.00001);
}
