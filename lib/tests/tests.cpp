#include <cmath>
#include <iostream>
#include <limits>
#include "gtest/gtest.h"
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
