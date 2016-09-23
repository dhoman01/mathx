#include <cmath>
#include <iostream>
#include <limits>
#include "gtest/gtest.h"
#include "vectors.hpp"

TEST(DotProductTest, CorrectDotProduct) {
  std::vector<double> v = {4, 5, 6};
  std::vector<double> w = {3.33, 7, 8};
  EXPECT_EQ(96.32, mathx::vectors::dotProduct(v, w));
  EXPECT_EQ(mathx::vectors::dotProduct(v, w), mathx::vectors::dotProduct(w, v));
}

TEST(EuclideanLengthTest, CorrectLength) {
  std::vector<double> v = {4, 5, 6};
  EXPECT_EQ(8.774964387392123, mathx::vectors::euclideanLength(v));
}

TEST(CrossProductTest, CorrectCrossProduct) {
  std::vector<double> v = {-1, 7, 4};
  std::vector<double> w = {-5, 8, 4};
  std::vector<double> vxw = mathx::vectors::crossProduct(v, w);
  std::vector<double> wxv = mathx::vectors::crossProduct(w, v);
  double expected[3] = {-4, -16, 27};
  for (int i = 0; i < 3; i++) {
    EXPECT_EQ(expected[i], vxw[i]);
    EXPECT_EQ(-1 * expected[i], wxv[i]);
  }
}

TEST(OneNormTest, CorrectOneNorm) {
  std::vector<double> v = {4, 5, 6};
  EXPECT_EQ(15, mathx::vectors::oneNorm(v));
}

TEST(MaxNormTest, CorrectMaxNorm) {
  std::vector<double> v = {4, 5, 6};
  EXPECT_EQ(6, mathx::vectors::maxNorm(v));
}
