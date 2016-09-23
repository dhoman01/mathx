#include <sstream>
#include "gtest/gtest.h"
#include "mathx.hpp"

TEST(ComplexStructTest, Constructor) {
  mathx::complex<double> w(3.0001, 4.0001);
  EXPECT_EQ(3.0001, w.real);
  EXPECT_EQ(4.0001, w.imaginary);
}

TEST(ComplexStructTest, OperatorPlus) {
  mathx::complex<double> w(3.0001, 4.0001);
  EXPECT_EQ(mathx::complex<double>(3.0001 + 3.0001, 4.0001 + 4.0001), w + w);
}

TEST(ComplexStructTest, OperatorMinus) {
  mathx::complex<double> w(3.0001, 4.0001);
  EXPECT_EQ(mathx::complex<double>(0, 0), w - w);
}

TEST(ComplexStructTest, OperatorOutstream) {
  mathx::complex<double> w(6, 8);
  std::stringstream ss;
  ss << w;
  EXPECT_EQ("6+8i", ss.str());
}
