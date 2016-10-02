#include "gtest/gtest.h"
#include "array.hpp"

using namespace mathx::roots;

TEST(RootsTest, BisectionMethodTest){
  auto f = [](double x){ return (x*x)+(3*x)+2; };
  double root = bisect(f, -1.5, 0.0, -.25, 2.0, std::pow(10,-16), 50);
  EXPECT_DOUBLE_EQ(-1.00, root);
}

TEST(RootsTest, FixedPointIterationTest){
  auto f = [](double x){ return (x*x)+(3*x)+2; };
  auto g = [](double x){ return -1*((x*x+2)/3);};
  double root = fixed_point_iter(g, f, -.75, std::pow(10,-20), 100);
  EXPECT_DOUBLE_EQ(-1.000000000000000000001, root);
}

TEST(RootsTest, NewtonsMethodTest){
  auto f = [](double x){ return (x*x)+(3*x)+2; };
  auto df = [](double x){ return 2*x + 3; };
  double root = mathx::roots::newtons_method(f, df, 1, std::pow(10, -8), 50);
  EXPECT_DOUBLE_EQ(-1.000000000000000000001, root);
}

TEST(RootsTest, SecantMethodTest){
  auto f = [](double x){ return (x*x)+(3*x)+2; };
  double root = secant_method(f, -1.5, -.75, std::pow(10, -16), 50);
  EXPECT_DOUBLE_EQ(-1.000000000000000000001, root);
}

TEST(RootsTest, HybridMethodTest){
  auto f = [](double x){ return (x*x)+(3*x)+2; };
  auto df = [](double x){ return 2*x + 3; };
  double root = hybrid_method(f, df, -1.999999, 3, std::pow(10,-16), 50);
  EXPECT_DOUBLE_EQ(-1.000000000000000000001, root);
}
