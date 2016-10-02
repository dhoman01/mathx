#include "MathxTest.hpp"
#include "UtilsTest.hpp"
#include "ArrayTest.hpp"
#include "VectorsTest.hpp"
#include "RootsTest.hpp"
#include "gtest/gtest.h"

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
