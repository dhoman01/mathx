#include "gtest/gtest.h"
#include "array.hpp"

using namespace mathx;

TEST(ArrayTest, ArrayInitializesWithSizeZero){
  array<int> arr;
  EXPECT_EQ(0, arr.size());
  EXPECT_EQ(0, arr.capacity());
}

TEST(ArrayTest, ArrayPushTest){
  array<int> arr;
  arr.push(1);
  EXPECT_EQ(1, arr.size());
  EXPECT_EQ(2, arr.capacity());
  EXPECT_EQ(1, arr[0]);
}

TEST(ArrayTest, ArrayGetTest){
  array<int> arr;
  arr.push(1);
  EXPECT_EQ(1, arr.get(0));
}

TEST(ArrayTest, ArrayPopTest){
  array<int> arr;
  arr.push(1);
  EXPECT_EQ(1, arr.pop());
  EXPECT_EQ(0, arr.size());
  EXPECT_EQ(0, arr.capacity());
}

TEST(ArrayTest, ArrayClearTest){
  array<int> arr;
  for(int i = 0; i < 100; i++)
    arr.push(i);

  arr.clear();
  EXPECT_EQ(0, arr.size());
  EXPECT_EQ(0, arr.capacity());
}

TEST(ArrayTest, ForEachLoopTest){
  array<double> arr;
  for(int i = 1; i <= 100; i++)
    arr.push(i);

  int sum = 0;
  for(double i : arr)
    sum += i;


  EXPECT_EQ(5050, sum);
}

TEST(ArrayTest, InitializerListTest){
  array<double> arr = {0,1,2};
  EXPECT_EQ(3, arr.size());
  EXPECT_EQ(3, arr.capacity());
  for(int i = 0; i < arr.size(); i++)
    EXPECT_EQ(i, arr[i]);
}
