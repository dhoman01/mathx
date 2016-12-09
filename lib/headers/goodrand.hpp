#ifndef GOOD_RAND_HPP
#define GOOD_RAND_HPP

#include <random>

namespace goodrand {
static std::random_device rd;
static std::mt19937 mt(rd());

// Get a random integer with
// uniform distribution in the
// range l to u
int get_rand(int l, int u) {
  std::uniform_int_distribution<> d(l, u);
  return d(mt);
}

// Get a random double with
// uniform distribution in the
// range l to u
double get_rand(double l, double u) {
  std::uniform_real_distribution<> d(l, u);
  return d(mt);
}
}

#endif
