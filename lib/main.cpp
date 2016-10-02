#include <cstdlib>  // Gets machines EXIT_SUCCESS
#include "mathx.hpp"

int main(){
  auto f = [](double x){ return (x*x)+(3*x)+2; };
  auto g = [](double x){ return -1*((x*x+2)/3); };
  auto df = [](double x){ return 2*x + 3; };
  double root = mathx::roots::bisect(f, -1.5, 0.0, -.25, 2.0, std::pow(10,-8), 50);
  std::cout << "root of X^2+3x+2 " << root << std::endl;
  std::cout << root << "^2 + 3" << root << " + 2 = " << f(root) << std::endl;
  std::cout << std::endl;
  double root2 = mathx::roots::fixed_point_iter(g, f, -.75, std::pow(10,-8), 50);
  std::cout << "root2 " << root2 << std::endl;
  std::cout << std::endl;
  double root3 = mathx::roots::newtons_method(f, df, -3, std::pow(10, -8), 50);
  std::cout << "root3 " << root3 << std::endl;
  std::cout << std::endl;
  double root4 = mathx::roots::secant_method(f, -1.5, -.75, std::pow(10, -8), 50);
  std::cout << "root4 " << root4 << std::endl;
  return EXIT_SUCCESS;
}
