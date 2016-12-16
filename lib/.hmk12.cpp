#include "mathx.hpp"
#include "goodrand.hpp"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <chrono>
#include <tuple>

void problem2(bool debug=false);
void problem3(bool debug=false);

void doPrintRunningTime(std::chrono::duration<double> child_runningtime);

int main(){
  problem2();
  problem3();

  // mathx::matrix<double> A = {{1,1,1,1,1,1},
  //                            {0,(1./5.),(2./5.),(3./5.),(4./5.),1},
  //                            {0,(1./25.),(4./25.),(9./25.),(16./25.),1},
  //                            {0,(1./125.),(8./125.),(27./125.),(64./125.),1},
  //                            {0,(1./625.),(16./625.),(81./625.),(256./625.),1},
  //                            {0,(1./3125.),(32./3125.),(243./3125.),(1024./3125.),1}};
  // mathx::array<double> b = {1,(1./2.),(1./3.),(1./4.),(1./5.),(1./6.)};

  mathx::matrix<double> A = {{1,1,1,1,1,1},{-1,-0.6,-0.2,0.2,0.6,1},{1,.36,.04,.04,.36,1},
                             {-1,-0.216,-0.008,0.008,0.216,1}, {1,0.1296,0.0016,0.0016,0.1296,1},
                             {-1,-0.07776,-0.00032,0.00032,0.07776,1}};
  mathx::array<double> b = {2,0,(2./3.),0,(2./5.),0};

  std::cout << "A" << std::endl;
  std::cout << A.to_string() << std::endl;
  std::cout << "b: " << b.to_string() << std::endl;

  mathx::array<double> coef = mathx::linsolv::solve(A,b,0);
  std::cout << "coef: " << std::endl;
  for(int i = 0; i < coef.size(); i++){
    std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << coef[i] << std::endl;
  }

  return EXIT_SUCCESS;
}

void problem2(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem Two                   " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  mathx::matrix<double> A = {{1.,1.,1.,1.,1.},{2.,1.,0,-1.,-2.},{2.,.5,0,.5,2.},{(4./3.),(1./6.),0,-(1./6.),-(4./3.)},{(2./3.),(1./24.),0,(1./24.),(2./3.)}};
  mathx::array<double> b = {0., 0., 0., 1., 0.};


  mathx::array<double> coeff = mathx::linsolv::solve(A,b,2);
  std::cout << "coeff: " << coeff.to_string() << std::endl;

  if(debug){
    std::cout << "A" << std::endl;
    std::cout << A.to_string() << std::endl;
    std::cout << "b: " << b.to_string() << std::endl;
  }

  /**********************************************/
  /***           TEST OF APPROX               ***/
  /**********************************************/
  auto f = [](double x){ return std::sin(x); };
  auto d3f = [](double x){ return -1 * std::cos(x); };
  auto approx = [&](double x, double h){
    double num =  f(x + (2 * h)) - ( 2 * f(x + h) ) + ( 2 * f(x - h) ) - f(x - (2 * h));
    double denom = 2 * std::pow(h, 3);
    return num / denom;
   };

  auto x = 4.;
  auto h = std::pow(10,-4);
  std::cout << "\n\n------- TEST --------" << std::endl;
  std::cout << "x       = " << x << std::endl;
  std::cout << "h       = " << h << std::endl;
  std::cout << "f(x)    = " << f(x) << std::endl;
  std::cout << "f'''(x) = " << d3f(x) << std::endl;
  std::cout << "f'''(x) ~ " << approx(x, h) << std::endl;
}

void problem3(bool debug){
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "                Problem Three                 " << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  const double PI = 3.14159265358979323846;
  auto f = [](double t){ return 3 * std::pow(t, 2) - 2 * std::pow(t, 3); };
  auto g = [PI](double x){ return std::sin(PI * x); };
  mathx::array<double> t;
  mathx::array<double> ft;
  mathx::array<double> x;
  mathx::array<double> fx;

  t.push(goodrand::get_rand(-1.0, 0.0));
  x.push(goodrand::get_rand(0.0, 1.0));
  ft.push(t[0]);
  fx.push(x[0]);

  for(double i = 1; i < 16; i++){
    double ti = goodrand::get_rand(t[i - 1], 2.0);
    double xi = goodrand::get_rand(t[i - 1], 3.0);
    t.push(ti);
    ft.push(f(ti));
    x.push(xi);
    fx.push(g(xi));
  }

  mathx::matrix<double> diff_table_t = mathx::interpolation::divided_differences(t, ft);
  mathx::matrix<double> diff_table_x = mathx::interpolation::divided_differences(x, fx);

  mathx::array<double> coeff_t = mathx::interpolation::newtons_coeff(diff_table_t);
  mathx::array<double> coeff_x = mathx::interpolation::newtons_coeff(diff_table_x);

  std::cout << "coeff t: " << coeff_t.to_string() << std::endl;
  std::cout << "coeff x: " << coeff_x.to_string() << std::endl;

  if(debug){
    std::cout << "t: " << t.to_string() << std::endl;
    std::cout << "x: " << x.to_string() << std::endl;
    std::cout << "diff_table_t" << std::endl;
    std::cout << diff_table_t.to_string() << std::endl;
    std::cout << "diff_table_x" << std::endl;
    std::cout << diff_table_x.to_string() << std::endl;
  }
}
