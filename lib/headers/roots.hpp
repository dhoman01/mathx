#ifndef ROOTS_HPP
#define ROOTS_HPP

#include <cmath>
#include <stdexcept>

namespace mathx {
  namespace roots {
    double bisect(function f, double a, double b, double fa, double fb, double tol, uint max){
      if(a == b || fa * fb > 0 || tol <= 0) throw std::runtime_error("check your parameters");
      if(fa * fb == 0){
        if(fa == 0) return a;
        if(fb == 0) return a;
      }
      if(a > b) {
        double tmp = a;
        a = b;
        b = tmp;
        tmp = fa;
        fa = fb;
        fb = tmp;
      }

      double c = b;
      for(uint k = 0; k <= max && std::abs(b - a) > tol; k++){
        c = (a + b) / 2;
        double fc = f(c);
        if(fa * fc < 0){
          b = c;
          fb = fc;
        } else {
          a = c;
          fa = fc;
        }
      }

      return c;
    }

    double fixed_point_iter(function g, function f, double x0, double tol, uint max){
      if(tol < 0) throw std::runtime_error("check your parameters");
      double xk = x0 -1;
      for(uint k = 0; k < max && std::abs(x0 - xk) > tol; k++){
        xk = x0;
        x0 = g(x0);
      }

      return x0;
    }

    double newtons_method(function f, function df, double x0, double tol, uint max){
      if(f(x0) == 0) return x0;
      if(tol < 0 || df(x0) == 0) throw std::runtime_error("check your parameters");
      double xk = x0 -1;
      for(uint k = 0; k < max && std::abs(x0 - xk) > tol; k++){
        xk = x0;
        if(df(xk) == 0) throw std::runtime_error("division by zero");
        x0 = xk - f(xk)/df(xk);
      }

      return x0;
    }

    double secant_method(function f, double x0, double x1, double tol, uint max){
      double fk = f(x0);
      double fk_1 = f(x1);
      if(fk == 0) return x0;
      if(fk_1 == 0) return x1;
      if(tol < 0) throw std::runtime_error("check your parameters");
      for(uint k = 0; k < max && std::abs(x0 - x1) > tol; k++){
        double tmp = x0;
        x0 = x0 -((fk*(x0-x1))/(fk-fk_1));
        x1 = tmp;
        fk_1 = fk;
        fk = f(x0);
      }

      return x0;
    }

    double hybrid_method(function f, function df, double a, double b, double tol, uint max){
      double fa = f(a);
      double fb = f(b);
      double x0 = bisect(f, a, b, fa, fb, tol, 5);
      return newtons_method(f, df, x0, tol, max);
    }
  }
}

#endif
