#ifndef ROOTS_HPP
#define ROOTS_HPP

#include <cmath>
#include <stdexcept>

namespace mathx {
  namespace roots {
    // clang-format off
    /**
    * @brief Finds roots of a polynomial in range [a,b] using bisection.
    * @details The bisection method requires a continuous function on the interval
    * of [a,b], with \f$f(a)\cdot f(b)<0\f$. Bisection will find the root by taking
    * the midpoint, c, between a and b and evaluating if \f$f(a)\cdot f(c)<0\f$ OR
    * \f$f(b)\cdot f(c)<0$. With c becoming the new b or a, respectivally.
    * Linear error convergence.
    * @param f A continuous function. Data type double(double). See mathx::typedef function.
    * @param a - a point at which \f$f(a) < 0\f$
    * @param b - a point at which \f$f(b) > 0\f$ and $a<b$
    * @param fa - \f$f(a)\f$. Used for convience to reduce the number of evaluations
    * @param fb - \f$f(b)\f$. Used for convience to reduce the number of evaluations
    * @param tol - the tolerance between a and b. Must be smaller than \f$b_{initial}-a_{initial}\f$
    * @param max - maximum iterations the algorithm is allowed to execute. Used to prevent infinite loops
    */
    // clang-format on
    double bisect(function f, double a, double b, double fa, double fb, double tol, uint max){
      // Ensure parameters lead to valid problem
      if(a == b || fa * fb > 0 || tol <= 0) throw std::runtime_error("check your parameters");

      // Check if a or b are roots
      if(fa * fb == 0){
        if(fa == 0) return a;
        if(fb == 0) return a;
      }

      // If user inputs a > b
      // swap a and b, fa and fb
      if(a > b) {
        std::swap(a, b);
        std::swap(fa, fb);
      }

      // Work loop
      double c = std::nan("1");
      for(uint k = 0; k <= max && std::abs(b - a) > tol; k++){
        c = (a + b) / 2;   // Assign c to be midpoint of a and b
        double fc = f(c);
        if(fa * fc < 0){   // If fa*fc < 0 then the root is in [a,c]
          b = c;           // c becomes new b
          fb = fc;
        } else {           // else fb*fc < 0 and the root is in [c,b]
          a = c;           // c becomes new a
          fa = fc;
        }
      }

      // Either the tolerance or max
      // iterations were exceeded, return c
      // If |b-a| < tol then c is an
      // approximation of the root
      return c;
    }

    // clang-format off
    /**
    * @brief Find a polynomial's roots using fixed point iteration.
    * @details Fixed point iteration requires a continuous function
    * and an initial guess sufficiently close to the root. Using
    * \f$g(x)\f$, a function that \f$x\f$ is equal to, the algorithm
    * finds the root by evaluating a sequence of itearees.
    * Linear error convergence.
    * @param g - A function derived from f, such that \f$x=g(x)\f$. Data type double(double). See mathx::typedef function.
    * @param f - A continuous function. Data type double(double). See mathx::typedef function.
    * @param x0 - An initial guess of the root of f.
    * @param tol - the tolerance between a and b. Must be greater than 0
    * @param max - maximum iterations the algorithm is allowed to execute. Used to prevent infinite loops
    */
    // clang-format on
    double fixed_point_iter(function g, function f, double x0, double tol, uint max){
      // Ensure tolerance is valid
      if(tol < 0) throw std::runtime_error("check your parameters");


      double xk = x0 - 1; // Initialize xk to be used in calculating error
      // Work loop
      for(uint k = 0; k < max && std::abs(x0 - xk) > tol; k++){
        xk = x0;     // Store the x_k iteration
        x0 = g(x0);  // Caclulate x_{k+1}
      }

      // Either the tolerance or max
      // iterations were exceeded, return x0
      // If |x0-xk| < tol then x0 is an
      // approximation of the root
      if(std::abs(x0 - xk) > tol) return std::nan("1");
      return x0;
    }

    // clang-format off
    /**
    * @brief Find a polynomial's roots using Netwon's method.
    * @details Newton's method is an iterative method used to
    * find roots. Newton's method requires that \f$f\in C^2[a,b].
    * Also the \f$k+1\f$ iterate is calculated by:
    * \f[x_k-\frac{f(x_k)}{f'(x_k)}, k=0,1,2,\cdots\f]
    * @param f - a function in \f$C^2[a,b]\f$. Data type double(double). See mathx::typedef function.
    * @param df - the derivative of f. Data type double(double). See mathx::typedef function.
    * @param x0 - an initial guess of the root
    * @param tol - the tolerance between a and b. Must be greater than 0
    * @param max - maximum iterations the algorithm is allowed to execute. Used to prevent infinite loops
    */
    double newtons_method(function f, function df, double x0, double tol, uint max){
      // Check if guess is root
      if(f(x0) == 0) return x0;

      // Ensure valid problem
      if(tol < 0 || df(x0) == 0) throw std::runtime_error("check your parameters");

      // Initialize xk to be used in error calculation
      double xk = x0 -1;
      // Work loop
      for(uint k = 0; k < max && std::abs(x0 - xk) > tol; k++){
        xk = x0;                // Store x_k

        double fk = f(xk);
        double dfk = df(xk);

        // If xk is a root return
        if(fk == 0) return xk;

        // If the derivative equals 0 retrun NaN
        if(dfk == 0) return std::nan("0");

        x0 = xk - fk/dfk; // Store x_{k+1}
      }

      // Either the tolerance or max
      // iterations were exceeded, return x0
      // If |x0-xk| < tol then x0 is an
      // approximation of the root
      return x0;
    }

    // clang-format off
    /**
    * @brief Find a polynomial's roots using the Secant method (variant of Netwon's method)
    * @details A variation of Netwon's Method, the Secant method eliviates the need for
    * the user to provide the derivative of f. The Secant method does this by approximating
    * the \f$k+1\f$ iterate by:
    * \f[x_{k+1}=x_k-\frac{f(x_k)(x_k-x_{k-1})}{f(x_k)-f(x_{k-1})}, k=0,1,2,\cdots \f]
    * However, the user must provide to guesses of the root
    * @param f - a function in \f$C^2[a,b]\f$. Data type double(double). See mathx::typedef function.
    * @param x0 - an initial guess of the root
    * @param x1 - an second guess of the root
    * @param tol - the tolerance between a and b. Must be greater than 0
    * @param max - maximum iterations the algorithm is allowed to execute. Used to prevent infinite loops
    */
    // clang-format on
    double secant_method(function f, double x0, double x1, double tol, uint max){
      // Initialize f_k and f_{k+1}
      double fk = f(x0);
      double fk_1 = f(x1);

      // test if either guess is root
      if(fk == 0) return x0;
      if(fk_1 == 0) return x1;

      // Ensure tolerance is valid
      if(tol < 0) throw std::runtime_error("check your parameters");

      // Work loop
      for(uint k = 0; k < max && std::abs(x0 - x1) > tol; k++){
        double tmp = x0;
        x0 = x0 -((fk*(x0-x1))/(fk-fk_1));    // Caclulate x_{k+1}
        x1 = tmp;                             // Store x_k
        fk_1 = fk;
        fk = f(x0);
      }

      // Either the tolerance or max
      // iterations were exceeded, return x0
      // If |x0-x1| < tol then x0 is an
      // approximation of the root
      return x0;
    }

    // clang-format off
    /**
    * @brief A globalization of Newton's method, using bisection to get a better guess.
    * @details Newton's method requires an initial guess sufficiently close to the root.
    * Whereas bisection just requires that \f$f\f$ changes signs on [a,b]. Thus we can
    * use bisection to arrive at a better guess for Newton's method.
    * @param f - a function in \f$C^2[a,b]\f$. Data type double(double). See mathx::typedef function.
    * @param df - the derivative of f. Data type double(double). See mathx::typedef function.
    * @param a - a point at which \f$f(a) < 0\f$
    * @param b - a point at which \f$f(b) > 0\f$ and $a<b$
    * @param tol - the tolerance between a and b. Must be greater than 0
    * @param max - maximum iterations the algorithm is allowed to execute. Used to prevent infinite loops
    */
    double hybrid_method(function f, function df, double a, double b, double tol, uint max){
      // Initialize variables
      double fa = f(a);
      double fb = f(b);
      double x0 = std::nan("1");

      // Work loop
      uint count = 0;
      while(x0 != x0 && count < max){       // x0 != x0 is IEEE test for NaN
        // Use bisection to get a "good" guess
        x0 = bisect(f, a, b, fa, fb, tol, 5 + count);
        // Use Newton's method to find root
        // If guess isn't "good" enough fail and
        // use bisection again
        x0 = newtons_method(f, df, x0, tol, max / 2);
        count++;
      }

      // Either the root was found or
      // max iterations were reached.
      return x0;
    }
  }
}

#endif
