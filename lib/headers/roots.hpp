#ifndef ROOTS_HPP
#define ROOTS_HPP

#include <cmath>
#include <stdexcept>

namespace mathx {
  /*! This namespace contians functions used to find roots of functions.\n\n
  *   Beginning with the bisection method, bisection is a robust algorithm as the only
  *   condition for convergence is that the function changes signs on the interval \f$[a,b]\f$. Bisection is a \f$O(n)\f$ method in converenge.\n\n
  *   Fixed point iteration (or functional iteration) does not need a bracketing interval. However, it does need a guess that is "sufficiently" close enough to the root and has some restrictions on \f$f\f$. For fixed point iteration to converge there must exist some \f$g\in C[a,b]\f$ with \f$a\leq g(x)\leq b\f$ for all \f$x\in [a,b]\f$. This implies that there is a fixed point \f$x^{*}\f$ in  the interval \f$[a,b]\f$. If, in addition, the derivative of \f$g\f$ exists and there is a constant \f$\rho<1\f$ such that \f[|g'(x)|\leq \rho\quad\quad\textrm{for all }x\in [a,b],\f] then the fixed point is unique @cite AscherGrief Fixed Point iteration's rate of convergence is dependant on the choice of \f$g(x)\f$.\n\n
  *   Newton's method restricts \f$f\f$ to the set \f$C^2[a,b]\f$ and defines \f[x_{k+1}=x_k-\frac{f(x_k)}{f'(x_k)}.\f]. Newton's method is \f$O(n^2)\f$. Some downsides of Newton's methods are that the derivative must exist and you must know how to evaluate it and you must also know the local nature of the method's convergence @cite AscherGrief\n\n
  *   The Secant method addresses the disadvantage of Newton's method requiring the derivative. This is accomplished by replacing \f$f'\f$ with its finite difference approximation \f[f'(x_k)\approx \frac{f(x_k)-f(x_{k-1})}{x_k-x_{k-1}}.\f]. The secant method converges superlinearly.\n\n
  *   As both Newton's method and the Secant method require a guess that is "sufficiently" close to the root, where "sufficiently" is problem dependant, a naive approach would be to combine two methods. The hybrid method, or globalization of the Secant method, attempts to find a root on the interval \f$[a,b]\f$ using the bisection method. The method will run for a few iterations thus reducing the size of the interval. After a few iterations, the Secant method will be attempted, using \f$a\f$ and \f$b\f$ from bisection as the initial guesses, and will continue as long as it is making improvement that is better than bisection (\f$|f_{k+1}|<0.5|f_k|\f$). This globalization of the Secant method alleviates the need to have knowledge of the local nature of the method's convergence. A similar method could be written for globalizing Newton's method.
  */
  namespace roots {
    /**
    * @brief Finds roots of a polynomial in range [a,b] using bisection.
    * @details The bisection method requires a continuous function on the interval
    * of [a,b], with \f$f(a)\cdot f(b)<0\f$. Bisection will find the root by taking
    * the midpoint, c, between a and b and evaluating if \f$f(a)\cdot f(c)<0\f$ OR
    * \f$f(b)\cdot f(c)<0\f$. With c becoming the new b or a, respectivally.
    * Linear error convergence.
    * @param f A continuous function. Data type double(double). See mathx::typedef function.
    * @param a - a point at which \f$f(a) < 0\f$
    * @param b - a point at which \f$f(b) > 0\f$ and \f$a<b\f$
    * @param fa - \f$f(a)\f$. Used for convience to reduce the number of evaluations
    * @param fb - \f$f(b)\f$. Used for convience to reduce the number of evaluations
    * @param tol - the tolerance between a and b. Must be smaller than \f$b_{initial}-a_{initial}\f$
    * @param max - maximum iterations the algorithm is allowed to execute. Used to prevent infinite loops
    * @returns r - a root of the function \f$f\f$
    */
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
    * @returns r - a root of the function \f$f\f$
    */
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

    /**
    * @brief Find a polynomial's roots using Netwon's method.
    * @details Newton's method is an iterative method used to
    * find roots. Newton's method requires that \f$f\in C^2[a,b]\f$.
    * Also the \f$k+1\f$ iterate is calculated by:
    * \f[x_k-\frac{f(x_k)}{f'(x_k)}, k=0,1,2,\cdots\f]
    * @param f - a function in \f$C^2[a,b]\f$. Data type double(double). See mathx::typedef function.
    * @param df - the derivative of f. Data type double(double). See mathx::typedef function.
    * @param x0 - an initial guess of the root
    * @param tol - the tolerance between a and b. Must be greater than 0
    * @param max - maximum iterations the algorithm is allowed to execute. Used to prevent infinite loops
    * @returns r - a root of the function \f$f\f$
    */
    double newtons_method(function f, function df, double x0, double tol, uint max){
      // Check if guess is root
      if(f(x0) == 0) return x0;

      // Ensure valid problem
      if(tol < 0 || df(x0) == 0) throw std::runtime_error("check your parameters");

      // Initialize xk to be used in error calculation
      double xk = x0 - 1;
      // Work loop
      for(uint k = 0; k < max && std::abs(x0 - xk) > tol; k++){
        xk = x0;                // Store x_k

        double fk = f(xk);
        double dfk = df(xk);

        // If xk is a root return
        if(fk == 0) return xk;

        // If the derivative equals 0 retrun NaN
        if(dfk == 0) return std::nan("0");

        x0 = xk - fk / dfk; // Store x_{k+1}
      }

      // Either the tolerance or max
      // iterations were exceeded, return x0
      // If |x0-xk| < tol then x0 is an
      // approximation of the root
      return x0;
    }

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
    * @returns r - a root of the function \f$f\f$
    */
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
        x0 = x0 -( ( fk * ( x0 - x1 ) ) / ( fk - fk_1 ) );    // Caclulate x_{k+1}
        x1 = tmp;                                             // Store x_k
        fk_1 = fk;
        fk = f(x0);
      }

      // Either the tolerance or max
      // iterations were exceeded, return x0
      // If |x0-x1| < tol then x0 is an
      // approximation of the root
      return x0;
    }

    /**
    * @brief A globalization of the Secant method, using bisection to get a better guess.
    * @details The Secant method requires two initial guesses sufficiently close to the root.
    * Whereas bisection just requires that \f$f\f$ changes signs on \f$[a,b]\f$. Thus we can
    * use bisection to arrive at a better guesses for the Secant method.
    * @param f - a function in \f$C^2[a,b]\f$. Data type double(double). See mathx::typedef function.
    * @param a - a point at which \f$f(a) < 0\f$
    * @param b - a point at which \f$f(b) > 0\f$ and \f$a<b\f$
    * @param tol - the tolerance between a and b. Must be greater than 0
    * @param max - maximum iterations the algorithm is allowed to execute. Used to prevent infinite loops
    * @returns r - a root of the function \f$f\f$
    */
    double hybrid_method(function f, double a, double b, double tol, uint max){
      // Ensure parameters lead to valid problem
      if(a == b || tol <= 0) throw std::runtime_error("check your parameters");

      // Initialize variables
      double fa = f(a);
      double fb = f(b);

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

      double c = std::nan("1");

      // Work loop
      // Do five steps of bisection
      // Then try Secant method until
      // one of the following:
      //   1. Max iter reached
      //   2. Tolerance exceeded
      //   3. Convergence is not optimal (|fkp1|<0.5|fk|)
      for(uint k = 0; k <= max && std::abs(b - a) > tol; k++){
        c = (a + b) / 2;   // Assign c to be midpoint of a and b
        double fc = f(c);

        // Reinitialize variables
        if(fa * fc < 0){   // If fa*fc < 0 then the root is in [a,c]
          b = c;           // c becomes new b
          fb = fc;
        } else {           // else fb*fc < 0 and the root is in [c,b]
          a = c;           // c becomes new a
          fa = fc;
        }

        // Try the Secant method on 5
        // iteration intervals
        if(k % 5 == 0 && k != 0){
          double x0 = a;
          double x1 = b;
          // Initialize f_k and f_{k+1}
          double fk = f(x0);
          double fk_1 = f(x1);

          // test if either guess is root
          if(fk == 0) return x0;
          if(fk_1 == 0) return x1;

          // Work loop
          for(uint n = 0; n < max && std::abs(x0 - x1) > tol && std::abs(fk) < 0.5 * std::abs(fk_1); n++){
            double tmp = x0;
            x0 = x0 -( (fk * ( x0 - x1 ) ) / ( fk - fk_1 ) );    // Caclulate x_{k+1}
            x1 = tmp;                                            // Store x_k
            fk_1 = fk;
            fk = f(x0);
          }


          // The tolerance was exceeded
          if(std::abs(x0 - x1) < tol)
            return x0;
        }
      }

      // Either the tolerance or max
      // iterations were exceeded, return c
      // If |b-a| < tol then c is an
      // approximation of the root
      return c;
    }
  }

  /** @example roots.cpp
  * This example shows how to find the roots of \f$f(x)=(x-1)^2-3\f$ using all five methods.
  */
}

#endif
