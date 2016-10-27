#include <iostream>
#include <iomanip>

namespace mathx {
  namespace linsolv {

    /**
    * @brief Perform back substitution to solve Ux=b
    * @details
    * @param - U an upper triangular matrix
    * @param - b a vector of values for the right-hand side of the equation
    */
    template<typename T>
    array<T> back_substitution(matrix<T>& U, array<T>& b){
      // Initialize solution vector
      array<T> x(b.size(), 0);
      // Back Substitution
      for(int k = b.size() - 1; k >= 0; k--){
        x[k] = b[k];
        for(int j = k+1; j < U.cols(); j++)
          x[k] -= U.get(k,j) * x[j];
        x[k] /= U.get(k, k);
      }

      // Return the solution
      return x;
    };

    /**
    * @brief Perform forward substitution to solve Lx=b
    * @details
    * @param - L a lower triangular matrix
    * @param - b a vector of values for the right-hand side of the equation
    */
    template<typename T>
    array<T> foward_subsitition(matrix<T>& L, array<T>& b){
      // Initialize solution vector
      array<T> x(b.size(), 0);

      // Forward Substitution
      for(int k = 0; k < b.size(); k++){
        x[k] = b[k];
        for(int j = 0; j < b.size() - 1; j++)
          x[k] -= L.get(k, j) * x[k];
        x[k] /= L.get(k, k);
      }

      // Return the solution
      return x;
    }

    /**
    * @brief Perform Guassian Elimination on a matrix given a solution vector
    * @details
    * @param A - the input matrix
    * @param b - the solution vector
    */
    template<typename T>
    void gaussian_elimination(matrix<T>& A, array<T>& b){
      for(int k = 0; k < A.rows() - 1; k++){
        for(int i = k + 1; i < A.rows(); i++){
          T l = A.get(i, k) / A.get(k, k);
          for(int j = k + 1; j < A.rows(); j++){
            T a_kj = A.get(i, j) - l * A.get(k, j);
            A.set(i, j, a_kj);
          }
          b[i] = b[i] - l * b[k];
        }
      }
    }

    /**
    * @brief Solve the linear system Ax=b
    * @details
    * @param A - input matrix
    * @param b - solution vector
    */
    template<typename T>
    array<T> solve(matrix<T> A, array<T> b){
      gaussian_elimination(A, b);
      return back_substitution(A, b);
    }

    /**
    * @brief Multiply a matrix by a vector
    * @details
    * @param A - input matrix
    * @param x - input vector
    */
    template<typename T>
    array<T> product(matrix<T>& A, array<T>& x){
      array<T> b(x.size(), 0);
      for(int i = 0; i < A.rows(); i++){
        for(int j = 0; j < A.cols(); j++){
          b[i] += A.get(i, j) * x[j];
        }
      }

      return b;
    }
  }
}
