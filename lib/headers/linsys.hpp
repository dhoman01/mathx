#include <cfloat>  // for std::isnan()
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <string>

namespace mathx {
  namespace linsolv {

    /**
    * @brief Multiply a matrix by a vector
    * @details
    * @param A - input matrix
    * @param x - input vector
    */
    template<typename T>
    array<T> product(matrix<T>& A, array<T>& x, bool a_trans = false){
      if(a_trans){
        array<T> b(A.cols(), 0);
        for(int j = 0; j < A.cols(); j++){
          for(int i = 0; i < A.rows(); i++){
            b[j] += A[i][j] * x[i];
          }
        }
        return b;
      } else {
        array<T> b(A.rows(), 0);
        for(int i = 0; i < A.rows(); i++){
          for(int j = 0; j < A.cols(); j++){
            b[i] += A[i][j] * x[j];
          }
        }
        return b;
      }
    }

    /**
    * @brief Multiply a tri-diagonal matrix by a vector
    * @details
    * @param A - input matrix
    * @param x - input vector
    */
    template<typename T>
    array<T> product(array<double> al, array<double> am, array<double>au, array<T>& x){
      array<T> b(x.size(), 0);
      b[0] = am[0] * x[0] + au[0] * x[1];
      for(int i = 1; i < x.size() - 1; i++){
        b[i] = (al[i] * x[i - 1]) + (am[i] * x[i]) + (au[i] * x[i + 1]);
      }

      b[x.size() - 1] = al[x.size() - 1] * x[x.size() - 2] + am[x.size() - 1] * x[x.size() - 1];

      return b;
    }

    /**
    * @brief Multiply two matrices
    * @details
    * @param A - input matrix
    * @param x - input vector
    */
    template<typename T>
    matrix<T> product(matrix<T>& A, matrix<T>& B){
      matrix<T> C(A.rows(), B.cols());
      for(int i = 0; i < A.rows(); i++){
        for(int j = 0; j < B.cols(); j++){
          for(int k = 0; k < A.cols(); k++){
            C[i][j] += A[i][k] * B[k][j];
          }
        }
      }

      return C;
    }

    /**
    * Returns transpose of a matrix
    * @param A - input matrix
    */
    template<typename T>
    matrix<T> transpose(matrix<T>& A){
      matrix<T> B(A.cols(), A.rows());
      for(int j = 0 ; j < A.cols(); j++){
        for(int i = 0; i < A.rows(); i++){
          B[j][i] = A[i][j];
        }
      }

      return B;
    }

    /**
    * @brief Multiply a matrix by its transpose (A^T)A
    * @details
    * @param A - input matrix
    **/
    template<typename T>
    matrix<T> mult_transpose(matrix<T>& A){
      matrix<T> B(A.cols(),A.cols(), (T)0);
      for(int i = 0; i < A.cols(); i++){
        for(int j = 0; j < A.cols(); j++){
          for(int k = 0; k < A.rows(); k++){
            B[i][j] += A[k][i]*A[k][j];
          }
        }
      }

      return B;
    }

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
          x[k] -= U[k][j] * x[j];
        x[k] /= U[k][k];
      }

      // Return the solution
      return x;
    };

    /**
    * @brief Perform forward substitution to solve Lx=b
    * @details
    * @param - L a lower triangular matrix
    * @param - b a vector of values for the right-hand side of the equation
    * @param - isLU a flag to interpret D as all ones
    */
    template<typename T>
    array<T> forward_substitution(matrix<T>& L, array<T>& b, bool isLU = false){
      // Initialize solution vector
      array<T> x(b.size(), 0);

      x[0] = b[0] / (isLU ? 1 : L[0][0]);
      // Forward Substitution
      for(int i = 1; i < b.size(); i++){
        x[i] = b[i];
        for(int j = 0; j < i; j++)
          x[i] -= L[i][j] * x[j];
        x[i] /= (isLU ? 1 : L[i][i]);
      }


      // Return the solution
      return x;
    }

    /**
    * @brief Perform Guassian Elimination on a matrix given a solution vector
    * @details
    * @param A - the input matrix
    * @param b - the solution vector
    * @param pstrategy - flag declaring the pivoting strategy
    *                    0 = no pivoting
    *                    1 = partial pivoting
    *                    2 = scaled partial pivoting
    */
    template<typename T>
    void gaussian_elimination(matrix<T>& A, array<T>& b, int pstrategy = 0){
      for(int k = 0; k < A.rows() - 1; k++){
        if(pstrategy > 0){
          int kpiv = pstrategy == 1 ? A.find_pivot(k) : A.find_scaled_pivot(k);
          A.swap(k, kpiv);
          std::swap(b[k], b[kpiv]);
        }

        for(int i = k + 1; i < A.rows(); i++){
          T l = A[i][k] / A[k][k];
          for(int j = k + 1; j < A.rows(); j++){
            T a_kj = A[i][j] - l * A[k][j];
            A[i][j] = a_kj;
          }
          b[i] = b[i] - l * b[k];
        }
      }
    }

    /**
    * @brief Factor a square matrix A into L and U
    * @details
    * @param A - input matrix
    * @param b - solution vector (used in pivoting)
    * @param pstrategy - flag declaring the pivoting strategy
    *                    0 = no pivoting
    *                    1 = partial pivoting
    *                    2 = scaled partial pivoting
    */
    template<typename T>
    void lu(matrix<T>& A, array<T>& b, int pstrategy = 0){
      for(int k = 0; k < A.rows() - 1; k++){
        if(pstrategy > 0){
          int kpiv = pstrategy == 1 ? A.find_pivot(k) : A.find_scaled_pivot(k);
          A.swap(k, kpiv);
          std::swap(b[k], b[kpiv]);
        }

        for(int i = k + 1; i < A.rows(); i++){
          T l = A[i][k] / A[k][k];
          for(int j = k + 1; j < A.cols(); j++){
            A[i][j] = A[i][j] - l * A[k][j];
          }
          A[i][k] = l;
        }
      }
    }

    /**
    * @brief Perform Cholesky decomposition of a square, symmetric matrix
    * @details
    * @param A - input matrix
    * @throws Runtime Error if matrix is not symmetric
    */
    template<typename T>
    void cholesky(matrix<T>& A){
      // Check if A is symmetric
      if(!A.is_symmetric())
        throw std::runtime_error("Matrix not symmetric in Cholesky Decomposition");

      // Perform decomposition into GG^T
      int n = A.rows();
      for(int k = 0; k < n - 1; k++){
        A[k][k] = std::sqrt(A[k][k]);
        for(int i = k + 1; i < n; i++)
          A[i][k] = A[i][k] / A[k][k];

        for(int j = k + 1; j < n; j++){
          for(int i = j; i < n; i++){
            A[i][j] = A[i][j] - A[i][k] * A[j][k];
            // If A is not s.p.d. A[i][j] = NaN
            if(std::isnan(A[i][j])) throw std::runtime_error("Matrix not positive definite in Cholesky Decomposition");
          }
        }
      }
      A[n - 1][n - 1] = std::sqrt(A[n - 1][n - 1]);

      // Reflect across diagonal
      for(int i = 0; i < A.rows(); i++)
        for(int j = i; j < A.cols(); j++)
          A[i][j] = A[j][i];

    }

    /**
    * @brief Check if matrix is s.p.d. using Cholesky Decomposition
    * @details
    * @param A - input matrix
    */
    template<typename T>
    bool is_spd(matrix<T> A){
      try{
        // cholesky() throws an exception
        cholesky(A);
        return true;
      } catch(...){
        // cholesky() threw an error
        // A is not s.p.d
        return false;
      }
    }

    /**
    * @brief Decompose A into QR
    * @details
    * @param A - input matrix
    */
    template<typename T>
    matrix<T> qr_factorization(matrix<T>& A){
      int n = A.rows();

      matrix<T> Q(n,n, (T) 0);

      for(int j = 0; j < n; j++){
        // q_j = a_j
        for(int i = 0; i < n; i++){
          Q[i][j] = A[i][j];
        }

        for(int i = 0; i < j; i++){
          // r_i,j = <q_j, q_i>
          T rij = 0;
          for(int k = 0; k < n; k++){
            rij += Q[k][j] * Q[k][i];
          }

          // q_j = q_j - r_i,j * q_i
          for(int k = 0; k < n; k++){
            Q[k][j] -= rij * Q[k][i];
          }
        }

        // r_j,j = ||q_j||
        T rjj = 0;
        for(int i = 0; i < n; i++){
          rjj += Q[i][j] * Q[i][j];
        }
        rjj = std::sqrt(rjj);

        // q_j = q_j / r_j,j
        for(int i = 0; i < n; i++){
          Q[i][j] /= rjj;
        }
      }

      return Q;
    }

    /**
    * @brief Find a lower bound of the condition number of a square matrix
    * @details
    * @param A - input matrix
    */
    template<typename T>
    double kappa(matrix<T> A){
      T normA = A.one_norm();
      int count = 0;
      matrix<T> invA(A.rows(), A.cols());
      for(int i = 0; i < A.rows(); i++){
        for(int j = 0; j < A.cols(); j++){
          invA[i][j] = A[j][i] / A[i][i];
          count++;
        }
      }

      std::cout << A.rows() << "," << count + 1 << std::endl;
      return normA * invA.one_norm();
    }

    /**
    * @brief Solve a tri-diagonal system of equations
    * @details
    * @param al - lower diagonal
    * @param am - main diagonal
    * @param au - upper diagonal
    * @param b - solution vector
    */
    template<typename T>
    mathx::array<T> solve(array<T> al, array<T> am, array<T> au, array<T> b){
      // Set n
      int n = b.size() - 1;

      // Factor first row
      au[0] /= am[0];
      b[0] /= am[0];

      // Factor interior rows
      for(int i = 1; i < n; i++){
        au[i] /= am[i] - al[i] * au[i - 1];
        b[i] = (b[i] - al[i] * b[i - 1]) / (am[i] - al[i] * au[i - 1]);
      }

      // Factor last row
      b[n] = (b[n] - al[n] * b[n - 1]) / (am[n] - al[n] * au[n - 1]);

      // Solve
      for(int i = n; i >= 0; i--){
        b[i] -= au[i] * b[i + 1];
      }

      return b;
    }

    /**
    * @brief Solve the linear system Ax=b where A is s.p.d.
    * @details
    * @param A - input matrix
    * @param b - solution vector
    */
    template<typename T>
    array<T> solve(matrix<T>& A, array<T>& b){
      cholesky(A);
      array<T> y = forward_substitution(A, b);
      return back_substitution(A, y);
    }

    /**
    * @brief Solve the linear system Ax=b
    * @details
    * @param A - input matrix
    * @param b - solution vector
    * @param strategy - flag for the method used to solve linear system
    *                   0 = GE no pivoting+ BS
    *                   1 = GE partial pivoting + BS
    *                   2 = GE scaled partial pivoting + BS
    */
    template<typename T>
    array<T> solve(matrix<T> A, array<T> b, int strategy){
      switch(strategy){
        case 0:
        case 1:
        case 2:
        {
          gaussian_elimination(A, b, strategy);
          return back_substitution(A, b);
        }
        default: {
          return array<T>();
        }
      }
    }

    /**
    * @brief Solve the linear system Ax=b
    * @details
    * @param A - input matrix
    * @param b - solution vector
    * @param LU - a reference to store the LU factorization for later use
    * @param strategy - flag for the method used to solve linear system
    *                   0 = LU no pivoting + FS & BS
    *                   1 = LU partial pivoting + FS & BS
    *                   2 = LU scaled pivoting + FS & BS
    */
    template<typename T>
    array<T> solve(matrix<T> A, array<T> b, matrix<T>& LU, int strategy){
      switch(strategy){
        case 0:
        case 1:
        case 2: {
          lu(A, b, strategy);
          LU = A;
          array<T> y = forward_substitution(LU, b, true);
          return back_substitution(LU, y);
        }
        default: {
          return array<T>();
        }
      }

    }

    /**
    * @brief Solve the Least Squares via Normal Equations
    * @details
    * @param A - input matrix
    * @param b - solution vector
    */
    template<typename T>
    array<T> least_squares(matrix<T>& A, array<T>& b){
      // Compute (A^T)A
      matrix<T> B = mult_transpose(A);

      // Compute (A^T)b
      array<T> y = product(A, b, true);

      // Use cholesky factorization to solve
      return solve(B, y);
    }


    /**
    * @brief Solve the Least Squares via QR Factorization
    * @details
    * @param A - input matrix
    * @param b - solution vector
    */
    template<typename T>
    array<T> least_squares_QR(matrix<T>& A, array<T>& b){
      // Factorize A into Q
      matrix<T> Q = qr_factorization(A);

      // Compute Q transpose
      matrix<T> qT = transpose(Q);

      // Compute R from Q transpose x A
      matrix<T> R = product(qT, A);

      // Compute C from Q transpose x b
      array<T> c = product(qT, b);
      
      // Use back substitution to solve Rx = c
      return back_substitution(R,c);
    }
  }
}
