#include <cfloat>  // for std::isnan()
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <string>

namespace mathx {
  /*! This namespace is the workhorse of the package. The focus of this namespace is solving systems of linear equations. However, several utility methods where included in this namespace to simplify internal access.\n\n
   *  The code is ordered in the following sections
   *     1. Utility Methods - Methods to perform matrix operations such as multiplication
   *     2. Factorizations - Methods to factor matrices
   *     3. Iterative Methods - Iterative methods for solving systems of linear equations
   *     4. Matrix Util Methods - Methods used to get information such as the inverse of a matrix
   *     5. Solve Wrappers - These are helper methods that wrap the meta-algorithms to solve linear systems in a single function call
   *     6. Least Squares Wrappers - These are helper methods that wrap the meta-algorithms to solve least squares problems into a single function call
   *  \n\n Now that we have covered what is in this namespace, let us have some discussion on the merits of the different approaches to solving linear systems. Beginning with direct methods such as Gaussian elimination and LU fractorization.\n\n
   *  Gaussian elimination is robust especially when combined with pivoting, however it has an asymptotic complexity of \f$O(n^3)\f$. If you are solving a generic system of equations with little restraint on the matrix \f$A\f$ and you only need to solve once then Gaussian elimination is a solid algorithm. However, if you are solving multiple equations with the same matrix, then LU factorization will allow you to speed up subsequent solutions after the intital decomposition. LU factorization also has an asymptotic complexity of \f$O(n^3)\f$.\n\n
   *  Direct methods work well on matrices that are small (\f$n<10000\f$), however as \f$n\f$ grows larger the more round-off error accumulates, thus rendering the results useless. To overcome this issue of round-off error accumulation, you should use an iterative approach. This namespace contains three, Jacobi iteration, Gauss-Seidel iteration, and the Conjugate Gradient Method (CGM). All three methods require strictly positive definite matrices. CGM converges much quicker than either Jacobi or Gauss-Seidel.
   */
  namespace linsolv {

    /********************************************/
    /****           UTILITY METHODS          ****/
    /********************************************/

    /**
    * @brief Multiply a matrix by a vector
    * @param A - input matrix
    * @param x - input vector
    * @returns b - an array<T> that is the product of the action of A on x
    */
    template<typename T>
    array<T> matmul(matrix<T>& A, array<T>& x, bool a_trans = false){
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
    * @param A - input matrix
    * @param x - input vector
    * @returns b - an array<T> that is the product of the action of A on x
    */
    template<typename T>
    array<T> mamtul(array<double> al, array<double> am, array<double>au, array<T>& x){
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
    * @param A - input matrix
    * @param B - input matrix
    * @returns C - a matrix<T> that is the product of AB
    */
    template<typename T>
    matrix<T> matmul(matrix<T>& A, matrix<T>& B){
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
    * @returns A^T - a matrix<T> that is the transpose of the input matrix A
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
    * @details Thus method uses the fact that \f$A^TA\f$ is simply the product of the symmetric entries.
    * @param A - input matrix
    * @returns B - a matrix<T> that is the product of A and its transpose
    */
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
    * @brief Perform backwards substitution to solve Ux=b
    * @details Backwards substitution uses an upper traingular matrix to solve \f$U\textbf{x}=\textbf{b}\f$, where \f[x_k=\frac{b_k-\sum_{j=k+1}^na_{kj}x_j}{a_{kk}}\quad@cite AscherGrief\f]
    * @param - U an upper triangular matrix
    * @param - b a vector of values for the right-hand side of the equation
    * @returns x - an array<T> that is the solution of Ux=b
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
    * @details Forward substitution uses a lower triangular matrix to solve \f$L\textbf{x}=\textbf{b}\f$, where \f[x_k=\frac{b_k-\sum_{j=1}^{k-1}a_{kj}x_j}{a_{kk}}\quad@cite AscherGrief\f]
    * @param - L a lower triangular matrix
    * @param - b a vector of values for the right-hand side of the equation
    * @param - isLU a flag to interpret D as all ones
    * @returns x - an array<T> that is the solution fo Lx=b
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


    /********************************************/
    /****           FACTORIZATIONS           ****/
    /********************************************/

    /**
    * @brief Factor a square matrix A into L and U
    * @details The process of Gaussian elimination factors a matrix into \f$L\f$ and \f$U\f$. This method returns that decomposition.
    * @param A - input matrix
    * @param b - solution vector (used in pivoting)
    * @param pstrategy - flag declaring the pivoting strategy
    *                    0 = no pivoting
    *                    1 = partial pivoting
    *                    2 = scaled partial pivoting
    * @returns LU - a matrix<T> that is the LU decompostion of A
    */
    template<typename T>
    matrix<T> lu(matrix<T>& A, array<T>& b, int pstrategy = 0){
      matrix<T> LU(A.rows(), A.cols());
      LU = A;
      int m = A.rows();
      int n = A.cols();
      for(int k = 0; k < m - 1; k++){
        if(pstrategy > 0){
          int kpiv = pstrategy == 1 ? LU.find_pivot(k) : LU.find_scaled_pivot(k);
          LU.swap_row(k, kpiv);
          std::swap(b[k], b[kpiv]);
        }

        for(int i = k + 1; i < m; i++){
          T l = LU[i][k] / LU[k][k];
          for(int j = k + 1; j < n; j++){
            LU[i][j] = LU[i][j] - l * LU[k][j];
          }
          LU[i][k] = l;
        }
      }

      return LU;
    }

    /**
    * @brief Perform Cholesky decomposition of a s.p.d matrix
    * @details Cholesky decomposition is defined as \f$A=GG^{T}\f$ where \f$G=LD^{1/2}\f$ @cite AscherGrief This method is destructive to A
    * @param A - input matrix
    * @throws Runtime Error if matrix is not symmetric
    * @returns - nothing as this method modifies A in place
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
    * @details A matrix \f$A\f$ is s.p.d. if \f$A\in R^{nxn}\f$ and \f$A_{i,j}=A_{j,i}\f$ and all eigenvalues of \f$A\f$ are positive. Computing eigenvalues is complex, however there is a simple test. If the matrix \f$A\f$ has a Cholesky factorization it is s.p.d.
    * @param A - input matrix
    * @returns is_spd - a boolean indicating if the input matrix is s.p.d.
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
    * @brief Decompose A into QR (where R = (Q^T)A) using MGS
    * @details Though the classical Gram-Shmidt algorithm is elegant it is also numerically unstable for columns of \f$A\f$ that are nearly linearly dependant @cite AscherGrief A simple fix is to use the already computed columns of \f$Q\f$ to find the jth column.
    * @param A - input matrix
    * @returns QR - a matrix<T> that is the QR factorization of the input matrix
    */
    template<typename T>
    matrix<T> qr_factorization_mgs(matrix<T>& A){
      int n = A.rows();

      matrix<T> Q(n,n, (T) 0);

      for(int j = 0; j < n; j++){
        // Set the jth column of
        // Q to the jth column of A
        for(int i = 0; i < n; i++){
          Q[i][j] = A[i][j];
        }

        for(int i = 0; i < j; i++){
          // r_i,j = equals the dot
          // product of the jth and
          // ith columns of Q
          T rij = 0;
          for(int k = 0; k < n; k++){
            rij += Q[k][j] * Q[k][i];
          }

          // q_j = q_j - r_i,j * q_i
          for(int k = 0; k < n; k++){
            Q[k][j] -= rij * Q[k][i];
          }
        }

        // Normalize the jth column of Q
        T rjj = 0;
        for(int i = 0; i < n; i++){
          rjj += Q[i][j] * Q[i][j];
        }
        rjj = std::sqrt(rjj);

        for(int i = 0; i < n; i++){
          Q[i][j] /= rjj;
        }
      }

      return Q;
    }


    /********************************************/
    /****         ITERATIVE METHODS          ****/
    /********************************************/

    /**
    * @brief Find solution to linear system using Jacobi Iteration
    * @details Jacobi iteration defines \f[\textbf{x}_{k+1} = \textbf{x}_k + D^{-1}\textbf{r}_k\quad@cite AscherGrief\f] Jacobi iteration belongs to relaxation methods. As such Jacobi iteration will only converge for strictly diagonally dominant matrices.
    * @param A - a strictly diagonally dominant matrix
    * @param b - solution vector
    * @param x0 - initial guess
    * @param tol - error tolerance
    * @param maxiter - maximum number of iterations to perform
    * @param debug (false) - flag to print debug info
    * @returns x - an array<T> that is the solution to Ax=b
    */
    template<typename T>
    array<T> jacobi(matrix<T>& A, array<T>& b, array<T>& x0, double tol, int maxiter, bool debug = false){
      // initialize variables;
      array<T> xkp1;
      array<T> xk = x0;
      int iter = 0;
      int n = A.cols();
      double error = tol * 10;

      // Perform iterations until stopping
      // criteria are met
      while(iter < maxiter && error > tol){
        // Compute x^(k+1)[i] for i in [0,n);
        xkp1 = b;
        for(int i = 0; i < n; i++){
          for(int j = 0; j < i; j++)
            xkp1[i] -= A[i][j] * xk[j];

          for(int j = i + 1; j < n; j++)
            xkp1[i] -= A[i][j] * xk[j];

          xkp1[i] /= A[i][i];
        }

        // Calculate error
        error = vectors::norm(xkp1 - xk);

        // Assign new variables
        xk = xkp1;
        iter++;
      }

      if(debug) std::cout << n << ", " << iter << std::endl;

      return xkp1;
    }

    /**
    * @brief Find solution of linear system using Gauss-Seidel
    * @details Gauss-Seidel iteration defines \f[\textbf{x}_{k+1}=\textbf{x}_k+E^{-1}\textbf{r}_k\quad@cite AscherGrief\f] Gauss-Seidel iteration belongs to relaxation methods. As such Jacobi iteration will only converge for strictly diagonally dominant matrices. Gauss-Seidel converges at twice the rate of Jacobi.
    * @param A - a strictly diagonally dominant matrix
    * @param b - solution vector
    * @param x0 - initial guess
    * @param tol - error tolerance
    * @param maxiter - maximum number of iterations to perform
    * @param debug (false) - flag to print debug info
    * @returns x - an array<T> that is the solution of Ax=b
    */
    template<typename T>
    array<T> gauss_seidel(matrix<T>& A, array<T>& b, array<T>& x0, double tol, int maxiter, bool debug = false){
      // initialize variables;
      array<T> xkp1;
      array<T> xk = x0;
      int iter = 0;
      int n = A.cols();
      double error = tol * 10;

      // Perform iterations until stopping
      // criteria are met
      while(iter < maxiter && error > tol){
        // Compute x^(k+1)[i] for i in [0,n);
        xkp1 = b;
        for(int i = 0; i < n; i++){
          for(int j = 0; j < i; j++)
            xkp1[i] -= A[i][j] * xkp1[j];

          for(int j = i + 1; j < n; j++)
            xkp1[i] -= A[i][j] * xk[j];

          xkp1[i] /= A[i][i];
        }

        // Calculate error
        error = vectors::norm(xkp1 - xk);

        // Assign new variables
        xk = xkp1;
        iter++;
      }

      if(debug) std::cout << n << ", " << iter << std::endl;

      return xkp1;
    }

    /**
    * @brief Solve linear system using Conjugate Gradient method
    * @details The Conjugate Gradient method (CGM) overcomes a weakness of stationary methods in that it uses information gathered throughout its iterations. CGM defines \f[\textbf{x}_{k+1}=\textbf{x}_k+\alpha\textbf{p}_k\quad@cite AscherGrief\f] Where the vector \f$\textbf{p}_k\f$ is the search direction and the scalar \f$\alpha\f$ is the step size @cite AscherGrief
    * @param A - a strictly diagonally dominant matrix
    * @param b - solution vector
    * @param x0 - initial guess
    * @param tol - error tolerance
    * @param maxiter - maximum number of iterations to perform
    * @returns x - an array<T> that is the solution of Ax=b
    */
    template<typename T>
    array<T> cgm(matrix<T>& A, array<T>& b, array<T>& x0, double tol, int maxiter){
      // Initialize x^(k+1) and x^(k)
      array<T> xkp1;
      array<T> xk = x0;

      // Initialize r^(k+1) and r^(k)
      array<T> rkp1;
      array<T> rk = b - matmul(A,x0);

      // Initialize p^(k+1) and p^(k)
      array<T> pkp1;
      array<T> pk = rk;

      // Initialize tolerance and delta
      tol = std::pow(tol, 2);
      double deltak = vectors::dot_product(rk, rk);
      double deltakp1;

      // Initialize b delta
      double bdelta = vectors::dot_product(b, b);

      int iter = 0;
      int n = A.cols();
      while(deltak > tol * bdelta && iter < maxiter){
        array<T> sk = matmul(A, pk);
        double alphak = deltak / vectors::dot_product(pk,sk);
        xkp1 = xk;
        rkp1 = rk;
        // Iterate to find x^(k+1)[i]
        // and r^(k+1)[i]
        for(int i = 0; i < n; i++){
          xkp1[i] += alphak * pk[i];
          rkp1[i] -= alphak * sk[i];
        }

        // Find delta k+1 and p^(k+1)
        deltakp1 = vectors::dot_product(rkp1, rkp1);
        pkp1 = rkp1;
        for(int i = 0; i < n; i++)
          pkp1[i] += (deltakp1 / deltak) * pk[i];

        // Assign new values
        xk = xkp1;
        pk = pkp1;
        rk = rkp1;
        deltak = deltakp1;

        iter++;
      }

      return xkp1;
    }

    /********************************************/
    /****        MATRIX UTIL METHODS         ****/
    /********************************************/

    /**
    * @brief Use the power method to find the largest eigenvalue and corresponding eigenvector of a matrix
    * @details The power method finds the largest eigenvalue and corresponding eigenvector via an iterative approach.
    * @param A - input matrix
    * @param v0 - initial guess
    * @param tol - error tolerance
    * @param maxiter - max iterations to perform
    * @param debug - Print debug info (default=false)
    * @returns \f$\lambda,\textbf{v}\f$ - a pair<T, array<T>> that is the pair of the largest eigenvalue of \f$A\f$ and its corresponding eigenvector
    */
    template<typename T>
    std::pair<T,array<T>> power_method(matrix<T>& A, array<T>& v0, double tol, int maxiter, bool debug=false){
      // Initialize variables
      array<T> vkm1 = v0;
      array<T> vk;
      T lambda = 0;
      T lambdakm1 = 10;
      int iter = 0;
      double error = 10 * tol;

      if(debug) std::cout << "Iterations, Error, n" << std::endl;

      // Moving the matmul outside
      // the loop improves performance
      // and overcomes an issue of
      // overflow when n and iter are large
      array<T> v = matmul(A, vkm1);
      while(iter++ < maxiter && error > tol){
        // Normalize v and assing to vk;
        vk = vectors::normalize(v);

        // Calculate lambda_k
        lambda = vectors::dot_product(vk, v);

        // Calculate error
        error = std::abs(lambda - lambdakm1);

        if(debug) std::cout << iter << "," << error << "," << A.cols() << std::endl;

        // Reinitialize values for
        // the next iteration
        v = matmul(A, vk);
        lambdakm1 = lambda;
        vkm1 = v;
      }

      return std::make_pair(lambda, vk);
    }

    /**
    * @brief Shift a matrix by alpha
    * @param A - matrix to shift
    * @param alpha - alpha
    * @returns A - a matrix<T> shifted by \f$\alpha\f$
    */
    template<typename T>
    matrix<T> shift(matrix<T>& A, T alpha){
      matrix<T> shifted = A;
      for(int i = 0; i < shifted.rows(); i++)
        shifted[i][i] -= alpha;

      return shifted;
    }

    /**
    * @brief Use the inverse power method to find smallest eigenvalue and corresponding eigenvector
    * @details The power method finds the smallest eigenvalue and corresponding eigenvector via an iterative approach.
    * @param A - input matrix
    * @param v0 - intital guess
    * @param alpha - shift value
    * @param tol - error tolerance
    * @param maxiter - max iterations to perform
    * @param debug - Print debug info (default=false)
    * @returns \f$\lambda,\textbf{v}\f$ - a pair<T, array<T>> that is the pair of the smallest eigenvalue of \f$A\f$ and its corresponding eigenvector
    */
    template<typename T>
    std::pair<T, array<T>> inverse_power_method(matrix<T>& A, array<T>& v0, double alpha, double tol, int maxiter, bool debug=false){
      // Initialize variables
      array<T> vkm1 = v0;
      array<T> vk;
      T lambda = 0;
      T lambdakm1 = 10;
      double error = 10 * tol;
      int iter = 0;

      // Shift A by alpha
      A = shift(A, alpha);

      // Factor A* into L and U
      matrix<T> LU = lu(A, vkm1, 0);

      if(debug) std::cout << "Iterations, Error, n" << std::endl;

      while(iter++ < maxiter && error > tol){
        // Solve Shifted * v = vkm1
        array<T> b = forward_substitution(LU, vkm1, true);
        array<T> v = back_substitution(LU, b);

        // Normalize v and assing to vk;
        vk = vectors::normalize(v);

        // Calculate lambda_k
        lambda = vectors::dot_product(vk, matmul(A, vk));

        // Calculate error
        error = std::abs(lambda - lambdakm1);

        if(debug) std::cout << iter << "," << error << "," << A.cols() << std::endl;

        // Reinitialize values for
        // the next iteration
        lambdakm1 = lambda;
        vkm1 = v;
      }

      return std::make_pair(lambda, vk);
    }

    /**
    * @brief Computes the inverse of a matrix
    * @details This method computes the inverse by factoring \f$A\f$ into \f$L\f$ and \f$U\f$. Where \f$L\f$ and \f$U\f$ are strictly lower and upper triangular, respectively. The resulting decomposition is used to solve each one-spot vector where the \f$k^{th}\f$ solution is the \f$k^{th}\f$ column of \f$A^{-1}\f$.
    * @param A - matrix to compute inverse
    * @returns \f$A^{-1}\f$ - a matrix<T> that is the inverse of the input matrix
    */
    template<typename T>
    matrix<T> inverse(matrix<T>& A){
      array<T> onespot(A.cols(), 0);

      // Decompose A into L & U
      // Passing one-spot as placeholder
      // as lu() is expecting a b for
      // pivoting
      matrix<T> LU = lu(A, onespot);

      // Calculate inverse of A
      matrix<T> Ainv(A.rows(), A.cols());

      // Use LU to solve for each one-spot
      for(int k = 0; k < A.cols(); k++){
        onespot = array<T>(A.cols(), 0);
        onespot[k] = 1;
        array<T> y = forward_substitution(A, onespot, true);
        array<T> x = back_substitution(A,y);

        // kth solution is kth column
        // of the inverse
        for(int i = 0; i < A.cols(); i++)
          Ainv[i][k] = x[i];
      }

      return Ainv;
    }

    /**
    * @brief Find a lower bound of the condition number of a square matrix
    * @details The condition number of a matrix is defined as \f$||A||\cdot||A^{-1}||\f$.
    * @param A - input matrix
    * @param norm_type - 0 -> one norm; 1 -> infinity norm
    * @returns k - the approximation of \f$k(A)\f$
    */
    template<typename T>
    double kappa(matrix<T>& A, int norm_type=0){
      // Find the inverse of A
      matrix<T> Ainv = inverse(A);
      // Return the condition number
      // k(A) = ||A||*||A^-1||
      switch(norm_type){
        case 0:
        return A.one_norm() * Ainv.one_norm();
        case 1:
        default:
        return A.infinity_norm() * Ainv.infinity_norm();
      }
    }

    /********************************************/
    /****         DIRECT METHODS          ****/
    /********************************************/

    /**
    * @brief Perform Guassian elimination on a matrix given a solution vector
    * @details Gaussian elimination uses elementry row opperations to reduce a given matrix into upper triangular form @cite AscherGrief This method is destructive to A and b.
    * @param A - the input matrix
    * @param b - the solution vector
    * @param pstrategy - flag declaring the pivoting strategy
    *                    0 = no pivoting
    *                    1 = partial pivoting
    *                    2 = scaled partial pivoting
    * @returns - nothing as A and b are changed in place
    */
    template<typename T>
    void gaussian_elimination(matrix<T>& A, array<T>& b, int pstrategy = 0){
      for(int k = 0; k < A.rows() - 1; k++){
        if(pstrategy > 0){
          int kpiv = pstrategy == 1 ? A.find_pivot(k) : A.find_scaled_pivot(k);
          A.swap_row(k, kpiv);
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

    /********************************************/
    /****           SOLVE WRAPPERS           ****/
    /********************************************/

    /**
    * @brief Solve a tri-diagonal system of equations
    * @details
    * @param al - lower diagonal
    * @param am - main diagonal
    * @param au - upper diagonal
    * @param b - solution vector
    * @returns x - an array<T> that is the solution to Ax=b where A is tri-diagonal
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
    * @returns x - an array<T> that is the solution to Ax=b where A is s.p.d.
    */
    template<typename T>
    array<T> solve(matrix<T>& A, array<T>& b){
      cholesky(A);
      array<T> y = forward_substitution(A, b);
      return back_substitution(A, y);
    }

    /**
    * @brief Solve the linear system Ax=b using Gaussian Elimination
    * @param A - input matrix
    * @param b - solution vector
    * @param strategy - flag for the method used to solve linear system
    *                   0 = GE no pivoting+ BS
    *                   1 = GE partial pivoting + BS
    *                   2 = GE scaled partial pivoting + BS
    * @returns x - an array<T> that is the solution to Ax=b
    */
    template<typename T>
    array<T> solve(matrix<T> A, array<T> b, int strategy){
      gaussian_elimination(A, b, strategy);
      return back_substitution(A, b);
    }

    /**
    * @brief Solve the linear system Ax=b
    * @param A - input matrix
    * @param b - solution vector
    * @param LU - a reference to store the LU factorization for later use
    * @param strategy - flag for the method used to solve linear system
    *                   0 = LU no pivoting + FS & BS
    *                   1 = LU partial pivoting + FS & BS
    *                   2 = LU scaled pivoting + FS & BS
    * @returns x - an array<T> that is the solution to Ax=b
    */
    template<typename T>
    array<T> solve(matrix<T>& A, array<T> b, matrix<T>& LU, int strategy){
      LU = lu(A, b, strategy);
      array<T> y = forward_substitution(LU, b, true);
      return back_substitution(LU, y);
    }

    /********************************************/
    /****       LEAST SQUARES METHODS        ****/
    /********************************************/

    /**
    * @brief Solve the Least Squares via Normal Equations
    * @param A - input matrix
    * @param b - solution vector
    * @returns x - an array<T> that is the solution to the least squares problem
    */
    template<typename T>
    array<T> least_squares(matrix<T>& A, array<T>& b){
      // Compute (A^T)A
      matrix<T> B = mult_transpose(A);

      // Compute (A^T)b
      array<T> y = matmul(A, b, true);

      // Use cholesky factorization to solve
      return solve(B, y);
    }

    /**
    * @brief Solve the Least Squares via QR Factorization
    * @param A - input matrix
    * @param b - solution vector
    * @returns x - an array<T> that is the solution to Ax=b
    */
    template<typename T>
    array<T> least_squares_QR(matrix<T>& A, array<T>& b){
      // Factorize A into Q
      matrix<T> Q = qr_factorization_mgs(A);

      // Compute Q transpose
      matrix<T> qT = transpose(Q);

      // Compute R from Q transpose x A
      matrix<T> R = matmul(qT, A);

      // Compute C from Q transpose x b
      array<T> c = matmul(qT, b);

      // Use back substitution to solve Rx = c
      return back_substitution(R,c);
    }
  }

  /** @example linsolv.cpp
  * This example demonstrates how to use the functions in the namespace ::linsolv.
  */
}
