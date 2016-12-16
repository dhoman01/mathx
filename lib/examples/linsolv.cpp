#include "linsolv.hpp"
#include <cstdlib>
#include <iostream>
#include <tuple>

mathx::matrix<double> generateSPD(int n);

int main(){
  /* Define A and x */
  mathx::matrix<double> A = generateSPD(5);
  mathx::array<double> x(5, 1);

  /* Define LU */
  mathx::matrix<double> LU;

  /* Use matmul() to create b*/
  mathx::array<double> b = mathx::linsolv::matmul(A, x);

  double kappa = mathx::linsolv::kappa(A);

  /* Solve system using Gaussian elimination
  *  without pivoting
  */
  mathx::array<double> guass_elim = mathx::linsolv::solve(A, b, 0);

  /* Solve system using LU factorization
  *  storing the factorization into LU
  *  with no pivoting
  */
  mathx::array<double> lu_fact = mathx::linsolv::solve(A, b, LU, 0);

  /* Solve system using Cholesky factorization */
  mathx::matrix<double> Acopy = A;
  mathx::array<double> cholesky = mathx::linsolv::solve(Acopy, b);

  /* Initialize iteration variables */
  int maxiter = 10000;
  double tol = std::pow(10, -16);
  mathx::array<double> x0(5, -1);

  /* Solve system using Jacobi iteration */
  mathx::array<double> jacobi = mathx::linsolv::jacobi(A, b, x0, tol, maxiter);

  /* Solve system using Gauss-Seidel */
  mathx::array<double> gauss_seidel = mathx::linsolv::gauss_seidel(A, b, x0, tol, maxiter);

  /* Solve system using Conjugate Gradient Method */
  mathx::array<double> cgm = mathx::linsolv::cgm(A, b, x0, tol, maxiter);

  /* Find max eigenvalue and corresponding eigenvector */
  std::pair<double, mathx::array<double>> max_eigen = mathx::linsolv::power_method(A, x0, tol, maxiter);

  /* Find min eigenvalue and corresponding eigenvector */
  std::pair<double, mathx::array<double>> min_eigen = mathx::linsolv::inverse_power_method(A, x0, .008, tol, maxiter);

  /* Compute inverse of A */
  mathx::matrix<double> Ainv = mathx::linsolv::inverse(A);

  /* Compute QR factorization of A */
  mathx::matrix<double> Q = mathx::linsolv::qr_factorization_mgs(A);

  /* Compute LU factorization of A */
  LU = mathx::linsolv::lu(A, x0);

  /* Print results of linear systems */
  std::cout << "Gaussian elimination error   = " << mathx::vectors::norm(guass_elim - x) << std::endl;
  std::cout << "LU factorization error       = " << mathx::vectors::norm(lu_fact - x) << std::endl;
  std::cout << "Cholesky factorization error = " << mathx::vectors::norm(cholesky - x) << std::endl;
  std::cout << "Jacobi iteration error       = " << mathx::vectors::norm(jacobi - x) << std::endl;
  std::cout << "Gauss-Seidel error           = " << mathx::vectors::norm(gauss_seidel - x) << std::endl;
  std::cout << "CGM error                    = " << mathx::vectors::norm(cgm - x) << std::endl;

  /* Print other results */
  std::cout << "k(A) = " << kappa << std::endl;
  std::cout << "max eigen " << max_eigen.first << " vector is " << max_eigen.second.to_string() << std::endl;
  std::cout << "min eigen " << min_eigen.first << " vector is " << max_eigen.second.to_string() << std::endl;

  /***** OUTPUTS ****/
  // Gaussian elimination error   = 7.35098e-15
  // LU factorization error       = 7.35098e-15
  // Cholesky factorization error = 1.14676e-14
  // Jacobi iteration error       = 4.78041e-15
  // Gauss-Seidel error           = 2.89724e-15
  // CGM error                    = 7.13313e-15
  // k(A) = 271210
  // max eigen 1714.25 vector is [ -0.0191586 -0.0127752 -0.0378184 -0.0760822 -0.996118  ]^T
  // min eigen 1.62375 vector is [ -0.0191586 -0.0127752 -0.0378184 -0.0760822 -0.996118  ]^T

  return EXIT_SUCCESS;
}

/* Method to generate a matrix
 * that is s.p.d.
 */
mathx::matrix<double> generateSPD(int n){
  mathx::matrix<double> A = mathx::matrix<double>(n, n, true);

  return mathx::linsolv::mult_transpose(A);
}
