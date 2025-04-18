/**
 * @file ab13ed.h
 * @brief C wrapper for SLICOT routine AB13ED
 *
 * This file provides a C interface to the SLICOT routine AB13ED,
 * which estimates the distance from a real matrix A to the nearest
 * complex matrix with an eigenvalue on the imaginary axis.
 */

 #ifndef AB13ED_H
 #define AB13ED_H

 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

 
 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Estimates distance from matrix A to nearest matrix with eigenvalue on imaginary axis.
  *
  * Uses bisection to compute lower (LOW) and upper (HIGH) bounds for beta(A),
  * the 2-norm distance from a real matrix A to the nearest complex matrix
  * with an eigenvalue on the imaginary axis.
  * If A is stable, beta(A) is the complex stability radius.
  *
  * @param[in] n         The order of the matrix A, n >= 0.
  * @param[in] a         Double array for matrix A. Dimension (lda, n) or (n, lda).
  * @param[in] lda       The leading dimension of array A. >= max(1,n).
  * @param[out] low      Lower bound for beta(A).
  * @param[out] high     Upper bound for beta(A).
  * @param[in] tol       Relative tolerance for the bounds LOW, HIGH. If tol < sqrt(eps), sqrt(eps) is used. Recommended: 9.0.
  * @param[in] row_major Integer flag:
  * = 0: Array a is column-major (Fortran style).
  * = 1: Array a is row-major (C style).
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * = 1: QR algorithm failed to converge.
  * Memory allocation errors may also be returned.
  */
 SLICOT_C_WRAPPER_API
 int slicot_ab13ed(int n, const double* a, int lda,
                   double* low, double* high, double tol,
                   int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* AB13ED_H */
 