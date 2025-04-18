/**
 * @file ab07nd.h
 * @brief C wrapper for SLICOT routine AB07ND
 *
 * This file provides a C interface to the SLICOT routine AB07ND,
 * which computes the inverse (Ai,Bi,Ci,Di) of a given system (A,B,C,D).
 */

 #ifndef AB07ND_H
 #define AB07ND_H
 
 #include <stddef.h> // For size_t
 
 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Computes the inverse of a given linear system (A,B,C,D).
  *
  * Calculates the state-space representation (Ai,Bi,Ci,Di) of the inverse
  * of the system (A,B,C,D), assuming D is invertible. The original
  * matrices A, B, C, D are overwritten with the result.
  *
  * @param[in] n         The order of the state matrix A, n >= 0.
  * @param[in] m         The number of system inputs and outputs (D is m-by-m), m >= 0.
  * @param[in,out] a     Double array, dimension (lda, n) or (n, lda).
  * On entry, the state matrix A. On exit, the inverse state matrix Ai.
  * @param[in] lda       The leading dimension of array A. >= max(1,n).
  * @param[in,out] b     Double array, dimension (ldb, m) or (n, ldb).
  * On entry, the input matrix B. On exit, the inverse input matrix Bi.
  * @param[in] ldb       The leading dimension of array B. >= max(1,n).
  * @param[in,out] c     Double array, dimension (ldc, n) or (m, ldc).
  * On entry, the output matrix C. On exit, the inverse output matrix Ci.
  * @param[in] ldc       The leading dimension of array C. >= max(1,m).
  * @param[in,out] d     Double array, dimension (ldd, m) or (m, ldd).
  * On entry, the feedthrough matrix D. On exit, the inverse feedthrough matrix Di.
  * @param[in] ldd       The leading dimension of array D. >= max(1,m).
  * @param[out] rcond    The estimated reciprocal condition number of the original matrix D.
  * @param[in] row_major Integer flag:
  * = 0: Arrays a, b, c, d are column-major (Fortran style).
  * = 1: Arrays a, b, c, d are row-major (C style).
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * = i: (1 <= i <= M) D is exactly singular, U(i,i)=0. RCOND=0.
  * = M+1: D is numerically singular (RCOND < EPS). Computation completed.
  * Memory allocation errors may also be returned.
  */
 int slicot_ab07nd(int n, int m,
                   double* a, int lda, double* b, int ldb,
                   double* c, int ldc, double* d, int ldd,
                   double* rcond, int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* AB07ND_H */
 