/**
 * @file tf01md.h
 * @brief C wrapper for SLICOT routine TF01MD
 *
 * This file provides a C interface to the SLICOT routine TF01MD,
 * which computes the output response sequence of a linear time-invariant
 * discrete-time system given its state-space model (A,B,C,D) and
 * an input sequence.
 */

 #ifndef TF01MD_H
 #define TF01MD_H
 
 #include <stddef.h> // For size_t
 
 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Computes the output sequence y(k) for k=1..NY for a discrete-time system.
  *
  * Calculates y(k) = C*x(k) + D*u(k) and x(k+1) = A*x(k) + B*u(k),
  * given the initial state x(1) and the input sequence u(1)...u(NY).
  *
  * @param[in] n         Order of matrix A, n >= 0.
  * @param[in] m         Number of system inputs, m >= 0.
  * @param[in] p         Number of system outputs, p >= 0.
  * @param[in] ny        Number of output vectors y(k) to compute, ny >= 0.
  * @param[in] a         Double array, dimension (lda, n) or (n, lda). State matrix A.
  * @param[in] lda       Leading dimension of A. >= max(1,n).
  * @param[in] b         Double array, dimension (ldb, m) or (n, ldb). Input matrix B.
  * @param[in] ldb       Leading dimension of B. >= max(1,n).
  * @param[in] c         Double array, dimension (ldc, n) or (p, ldc). Output matrix C.
  * @param[in] ldc       Leading dimension of C. >= max(1,p).
  * @param[in] d         Double array, dimension (ldd, m) or (p, ldd). Direct link matrix D.
  * @param[in] ldd       Leading dimension of D. >= max(1,p).
  * @param[in] u         Double array, dimension (ldu, ny) or (m, ldu). Input sequence U = [u(1) u(2) ... u(NY)].
  * @param[in] ldu       Leading dimension of U. >= max(1,m).
  * @param[in,out] x     Double array, dimension (n).
  * On entry, the initial state vector x(1).
  * On exit, the final state vector x(NY+1).
  * @param[out] y        Double array, dimension (ldy, ny) or (p, ldy). Output sequence Y = [y(1) y(2) ... y(NY)].
  * @param[in] ldy       Leading dimension of Y. >= max(1,p).
  * @param[in] row_major Integer flag: 0 for column-major, 1 for row-major.
  * Affects interpretation of lda, ldb, ldc, ldd, ldu, ldy.
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value
  * Memory allocation errors may also be returned.
  */
 int slicot_tf01md(int n, int m, int p, int ny,
                   const double* a, int lda, const double* b, int ldb,
                   const double* c, int ldc, const double* d, int ldd,
                   const double* u, int ldu,
                   double* x, double* y, int ldy, int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* TF01MD_H */
 