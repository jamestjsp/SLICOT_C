/**
 * @file sb10jd.h
 * @brief C wrapper for SLICOT routine SB10JD
 *
 * This file provides a C interface to the SLICOT routine SB10JD,
 * which converts a descriptor state-space system into regular
 * state-space form.
 */

 #ifndef SB10JD_H
 #define SB10JD_H
 
 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Converts a descriptor system E*dx/dt = A*x + B*u, y = C*x + D*u to state-space form.
  *
  * Converts the descriptor system (E, A, B, C, D) into an equivalent
  * regular state-space system (Ad, Bd, Cd, Dd) where Ad = inv(E)*A, etc.
  * The input matrices A, B, C, D are overwritten with the resulting Ad, Bd, Cd, Dd.
  * The input matrix E is destroyed.
  *
  * @param[in] n         The order of the descriptor system (order of A, E), n >= 0.
  * @param[in] m         The number of inputs (columns of B, D), m >= 0.
  * @param[in] np        The number of outputs (rows of C, D), np >= 0.
  * @param[in,out] a     Double array, dimension (lda, n) or (n, lda).
  * On entry, the matrix A. On exit, the state matrix Ad of the converted system.
  * @param[in] lda       Leading dimension of A. >= max(1,n).
  * @param[in,out] b     Double array, dimension (ldb, m) or (n, ldb).
  * On entry, the matrix B. On exit, the input matrix Bd of the converted system.
  * @param[in] ldb       Leading dimension of B. >= max(1,n).
  * @param[in,out] c     Double array, dimension (ldc, n) or (np, ldc).
  * On entry, the matrix C. On exit, the output matrix Cd of the converted system.
  * @param[in] ldc       Leading dimension of C. >= max(1,np).
  * @param[in,out] d     Double array, dimension (ldd, m) or (np, ldd).
  * On entry, the matrix D. On exit, the feedthrough matrix Dd of the converted system.
  * @param[in] ldd       Leading dimension of D. >= max(1,np).
  * @param[in,out] e     Double array, dimension (lde, n) or (n, lde).
  * On entry, the descriptor matrix E. On exit, E is destroyed.
  * @param[in] lde       Leading dimension of E. >= max(1,n).
  * @param[out] nsys     The order of the converted state-space system (should equal N if E is nonsingular).
  * @param[in] row_major Integer flag:
  * = 0: Arrays a, b, c, d, e are column-major (Fortran style).
  * = 1: Arrays a, b, c, d, e are row-major (C style).
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * = 1: SVD algorithm failed to converge.
  * Memory allocation errors may also be returned.
  */
 SLICOT_EXPORT
 int slicot_sb10jd(int n, int m, int np,
                   double* a, int lda, double* b, int ldb,
                   double* c, int ldc, double* d, int ldd,
                   double* e, int lde, int* nsys,
                   int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* SB10JD_H */
 