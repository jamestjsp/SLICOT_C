/**
 * @file sb04qd.h
 * @brief C wrapper for SLICOT routine SB04QD
 *
 * This file provides a C interface to the SLICOT routine SB04QD,
 * which solves the discrete-time Sylvester equation X + AXB = C
 * using the Hessenberg-Schur method.
 */

 #ifndef SB04QD_H
 #define SB04QD_H
 
 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Solves the discrete-time Sylvester equation X + AXB = C.
  *
  * Uses the Hessenberg-Schur method. A is reduced to upper Hessenberg form H,
  * B' is reduced to real Schur form S. The transformed equation Y + HYS' = F
  * is solved for Y, and then X = UYZ' is computed.
  * Input matrices A, B, C are overwritten.
  *
  * @param[in] n         The order of matrix A (rows of C, X), n >= 0.
  * @param[in] m         The order of matrix B (columns of C, X), m >= 0.
  * @param[in,out] a     Double array, dimension (lda, n) or (n, lda).
  * On entry, the matrix A.
  * On exit, contains Hessenberg form H and reflector data for U.
  * @param[in] lda       Leading dimension of A. >= max(1,n).
  * @param[in,out] b     Double array, dimension (ldb, m) or (m, ldb).
  * On entry, the matrix B.
  * On exit, contains the Schur form S of B'.
  * @param[in] ldb       Leading dimension of B. >= max(1,m).
  * @param[in,out] c     Double array, dimension (ldc, m) or (n, ldc).
  * On entry, the matrix C.
  * On exit, the solution matrix X.
  * @param[in] ldc       Leading dimension of C. >= max(1,n).
  * @param[out] z        Double array, dimension (ldz, m) or (m, ldz).
  * The orthogonal matrix Z from the Schur factorization of B'.
  * @param[in] ldz       Leading dimension of Z. >= max(1,m).
  * @param[in] row_major Integer flag:
  * = 0: Arrays a, b, c, z are column-major (Fortran style).
  * = 1: Arrays a, b, c, z are row-major (C style).
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * > 0: (1 <= info <= M) QR algorithm failed for B'.
  * > M: Singular matrix encountered solving for column (info-m) of X.
  * Memory allocation errors may also be returned.
  */
 SLICOT_C_WRAPPER_API
 int slicot_sb04qd(int n, int m,
                   double* a, int lda, double* b, int ldb,
                   double* c, int ldc, double* z, int ldz,
                   int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* SB04QD_H */
 