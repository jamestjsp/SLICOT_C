/**
 * @file ab04md.h
 * @brief C wrapper for SLICOT routine AB04MD
 *
 * This file provides a C interface to the SLICOT routine AB04MD,
 * which performs a bilinear transformation on the parameters (A,B,C,D)
 * of a system, equivalent to transforming between discrete-time and
 * continuous-time representations.
 */

 #ifndef AB04MD_H
 #define AB04MD_H
 
 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Performs discrete-time <-> continuous-time conversion via bilinear transformation.
  *
  * Transforms the system matrices (A,B,C,D) between discrete-time and
  * continuous-time representations using a specified bilinear transformation
  * defined by ALPHA and BETA.
  *
  * @param[in] type      Specifies the transformation type:
  * = 'D': Discrete-time   -> Continuous-time.
  * = 'C': Continuous-time -> Discrete-time.
  * @param[in] n         The order of the matrix A, n >= 0.
  * @param[in] m         The number of system inputs (columns of B, D), m >= 0.
  * @param[in] p         The number of system outputs (rows of C, D), p >= 0.
  * @param[in] alpha     Transformation parameter ALPHA (must be non-zero).
  * @param[in] beta      Transformation parameter BETA (must be non-zero).
  * @param[in,out] a     Double array, dimension (lda, n) (col-major) or (n, lda) (row-major).
  * On entry, the state matrix A. On exit, the transformed state matrix.
  * @param[in] lda       The leading dimension of array A.
  * If row_major=0, lda >= max(1,n).
  * If row_major=1, lda >= max(1,n) (number of columns).
  * @param[in,out] b     Double array, dimension (ldb, m) (col-major) or (n, ldb) (row-major).
  * On entry, the input matrix B. On exit, the transformed input matrix.
  * @param[in] ldb       The leading dimension of array B.
  * If row_major=0, ldb >= max(1,n).
  * If row_major=1, ldb >= max(1,m) (number of columns).
  * @param[in,out] c     Double array, dimension (ldc, n) (col-major) or (p, ldc) (row-major).
  * On entry, the output matrix C. On exit, the transformed output matrix.
  * @param[in] ldc       The leading dimension of array C.
  * If row_major=0, ldc >= max(1,p).
  * If row_major=1, ldc >= max(1,n) (number of columns).
  * @param[in,out] d     Double array, dimension (ldd, m) (col-major) or (p, ldd) (row-major).
  * On entry, the input/output matrix D. On exit, the transformed input/output matrix.
  * @param[in] ldd       The leading dimension of array D.
  * If row_major=0, ldd >= max(1,p).
  * If row_major=1, ldd >= max(1,m) (number of columns).
  * @param[in] row_major Integer flag:
  * = 0: Arrays a, b, c, d are column-major (Fortran style).
  * = 1: Arrays a, b, c, d are row-major (C style).
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * = 1: if the matrix (ALPHA*I + A) is singular (TYPE='D').
  * = 2: if the matrix (BETA*I - A) is singular (TYPE='C').
  * Memory allocation errors may also be returned.
  */
 int slicot_ab04md(char type, int n, int m, int p,
                   double alpha, double beta,
                   double* a, int lda,
                   double* b, int ldb,
                   double* c, int ldc,
                   double* d, int ldd,
                   int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* AB04MD_H */
 