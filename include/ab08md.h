/**
 * @file ab08md.h
 * @brief C wrapper for SLICOT routine AB08MD
 *
 * This file provides a C interface to the SLICOT routine AB08MD,
 * which computes the normal rank of the transfer-function matrix of a
 * state-space model (A,B,C,D).
 */

 #ifndef AB08MD_H
 #define AB08MD_H
 
 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Computes the normal rank of the transfer-function matrix G(s) = C(sI-A)^-1 B + D.
  *
  * This function computes the normal rank of the transfer-function matrix
  * corresponding to the state-space model (A,B,C,D). Input matrices might be
  * scaled in place if balancing is requested.
  *
  * @param[in] equil     Specifies whether to balance the compound matrix:
  * - 'S': Perform balancing (scaling). Input matrices a,b,c,d might be modified.
  * - 'N': Do not perform balancing. Input matrices are read-only.
  * @param[in] n         The order of matrix A (number of state variables), n >= 0.
  * @param[in] m         The number of system inputs, m >= 0.
  * @param[in] p         The number of system outputs, p >= 0.
  * @param[in,out] a     Double array, dimension (lda, n) (col-major) or (n, lda) (row-major).
  * On entry, the state dynamics matrix A.
  * On exit, if equil='S', it may be overwritten by the scaled matrix.
  * @param[in] lda       The leading dimension of array A.
  * If row_major=0, lda >= max(1,n).
  * If row_major=1, lda >= max(1,n) (number of columns).
  * @param[in,out] b     Double array, dimension (ldb, m) (col-major) or (n, ldb) (row-major).
  * On entry, the input/state matrix B.
  * On exit, if equil='S', it may be overwritten by the scaled matrix.
  * @param[in] ldb       The leading dimension of array B.
  * If row_major=0, ldb >= max(1,n).
  * If row_major=1, ldb >= max(1,m) (number of columns).
  * @param[in,out] c     Double array, dimension (ldc, n) (col-major) or (p, ldc) (row-major).
  * On entry, the state/output matrix C.
  * On exit, if equil='S', it may be overwritten by the scaled matrix.
  * @param[in] ldc       The leading dimension of array C.
  * If row_major=0, ldc >= max(1,p).
  * If row_major=1, ldc >= max(1,n) (number of columns).
  * @param[in,out] d     Double array, dimension (ldd, m) (col-major) or (p, ldd) (row-major).
  * On entry, the direct transmission matrix D.
  * On exit, if equil='S', it may be overwritten by the scaled matrix.
  * @param[in] ldd       The leading dimension of array D.
  * If row_major=0, ldd >= max(1,p).
  * If row_major=1, ldd >= max(1,m) (number of columns).
  * @param[out] rank     The computed normal rank of the transfer-function matrix.
  * @param[in] tol       Tolerance used for rank decisions. If tol is less than
  * sqrt((n+p)*(n+m))*eps, a default value based on machine
  * precision is used internally by the Fortran routine.
  * @param[in] row_major Integer flag:
  * - 0: Arrays a, b, c, d are column-major (Fortran style).
  * - 1: Arrays a, b, c, d are row-major (C style).
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * Common SLICOT errors related to memory allocation might also be returned
  * (e.g., SLICOT_MEMORY_ERROR from slicot_utils.h).
  */
 SLICOT_C_WRAPPER_API
 int slicot_ab08md(char equil, int n, int m, int p,
                   double* a, int lda,
                   double* b, int ldb,
                   double* c, int ldc,
                   double* d, int ldd,
                   int* rank, double tol, int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* AB08MD_H */
 