/**
 * @file tb01pd.h
 * @brief C wrapper for SLICOT routine TB01PD
 *
 * This file provides a C interface to the SLICOT routine TB01PD,
 * which finds a reduced (minimal, controllable, or observable)
 * state-space representation (Ar,Br,Cr) in upper block Hessenberg form
 * for a given system (A,B,C).
 */

 #ifndef TB01PD_H
 #define TB01PD_H
 
 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Computes minimal, controllable, or observable block Hessenberg realization.
  *
  * Finds a reduced state-space representation (Ar,Br,Cr) for a given
  * system (A,B,C). Depending on JOB, the realization is minimal,
  * controllable, or observable. The resulting state matrix Ar is in
  * upper block Hessenberg form. Input matrices A, B, C are overwritten.
  *
  * @param[in] job       Specifies the type of reduction:
  * = 'M': Minimal realization (controllable and observable).
  * = 'C': Controllable realization.
  * = 'O': Observable realization.
  * @param[in] equil     Specifies whether to balance the triplet (A,B,C) first:
  * = 'S': Perform balancing (scaling).
  * = 'N': Do not perform balancing.
  * @param[in] n         Order of the original system A, n >= 0.
  * @param[in] m         Number of inputs (columns of B), m >= 0.
  * @param[in] p         Number of outputs (rows of C), p >= 0.
  * @param[in,out] a     Double array, dimension (lda, n) or (n, lda).
  * On entry, the matrix A. On exit, the leading nr-by-nr part contains Ar.
  * @param[in] lda       Leading dimension of A. >= max(1,n).
  * @param[in,out] b     Double array. Dimension (ldb, m) or (n, ldb) if job='C'.
  * Dimension (ldb, max(m,p)) or (n, ldb) if job='M' or 'O'.
  * On entry, the leading n-by-m part contains B. Remainder used as workspace if needed.
  * On exit, the leading nr-by-m part contains Br.
  * @param[in] ldb       Leading dimension of B. >= max(1,n).
  * @param[in,out] c     Double array, dimension (ldc, n) or (max(m,p), ldc).
  * On entry, the leading p-by-n part contains C. Remainder used as workspace if needed.
  * On exit, the leading p-by-nr part contains Cr.
  * @param[in] ldc       Leading dimension of C. >=max(1, max(m,p)) if n>0, else >=1.
  * @param[out] nr       The order of the reduced state-space representation (Ar,Br,Cr).
  * @param[in] tol       Tolerance for rank determination. If tol<=0, default used.
  * @param[out] iwork    Integer array, dimension (n + max(m,p)).
  * On exit, the first nonzero elements contain the orders of diagonal blocks of Ar.
  * @param[in] row_major Integer flag:
  * = 0: Arrays a, b, c are column-major (Fortran style).
  * = 1: Arrays a, b, c are row-major (C style).
  * (Array iwork is 1D).
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * Memory allocation errors may also be returned.
  */
 SLICOT_C_WRAPPER_API
 int slicot_tb01pd(char job, char equil, int n, int m, int p,
                   double* a, int lda, double* b, int ldb,
                   double* c, int ldc, int* nr, double tol,
                   int* iwork, int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* TB01PD_H */
 