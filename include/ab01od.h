/**
 * @file ab01od.h
 * @brief C wrapper for SLICOT routine AB01OD
 *
 * This file provides a C interface to the SLICOT routine AB01OD,
 * which reduces a state-space system (A,B) to staircase form using
 * orthogonal state and input transformations.
 */

 #ifndef AB01OD_H
 #define AB01OD_H

 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Reduces the system (A,B) to staircase form using orthogonal transformations.
  *
  * Computes orthogonal matrices U and V such that Ac = U'*A*U and Bc = U'*B*V
  * are in upper staircase form. Determines the controllable subspace order NCONT,
  * the controllability index INDCON, and the dimensions of the stairs KSTAIR.
  * The reduction can be performed in stages (forward, backward, or all).
  *
  * @param[in] stages    Specifies the reduction stages to perform:
  * = 'F': Forward stage only (orthogonal canonical form).
  * = 'B': Backward stage only (triangularizes subdiagonals of existing form).
  * = 'A': All stages (forward then backward).
  * @param[in] jobu      Specifies whether to compute the state transformation U:
  * = 'N': Do not compute U.
  * = 'I': Compute U explicitly. If stages='B', U is updated; otherwise, it's computed from scratch.
  * @param[in] jobv      Specifies whether to compute the input transformation V:
  * = 'N': Do not compute V (not referenced if stages='F').
  * = 'I': Compute V explicitly (not referenced if stages='F').
  * @param[in] n         The order of the matrix A, n >= 0.
  * @param[in] m         The number of system inputs (columns of B), m >= 0.
  * @param[in,out] a     Double array, dimension (lda, n) (col-major) or (n, lda) (row-major).
  * On entry, the state matrix A. If stages='B', A must be in orthogonal canonical form.
  * On exit, contains the transformed matrix Ac = U'*A*U in staircase form.
  * @param[in] lda       The leading dimension of array A.
  * If row_major=0, lda >= max(1,n).
  * If row_major=1, lda >= max(1,n) (number of columns).
  * @param[in,out] b     Double array, dimension (ldb, m) (col-major) or (n, ldb) (row-major).
  * On entry, the input matrix B. If stages='B', B must be in orthogonal canonical form.
  * On exit, contains the transformed matrix Bc = U'*B*V (if stages!='F') or U'*B (if stages='F').
  * @param[in] ldb       The leading dimension of array B.
  * If row_major=0, ldb >= max(1,n).
  * If row_major=1, ldb >= max(1,m) (number of columns).
  * @param[in,out] u     Double array, dimension (ldu, n) (col-major) or (n, ldu) (row-major).
  * If jobu='I' and stages='B', U on entry must be the transformation from the forward stage.
  * On exit, if jobu='I', contains the orthogonal state transformation matrix U.
  * If jobu='N', it is not referenced (can be NULL, ldu ignored).
  * @param[in] ldu       The leading dimension of array U.
  * If jobu='I':
  * If row_major=0, ldu >= max(1,n).
  * If row_major=1, ldu >= max(1,n) (number of columns).
  * If jobu='N', ldu >= 1 (but ignored).
  * @param[out] v        Double array, dimension (ldv, m) (col-major) or (m, ldv) (row-major).
  * On exit, if jobv='I' and stages!='F', contains the M-by-M orthogonal input transformation matrix V.
  * If jobv='N' or stages='F', it is not referenced (can be NULL, ldv ignored).
  * @param[in] ldv       The leading dimension of array V.
  * If jobv='I' and stages!='F':
  * If row_major=0, ldv >= max(1,m).
  * If row_major=1, ldv >= max(1,m) (number of columns).
  * If jobv='N' or stages='F', ldv >= 1 (but ignored).
  * @param[in,out] ncont On entry (if stages='B'), the order of the controllable part.
  * On exit (if stages!='B'), the computed order of the controllable part.
  * @param[in,out] indcon On entry (if stages='B'), the controllability index.
  * On exit (if stages!='B'), the computed controllability index.
  * @param[in,out] kstair Integer array, dimension (n).
  * On entry (if stages='B'), the dimensions of the stairs.
  * On exit (if stages!='B'), the computed dimensions of the stairs (first indcon elements).
  * @param[in] tol       Tolerance for rank determination (not referenced if stages='B').
  * If tol > 0, used as lower bound for reciprocal condition number.
  * If tol <= 0, a default tolerance is computed internally.
  * @param[in] row_major Integer flag:
  * = 0: Arrays a, b, u, v are column-major (Fortran style).
  * = 1: Arrays a, b, u, v are row-major (C style).
  * (Array kstair is 1D and unaffected).
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * Memory allocation errors may also be returned.
  */
 SLICOT_C_WRAPPER_API
 int slicot_ab01od(char stages, char jobu, char jobv,
                   int n, int m,
                   double* a, int lda,
                   double* b, int ldb,
                   double* u, int ldu,
                   double* v, int ldv,
                   int* ncont, int* indcon, int* kstair,
                   double tol, int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* AB01OD_H */
 