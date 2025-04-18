/**
 * @file tc04ad.h
 * @brief C wrapper for SLICOT routine TC04AD
 *
 * This file provides a C interface to the SLICOT routine TC04AD,
 * which finds a state-space representation (A,B,C,D) equivalent to
 * a given left or right polynomial matrix representation.
 */

 #ifndef TC04AD_H
 #define TC04AD_H
 
 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Finds a state-space representation for a polynomial matrix fraction.
  *
  * Computes (A,B,C,D) such that C*inv(sI-A)*B + D has the same transfer
  * matrix as inv(P(s))*Q(s) (left) or Q(s)*inv(P(s)) (right).
  *
  * @param[in] leri      Specifies input type: 'L' (left fraction) or 'R' (right fraction).
  * @param[in] m         Number of system inputs, m >= 0.
  * @param[in] p         Number of system outputs, p >= 0.
  * @param[in] index     Integer array, dimension (max(m,p)).
  * If leri='L', index[0..p-1] contains row degrees of P(s).
  * If leri='R', index[0..m-1] contains column degrees of P(s).
  * @param[in] pcoeff    Double array, dimension (ldpco1, ldpco2, kpcoef). kpcoef=max(index)+1.
  * Coefficients of the denominator matrix P(s).
  * Dimensions depend on leri: (p,p,kpcoef) if 'L', (m,m,kpcoef) if 'R'.
  * Modified if leri='R', but restored on exit.
  * @param[in] ldpco1    First leading dimension of pcoeff. >= max(1,p) if 'L', >= max(1,m) if 'R'.
  * @param[in] ldpco2    Second leading dimension of pcoeff. >= max(1,p) if 'L', >= max(1,m) if 'R'.
  * @param[in] qcoeff    Double array, dimension (ldqco1, ldqco2, kpcoef).
  * Coefficients of the numerator matrix Q(s).
  * Dimensions depend on leri: (p,m,kpcoef) if 'L', (p,m,kpcoef) if 'R'.
  * Modified if leri='R', but restored on exit.
  * @param[in] ldqco1    First leading dimension of qcoeff. >= max(1,p) if 'L', >= max(1,m,p) if 'R'.
  * @param[in] ldqco2    Second leading dimension of qcoeff. >= max(1,m) if 'L', >= max(1,m,p) if 'R'.
  * @param[out] n        Order of the resulting state-space representation (sum of index).
  * @param[out] rcond    Estimated reciprocal condition number of the leading coefficient matrix of P(s).
  * @param[out] a        Double array, dimension (lda, n) or (n, lda). State matrix A.
  * @param[in] lda       Leading dimension of A. >= max(1,n).
  * @param[out] b        Double array, dimension (ldb, max(m,p)) or (n, ldb). Input matrix B. Workspace needed.
  * @param[in] ldb       Leading dimension of B. >= max(1,n).
  * @param[out] c        Double array, dimension (ldc, n) or (max(m,p), ldc). Output matrix C. Workspace needed.
  * @param[in] ldc       Leading dimension of C. >= max(1,m,p).
  * @param[out] d        Double array, dimension (ldd, max(m,p)) or (max(m,p), ldd). Direct transmission matrix D. Workspace needed.
  * @param[in] ldd       Leading dimension of D. >= max(1,m,p).
  * @param[in] row_major Integer flag: 0 for column-major, 1 for row-major.
  * Affects interpretation of lda, ldb, ldc, ldd, ldpco*, ldqco*.
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value
  * = 1: P(s) is not row/column proper (leading coefficient matrix is singular).
  * Memory allocation errors may also be returned.
  */
 SLICOT_C_WRAPPER_API
 int slicot_tc04ad(char leri, int m, int p, const int* index,
                   const double* pcoeff, int ldpco1, int ldpco2,
                   const double* qcoeff, int ldqco1, int ldqco2,
                   int* n, double* rcond,
                   double* a, int lda, double* b, int ldb,
                   double* c, int ldc, double* d, int ldd,
                   int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* TC04AD_H */
 