/**
 * @file tb03ad.h
 * @brief C wrapper for SLICOT routine TB03AD
 *
 * This file provides a C interface to the SLICOT routine TB03AD,
 * which finds a relatively prime left or right polynomial matrix
 * representation equivalent to a given state-space representation.
 */

 #ifndef TB03AD_H
 #define TB03AD_H
 
 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Finds a left/right polynomial matrix representation for a state-space system.
  *
  * Computes inv(P(s))*Q(s) (left) or Q(s)*inv(P(s)) (right) such that
  * it has the same transfer matrix as T(s) = C*inv(s*I-A)*B + D.
  *
  * @param[in] leri      Specifies representation type: 'L' (left) or 'R' (right).
  * @param[in] equil     Specifies balancing: 'S' (scale) or 'N' (no scaling).
  * @param[in] n         Order of the original state-space system, n >= 0.
  * @param[in] m         Number of system inputs, m >= 0.
  * @param[in] p         Number of system outputs, p >= 0.
  * @param[in,out] a     Double array, dimension (lda, n) or (n, lda).
  * On entry, the state matrix A.
  * On exit, the NR-by-NR state matrix Amin of a minimal realization.
  * @param[in] lda       Leading dimension of A. >= max(1,n).
  * @param[in,out] b     Double array, dimension (ldb, max(m,p)) or (n, ldb).
  * On entry, the N-by-M input matrix B. Workspace needed if p > m.
  * On exit, the NR-by-M input matrix Bmin.
  * @param[in] ldb       Leading dimension of B. >= max(1,n).
  * @param[in,out] c     Double array, dimension (ldc, n) or (max(m,p), ldc).
  * On entry, the P-by-N output matrix C. Workspace needed if m > p.
  * On exit, the P-by-NR output matrix Cmin.
  * @param[in] ldc       Leading dimension of C. >= max(1,m,p).
  * @param[in] d         Double array, dimension (ldd, max(m,p)) or (max(m,p), ldd).
  * The P-by-M direct transmission matrix D. Workspace needed.
  * @param[in] ldd       Leading dimension of D. >= max(1,m,p).
  * @param[out] nr       Order of the minimal state-space representation (Amin, Bmin, Cmin).
  * @param[out] index    Integer array, dimension (p) if leri='L', (m) if leri='R'.
  * Contains polynomial degrees for rows/columns of P(s).
  * @param[out] pcoeff   Double array, dimension (ldpco1, ldpco2, n+1).
  * Coefficients of the denominator matrix P(s).
  * Dimensions depend on leri: (p,p,n+1) if 'L', (m,m,n+1) if 'R'.
  * @param[in] ldpco1    First leading dimension of pcoeff. >= max(1,p) if 'L', >= max(1,m) if 'R'.
  * @param[in] ldpco2    Second leading dimension of pcoeff. >= max(1,p) if 'L', >= max(1,m) if 'R'.
  * @param[out] qcoeff   Double array, dimension (ldqco1, ldqco2, n+1).
  * Coefficients of the numerator matrix Q(s).
  * Dimensions depend on leri: (p,m,n+1) if 'L', (p,m,n+1) if 'R'.
  * @param[in] ldqco1    First leading dimension of qcoeff. >= max(1,p) if 'L', >= max(1,m,p) if 'R'.
  * @param[in] ldqco2    Second leading dimension of qcoeff. >= max(1,m) if 'L', >= max(1,m,p) if 'R'.
  * @param[out] vcoeff   Double array, dimension (ldvco1, ldvco2, n+1).
  * Coefficients of the intermediate matrix V(s).
  * Dimensions depend on leri: (p,nr,n+1) if 'L', (m,nr,n+1) if 'R'.
  * @param[in] ldvco1    First leading dimension of vcoeff. >= max(1,p) if 'L', >= max(1,m) if 'R'.
  * @param[in] ldvco2    Second leading dimension of vcoeff. >= max(1,n).
  * @param[in] tol       Tolerance for rank determination. <= 0 for default.
  * @param[in] row_major Integer flag: 0 for column-major, 1 for row-major.
  * Affects interpretation of lda, ldb, ldc, ldd, ldpco*, ldqco*, ldvco*.
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value
  * = 1: singular matrix encountered during V(s) computation
  * = 2: singular matrix encountered during P(s) computation
  * Memory allocation errors may also be returned.
  */
 SLICOT_EXPORT
 int slicot_tb03ad(char leri, char equil, int n, int m, int p,
                   double* a, int lda, double* b, int ldb,
                   double* c, int ldc, double* d, int ldd,
                   int* nr, int* index,
                   double* pcoeff, int ldpco1, int ldpco2,
                   double* qcoeff, int ldqco1, int ldqco2,
                   double* vcoeff, int ldvco1, int ldvco2,
                   double tol, int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* TB03AD_H */
 