/**
 * @file tb04ad.h
 * @brief C wrapper for SLICOT routine TB04AD
 *
 * This file provides a C interface to the SLICOT routine TB04AD,
 * which finds the transfer matrix T(s) of a given state-space
 * representation (A,B,C,D), expressed as row or column polynomial
 * vectors over monic least common denominators.
 */

 #ifndef TB04AD_H
 #define TB04AD_H
 
 #include <stddef.h> // For size_t
 
 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Finds the transfer matrix T(s) = C(sI-A)^-1 B + D.
  *
  * Computes T(s) expressed as row or column polynomial vectors over
  * monic least common denominators.
  *
  * @param[in] rowcol    Specifies output format: 'R' (rows over common den.) or 'C' (cols over common den.).
  * @param[in] n         Order of the original state-space system, n >= 0.
  * @param[in] m         Number of system inputs, m >= 0.
  * @param[in] p         Number of system outputs, p >= 0.
  * @param[in,out] a     Double array, dimension (lda, n) or (n, lda).
  * On entry, the state matrix A.
  * On exit, the NR-by-NR state matrix of a transformed system.
  * @param[in] lda       Leading dimension of A. >= max(1,n).
  * @param[in,out] b     Double array, dimension (ldb, max(m,p)) or (n, ldb).
  * On entry, the N-by-M input matrix B. Workspace needed if rowcol='C' and p > m.
  * On exit, the NR-by-M transformed input matrix B.
  * @param[in] ldb       Leading dimension of B. >= max(1,n).
  * @param[in,out] c     Double array, dimension (ldc, n) or (max(m,p), ldc).
  * On entry, the P-by-N output matrix C. Workspace needed if rowcol='C' and m > p.
  * On exit, the P-by-NR transformed output matrix C.
  * @param[in] ldc       Leading dimension of C. >= max(1,p) if 'R', >= max(1,m,p) if 'C'.
  * @param[in] d         Double array, dimension (ldd, max(m,p)) or (max(m,p), ldd).
  * The P-by-M direct transmission matrix D. Workspace needed if rowcol='C'.
  * @param[in] ldd       Leading dimension of D. >= max(1,p) if 'R', >= max(1,m,p) if 'C'.
  * @param[out] nr       Order of the transformed state-space representation.
  * @param[out] index    Integer array, dimension (p) if rowcol='R', (m) if rowcol='C'.
  * Contains denominator polynomial degrees.
  * @param[out] dcoeff   Double array, dimension (lddcoe, n+1).
  * Coefficients of denominator polynomials. Size p x (max_deg+1) or m x (max_deg+1).
  * @param[in] lddcoe    Leading dimension of dcoeff. >= max(1,p) if 'R', >= max(1,m) if 'C'.
  * @param[out] ucoeff   Double array, dimension (lduco1, lduco2, n+1).
  * Coefficients of the numerator matrix U(s).
  * Dimensions depend on rowcol: (p,m,n+1) if 'R', (p,m,n+1) if 'C'.
  * @param[in] lduco1    First leading dimension of ucoeff. >= max(1,p) if 'R', >= max(1,m) if 'C'.
  * @param[in] lduco2    Second leading dimension of ucoeff. >= max(1,m) if 'R', >= max(1,p) if 'C'.
  * @param[in] tol1      Tolerance for determining row/column of T(s). <= 0 for default.
  * @param[in] tol2      Tolerance for rank determination (controllability/observability). <= 0 for default.
  * @param[in] row_major Integer flag: 0 for column-major, 1 for row-major.
  * Affects interpretation of lda, ldb, ldc, ldd, lddcoe, lduco*.
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value
  * Memory allocation errors may also be returned.
  */
 int slicot_tb04ad(char rowcol, int n, int m, int p,
                   double* a, int lda, double* b, int ldb,
                   double* c, int ldc, double* d, int ldd,
                   int* nr, int* index,
                   double* dcoeff, int lddcoe,
                   double* ucoeff, int lduco1, int lduco2,
                   double tol1, double tol2, int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* TB04AD_H */
 