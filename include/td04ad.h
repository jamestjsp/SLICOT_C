/**
 * @file td04ad.h
 * @brief C wrapper for SLICOT routine TD04AD
 *
 * This file provides a C interface to the SLICOT routine TD04AD,
 * which finds a minimal state-space representation (A,B,C,D) for a
 * proper transfer matrix T(s) given as row or column polynomial
 * vectors over denominator polynomials.
 */

 #ifndef TD04AD_H
 #define TD04AD_H
 
 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Finds a minimal state-space representation for a transfer matrix fraction.
  *
  * Computes (A,B,C,D) equivalent to T(s) given as U(s)/D(s) where D(s)
  * contains denominator polynomials for rows ('R') or columns ('C').
  *
  * @param[in] rowcol    Specifies input format: 'R' (rows over common den.) or 'C' (cols over common den.).
  * @param[in] m         Number of system inputs, m >= 0.
  * @param[in] p         Number of system outputs, p >= 0.
  * @param[in] index     Integer array, dimension (p) if rowcol='R', (m) if rowcol='C'.
  * Contains denominator polynomial degrees.
  * @param[in] dcoeff    Double array, dimension (lddcoe, kdcoef). kdcoef=max(index)+1.
  * Coefficients of denominator polynomials D(s). Size p x kdcoef or m x kdcoef.
  * @param[in] lddcoe    Leading dimension of dcoeff. >= max(1,p) if 'R', >= max(1,m) if 'C'.
  * @param[in] ucoeff    Double array, dimension (lduco1, lduco2, kdcoef).
  * Coefficients of the numerator matrix U(s).
  * Dimensions depend on rowcol: (p,m,kdcoef) if 'R', (p,m,kdcoef) if 'C'.
  * Modified if rowcol='C', but restored on exit.
  * @param[in] lduco1    First leading dimension of ucoeff. >= max(1,p) if 'R', >= max(1,m,p) if 'C'.
  * @param[in] lduco2    Second leading dimension of ucoeff. >= max(1,m) if 'R', >= max(1,m,p) if 'C'.
  * @param[out] nr       Order of the resulting minimal realization.
  * @param[out] a        Double array, dimension (lda, nr) or (nr, lda). State matrix A.
  * @param[in] lda       Leading dimension of A. >= max(1,nr). Note: Max possible N = sum(index) needed for allocation if row_major.
  * @param[out] b        Double array, dimension (ldb, max(m,p)) or (nr, ldb). Input matrix B. Workspace needed.
  * @param[in] ldb       Leading dimension of B. >= max(1,nr). Note: Max possible N needed for allocation if row_major.
  * @param[out] c        Double array, dimension (ldc, nr) or (max(m,p), ldc). Output matrix C. Workspace needed.
  * @param[in] ldc       Leading dimension of C. >= max(1,m,p).
  * @param[out] d        Double array, dimension (ldd, max(m,p)) or (max(m,p), ldd). Direct transmission matrix D. Workspace needed.
  * @param[in] ldd       Leading dimension of D. >= max(1,p) if 'R', >= max(1,m,p) if 'C'.
  * @param[in] tol       Tolerance for rank determination. <= 0 for default.
  * @param[in] row_major Integer flag: 0 for column-major, 1 for row-major.
  * Affects interpretation of lda, ldb, ldc, ldd, lddcoe, lduco*.
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value
  * > 0: if info = i, the i-th leading coefficient in DCOEFF is near zero (potential overflow).
  * Memory allocation errors may also be returned.
  */
 SLICOT_C_WRAPPER_API
 int slicot_td04ad(char rowcol, int m, int p, const int* index,
                   const double* dcoeff, int lddcoe,
                   const double* ucoeff, int lduco1, int lduco2,
                   int* nr,
                   double* a, int lda, double* b, int ldb,
                   double* c, int ldc, double* d, int ldd,
                   double tol, int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* TD04AD_H */
 