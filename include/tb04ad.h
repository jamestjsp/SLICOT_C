/**
 * @file tb04ad.h
 * @brief Header for C wrapper of SLICOT routine TB04AD.
 */

#ifndef SLICOT_WRAPPER_TB04AD_H
#define SLICOT_WRAPPER_TB04AD_H

#include "slicot_utils.h" // Provides SLICOT_EXPORT macro

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Converts a state-space system (A,B,C,D) to a transfer matrix T(s).
 * @details T(s) is expressed as either row or column polynomial vectors over monic least common denominator polynomials.
 * This is a C wrapper for the SLICOT Fortran routine TB04AD.
 * **Workspace is allocated internally.**
 *
 * @param rowcol Character indicating if T(s) is required as rows or columns over common denominators:
 * 'R': T(s) as rows over common denominators.
 * 'C': T(s) as columns over common denominators.
 * @param n      (input) The order of the state-space representation (order of A). n >= 0.
 * @param m      (input) The number of system inputs. m >= 0.
 * @param p      (input) The number of system outputs. p >= 0.
 * @param a      (input/output) Matrix A of the state-space system. Dimensions (n x n).
 * On entry, contains the original state dynamics matrix.
 * On exit, contains the transformed state dynamics matrix of order NR.
 * Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param lda    Leading dimension of the C array storing A.
 * If row_major=0 (column-major), lda >= max(1, n).
 * If row_major=1 (row-major), lda >= max(1, n) (number of columns).
 * @param b      (input/output) Matrix B of the state-space system. Dimensions (n x m).
 * On entry, contains the original input matrix.
 * On exit, contains the transformed input matrix of dimension (NR x m).
 * Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param ldb    Leading dimension of the C array storing B.
 * If row_major=0 (column-major), ldb >= max(1, n).
 * If row_major=1 (row-major), ldb >= max(1, m) (number of columns).
 * @param c      (input/output) Matrix C of the state-space system. Dimensions (p x n).
 * On entry, contains the original output matrix.
 * On exit, contains the transformed output matrix of dimension (p x NR).
 * Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param ldc    Leading dimension of the C array storing C.
 * If row_major=0 (column-major), ldc >= max(1, p) (if ROWCOL='R') or ldc >= max(1,m,p) (if ROWCOL='C').
 * If row_major=1 (row-major), ldc >= max(1, n) (number of columns).
 * @param d      (input) Matrix D of the state-space system. Dimensions (p x m).
 * Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param ldd    Leading dimension of the C array storing D.
 * If row_major=0 (column-major), ldd >= max(1, p) (if ROWCOL='R') or ldd >= max(1,m,p) (if ROWCOL='C').
 * If row_major=1 (row-major), ldd >= max(1, m) (number of columns).
 * @param[out] nr Pointer to store the order of the transformed state-space representation.
 * @param[out] index Integer array, dimension (porm), where porm = p if rowcol='R', else m.
 * Stores the degrees of the denominator polynomials. Can be NULL if porm is 0.
 * @param[out] dcoeff Double precision array, dimension (porm x (max_deg_val+1)).
 * Stores the coefficients of each denominator polynomial.
 * DCOEFF(i,k) is the coefficient in s**(INDEX(i)-k) of the i-th denominator polynomial.
 * Must be allocated by the caller with dimensions (lddcoe_c x (max_n_expected+1)). Can be NULL if porm is 0.
 * @param lddcoe_c Leading dimension of the C array dcoeff (number of rows, i.e., porm). Ignored if dcoeff_out is NULL.
 * @param[out] ucoeff Double precision array, dimension (porm x porp x (max_deg_val+1)).
 * Stores the coefficients of the numerator matrix U(s).
 * UCOEFF(i,j,k) is the coefficient in s**(INDEX(iorj)-k) of polynomial (i,j) of U(s).
 * Must be allocated by the caller. Can be NULL if porm or porp is 0.
 * @param lduco1_c First dimension of the C array ucoeff (porm). Ignored if ucoeff_out is NULL.
 * @param lduco2_c Second dimension of the C array ucoeff (porp). Ignored if ucoeff_out is NULL.
 * @param tol1   (input) Tolerance for determining transfer function coefficients. Default used if <= 0.
 * @param tol2   (input) Tolerance for separating controllable/observable subsystem. Default used if <= 0.
 * @param row_major Specifies matrix storage for A, B, C, D:
 * 0 for column-major (Fortran default),
 * 1 for row-major (C default).
 *
 * @return info Error indicator:
 * = 0:  successful exit
 * < 0:  if info = -i, the i-th argument had an illegal value (wrapper or Fortran validation)
 * = SLICOT_MEMORY_ERROR (-1010): internal memory allocation failed.
 */
SLICOT_EXPORT
int slicot_tb04ad(const char* rowcol, int n, int m, int p,
                  double* a, int lda,
                  double* b, int ldb,
                  double* c, int ldc,
                  const double* d, int ldd,
                  int* nr, int* index, double* dcoeff, int lddcoe_c,
                  double* ucoeff, int lduco1_c, int lduco2_c,
                  double tol1, double tol2,
                  int row_major);

#ifdef __cplusplus
}
#endif

#endif // SLICOT_WRAPPER_TB04AD_H
