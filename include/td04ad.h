/**
 * @file td04ad.h
 * @brief Header for C wrapper of SLICOT routine TD04AD.
 * @details This routine computes a minimal state-space representation (A,B,C,D)
 * for a proper transfer matrix T(s) given as either row or column
 * polynomial vectors over denominator polynomials.
 */

#ifndef SLICOT_WRAPPER_TD04AD_H
#define SLICOT_WRAPPER_TD04AD_H

#include "slicot_utils.h" // Provides SLICOT_EXPORT macro

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Computes a minimal state-space representation for a proper transfer matrix.
 * @details This is a C wrapper for the SLICOT Fortran routine TD04AD.
 * Workspace (IWORK, DWORK) is allocated internally by this wrapper.
 *
 * @param rowcol Specifies how T(s) is given:
 * 'R': T(s) as rows over common denominators.
 * 'C': T(s) as columns over common denominators.
 * @param m The number of system inputs. m >= 0.
 * @param p The number of system outputs. p >= 0.
 * @param index Input INTEGER array.
 * If rowcol='R', dimension p. index[i] is degree of i-th row's denominator.
 * If rowcol='C', dimension m. index[j] is degree of j-th col's denominator.
 * All elements of index must be >= 0.
 * @param dcoeff Input DOUBLE PRECISION array of denominator coefficients.
 * Let kdcoef = max(index[k]) + 1.
 * If rowcol='R', dimension (lddcoe, kdcoef). dcoeff[i*lddcoe_f + k] (Fortran)
 * is coeff of s**(index[i]-k) of i-th denom.
 * If rowcol='C', dimension (lddcoe, kdcoef). dcoeff[j*lddcoe_f + k] (Fortran)
 * is coeff of s**(index[j]-k) of j-th denom.
 * Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param lddcoe Leading dimension of the C array storing dcoeff.
 * If row_major=0 (column-major):
 * If rowcol='R', lddcoe >= max(1, p).
 * If rowcol='C', lddcoe >= max(1, m).
 * If row_major=1 (row-major):
 * If rowcol='R', lddcoe >= max(1, kdcoef) (number of columns for p rows).
 * If rowcol='C', lddcoe >= max(1, kdcoef) (number of columns for m rows).
 * @param ucoeff Input DOUBLE PRECISION array of numerator coefficients.
 * Let kdcoef = max(index[k]) + 1.
 * If rowcol='R', conceptual dimensions P x M x kdcoef.
 * ucoeff[i][j][k] is coeff of s**(index[i]-k) of U(i,j)(s).
 * If rowcol='C', conceptual dimensions P x M x kdcoef for input part.
 * The Fortran routine uses a MAX(M,P) x MAX(M,P) x kdcoef array,
 * where the PxM part is input, and rest is workspace.
 * ucoeff[i][j][k] is coeff of s**(index[j]-k) of U(i,j)(s).
 * If rowcol='C', this array is modified by Fortran but restored on exit.
 * Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param lduco1 First leading dimension of the C array storing ucoeff.
 * If row_major=0 (column-major):
 * If rowcol='R', lduco1 >= max(1, p).
 * If rowcol='C', lduco1 >= max(1, p, m).
 * If row_major=1 (row-major):
 * If rowcol='R', lduco1 >= max(1, m) (cols of P-planes, assuming PxMxK layout).
 * If rowcol='C', lduco1 >= max(1, m) (cols of P-planes, assuming PxMxK layout for the PxM input part).
 * @param lduco2 Second leading dimension of the C array storing ucoeff.
 * If row_major=0 (column-major):
 * If rowcol='R', lduco2 >= max(1, m).
 * If rowcol='C', lduco2 >= max(1, p, m).
 * If row_major=1 (row-major):
 * If rowcol='R', lduco2 >= max(1, kdcoef) (cols of M-vectors).
 * If rowcol='C', lduco2 >= max(1, kdcoef) (cols of M-vectors).
 * @param tol Tolerance for rank determination. If tol <= 0, a default value is used.
 * @param[out] nr The order of the resulting minimal realization.
 * @param[out] a Output DOUBLE PRECISION array, state dynamics matrix (NR x NR).
 * Caller must allocate space for N_sum x N_sum, where N_sum = sum(index).
 * Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param lda Leading dimension of the C array storing A.
 * If row_major=0, lda >= max(1, N_sum).
 * If row_major=1, lda >= max(1, N_sum) (number of columns).
 * @param[out] b Output DOUBLE PRECISION array, input matrix (NR x M).
 * Caller must allocate space for N_sum x M.
 * Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param ldb Leading dimension of the C array storing B.
 * If row_major=0, ldb >= max(1, N_sum).
 * If row_major=1, ldb >= max(1, M) (number of columns).
 * @param[out] c Output DOUBLE PRECISION array, output matrix (P x NR).
 * Caller must allocate space for P x N_sum.
 * Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param ldc Leading dimension of the C array storing C.
 * If row_major=0, ldc >= max(1, p).
 * If row_major=1, ldc >= max(1, N_sum) (number of columns).
 * @param[out] d Output DOUBLE PRECISION array, direct transmission matrix (P x M).
 * Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param ldd Leading dimension of the C array storing D.
 * If row_major=0, ldd >= max(1, p).
 * If row_major=1, ldd >= max(1, m) (number of columns).
 * @param row_major Specifies matrix storage for input/output arrays:
 * 0 for column-major (Fortran default).
 * 1 for row-major (C default).
 *
 * @return info Error indicator:
 * = 0: successful exit.
 * < 0: if info = -i, the i-th argument had an illegal value.
 * Fortran argument indices:
 * 1:ROWCOL, 2:M, 3:P, 4:INDEX, 5:DCOEFF, 6:LDDCOE,
 * 7:UCOEFF, 8:LDUCO1, 9:LDUCO2, 10:NR, 11:A, 12:LDA,
 * 13:B, 14:LDB, 15:C, 16:LDC, 17:D, 18:LDD, 19:TOL,
 * 20:IWORK, 21:DWORK, 22:LDWORK.
 * > 0: if info = i, then i is the first 1-based index in DCOEFF
 * for which abs(DCOEFF(i,1)) is too small (leading coefficient
 * of a denominator polynomial is nearly zero).
 * = SLICOT_MEMORY_ERROR (-1010): internal memory allocation failed.
 */
SLICOT_EXPORT
int slicot_td04ad(char rowcol, int m, int p, const int* index,
                  const double* dcoeff, int lddcoe,
                  const double* ucoeff, int lduco1, int lduco2,
                  double tol,
                  int* nr, double* a, int lda,
                  double* b, int ldb, double* c, int ldc,
                  double* d, int ldd, int row_major);

#ifdef __cplusplus
}
#endif

#endif // SLICOT_WRAPPER_TD04AD_H
