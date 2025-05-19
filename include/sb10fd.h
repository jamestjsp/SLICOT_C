/**
 * @file sb10fd.h
 * @brief Header for C wrapper of SLICOT routine SB10FD.
 */

#ifndef SB10FD_H
#define SB10FD_H

#include "slicot_utils.h" // Provides SLICOT_EXPORT macro

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Computes the matrices of an H-infinity (sub)optimal n-state controller.
 * @details This is a C wrapper for the SLICOT Fortran routine SB10FD.
 * Workspace is allocated internally.
 *
 * @param n      (input) The order of the system. N >= 0.
 * @param m      (input) The column size of the matrix B. M >= 0.
 * @param np     (input) The row size of the matrix C. NP >= 0.
 * @param ncon   (input) The number of control inputs (M2). M >= NCON >= 0, NP-NMEAS >= NCON.
 * @param nmeas  (input) The number of measurements (NP2). NP >= NMEAS >= 0, M-NCON >= NMEAS.
 * @param gamma_val (input) The value of gamma. GAMMA >= 0.
 * @param a      (input) System state matrix A. Dimensions (N x N).
 * Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param lda    Leading dimension of the C array storing A.
 * If row_major=0, lda >= max(1, N).
 * If row_major=1, lda >= max(1, N) (number of columns).
 * @param b      (input) System input matrix B. Dimensions (N x M).
 * @param ldb    Leading dimension of B. If row_major=0, ldb >= max(1,N). If row_major=1, ldb >= max(1,M).
 * @param c      (input) System output matrix C. Dimensions (NP x N).
 * @param ldc    Leading dimension of C. If row_major=0, ldc >= max(1,NP). If row_major=1, ldc >= max(1,N).
 * @param d      (input) System input/output matrix D. Dimensions (NP x M).
 * @param ldd    Leading dimension of D. If row_major=0, ldd >= max(1,NP). If row_major=1, ldd >= max(1,M).
 * @param[out] ak Controller state matrix AK. Dimensions (N x N).
 * @param ldak   Leading dimension of AK. If row_major=0, ldak >= max(1,N). If row_major=1, ldak >= max(1,N).
 * @param[out] bk Controller input matrix BK. Dimensions (N x NMEAS).
 * @param ldbk   Leading dimension of BK. If row_major=0, ldbk >= max(1,N). If row_major=1, ldbk >= max(1,NMEAS).
 * @param[out] ck Controller output matrix CK. Dimensions (NCON x N).
 * @param ldck   Leading dimension of CK. If row_major=0, ldck >= max(1,NCON). If row_major=1, ldck >= max(1,N).
 * @param[out] dk Controller input/output matrix DK. Dimensions (NCON x NMEAS).
 * @param lddk   Leading dimension of DK. If row_major=0, lddk >= max(1,NCON). If row_major=1, lddk >= max(1,NMEAS).
 * @param[out] rcond Reciprocal condition numbers (array of size 4).
 * @param tol    (input) Tolerance for rank determination. If tol <= 0, sqrt(machine_precision) is used.
 * @param row_major Specifies matrix storage for A, B, C, D, AK, BK, CK, DK:
 * 0 for column-major (Fortran default),
 * 1 for row-major (C default).
 *
 * @return info Error indicator:
 * = 0:  successful exit
 * < 0:  if info = -i, the i-th argument had an illegal value (wrapper or Fortran validation)
 * (Note: argument indices in this wrapper may differ from Fortran)
 * = 1:  if the matrix | A-j*omega*I  B2  | had not full column rank
 * |    C1        D12 |
 * = 2:  if the matrix | A-j*omega*I  B1  |  had not full row rank
 * |    C2        D21 |
 * = 3:  if the matrix D12 had not full column rank
 * = 4:  if the matrix D21 had not full row rank
 * = 5:  if SVD algorithm did not converge
 * = 6:  if the controller is not admissible (gamma too small)
 * = 7:  if X-Riccati equation was not solved successfully
 * = 8:  if Y-Riccati equation was not solved successfully
 * = 9:  if the determinant of Im2 + Tu*D11HAT*Ty*D22 is zero
 * = SLICOT_MEMORY_ERROR (-1010): internal memory allocation failed.
 */
SLICOT_EXPORT
int slicot_sb10fd(int n, int m, int np, int ncon, int nmeas,
                  double gamma_val,
                  double* a, int lda,
                  double* b, int ldb,
                  double* c, int ldc,
                  double* d, int ldd,
                  double* ak, int ldak,
                  double* bk, int ldbk,
                  double* ck, int ldck,
                  double* dk, int lddk,
                  double* rcond, double tol,
                  int row_major);

#ifdef __cplusplus
}
#endif

#endif /* SB10FD_H */
