/**
 * @file sb10hd.h
 * @brief Header for C wrapper of SLICOT routine SB10HD.
 */

#ifndef SB10HD_H
#define SB10HD_H

#include "slicot_utils.h" // Provides SLICOT_EXPORT macro

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Computes the matrices of the H2 optimal n-state controller.
 * @details This is a C wrapper for the SLICOT Fortran routine SB10HD.
 * Workspace is allocated internally. Assumes D11 block of D is zero.
 *
 * @param n_param (input) The order of the system. N >= 0.
 * @param m_param (input) The column size of the matrix B. M >= 0.
 * @param np_param (input) The row size of the matrix C. NP >= 0.
 * @param ncon_param (input) The number of control inputs (M2). M >= NCON >= 0, NP-NMEAS >= NCON.
 * @param nmeas_param (input) The number of measurements (NP2). NP >= NMEAS >= 0, M-NCON >= NMEAS.
 * @param a      (input) System state matrix A. Dimensions (N x N).
 * Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param lda    Leading dimension of A. If row_major=0, lda >= max(1,N). If row_major=1, lda >= max(1,N) (cols).
 * @param b      (input) System input matrix B. Dimensions (N x M).
 * @param ldb    Leading dimension of B. If row_major=0, ldb >= max(1,N). If row_major=1, ldb >= max(1,M) (cols).
 * @param c      (input) System output matrix C. Dimensions (NP x N).
 * @param ldc    Leading dimension of C. If row_major=0, ldc >= max(1,NP). If row_major=1, ldc >= max(1,N) (cols).
 * @param d      (input) System input/output matrix D. Dimensions (NP x M). Assumed D11 = 0.
 * @param ldd    Leading dimension of D. If row_major=0, ldd >= max(1,NP). If row_major=1, ldd >= max(1,M) (cols).
 * @param[out] ak Controller state matrix AK. Dimensions (N x N).
 * @param ldak   Leading dimension of AK. If row_major=0, ldak >= max(1,N). If row_major=1, ldak >= max(1,N) (cols).
 * @param[out] bk Controller input matrix BK. Dimensions (N x NMEAS).
 * @param ldbk   Leading dimension of BK. If row_major=0, ldbk >= max(1,N). If row_major=1, ldbk >= max(1,NMEAS) (cols).
 * @param[out] ck Controller output matrix CK. Dimensions (NCON x N).
 * @param ldck   Leading dimension of CK. If row_major=0, ldck >= max(1,NCON). If row_major=1, ldck >= max(1,N) (cols).
 * @param[out] dk Controller input/output matrix DK. Dimensions (NCON x NMEAS).
 * @param lddk   Leading dimension of DK. If row_major=0, lddk >= max(1,NCON). If row_major=1, lddk >= max(1,NMEAS) (cols).
 * @param[out] rcond Reciprocal condition numbers (array of size 4).
 * @param tol    (input) Tolerance. If tol <= 0, sqrt(machine_precision) is used.
 * @param row_major Specifies matrix storage: 0 for column-major, 1 for row-major.
 *
 * @return info Error indicator:
 * = 0:  successful exit
 * < 0:  if info = -i, the i-th argument had an illegal value
 * = 1:  if D12 had not full column rank
 * = 2:  if D21 had not full row rank
 * = 3:  if SVD algorithm did not converge
 * = 4:  if X-Riccati equation was not solved successfully
 * = 5:  if Y-Riccati equation was not solved successfully
 * = SLICOT_MEMORY_ERROR (-1010): internal memory allocation failed.
 */
SLICOT_EXPORT
int slicot_sb10hd(int n_param, int m_param, int np_param, int ncon_param, int nmeas_param,
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

#endif /* SB10HD_H */
