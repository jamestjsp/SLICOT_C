/**
 * @file ib01cd.h
 * @brief Header for C wrapper of SLICOT routine IB01CD.
 */

#ifndef SLICOT_WRAPPER_IB01CD_H
#define SLICOT_WRAPPER_IB01CD_H

#include "slicot_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Estimate the initial state and system matrices B and D using A, B, and input-output data.
 * @details This routine estimates the initial state and, optionally, the system matrices B and D
 * of a linear time-invariant discrete-time system, given the system matrices (A,B,C,D), or 
 * (when B and D are estimated) only the matrix pair (A,C), and the input and output trajectories.
 * This is a C wrapper for the SLICOT Fortran routine IB01CD.
 * **Workspace is allocated internally.**
 * **Note on Zero Dimensions:** Behavior for zero dimensions (e.g., N=0, M=0) aligns with the 
 * underlying Fortran routine. Many SLICOT routines handle zero dimensions gracefully with INFO=0.
 *
 * @param jobx0 [in] Specifies whether to compute initial state:
 *        'X': compute the initial state x(0)
 *        'N': do not compute the initial state
 * @param comuse [in] Specifies whether system matrices B and D should be computed or used:
 *        'C': compute the system matrices B and D
 *        'U': use the system matrices B and D
 *        'N': do not compute/use the matrices B and D
 * @param job [in] Specifies which system matrices to compute/use (if COMUSE='C' or 'U'):
 *        'B': compute/use matrix B only (D is zero)
 *        'D': compute/use matrices B and D
 * @param n [in] The order of the system (N >= 0).
 * @param m [in] The number of system inputs (M >= 0).
 * @param l [in] The number of system outputs (L > 0).
 * @param nsmp [in] The number of samples (NSMP >= 0, with specific constraints based on other parameters).
 * @param a [in] System state matrix A (N x N). Can be NULL if N=0 or not required by JOBX0/COMUSE.
 *        Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param lda [in] Leading dimension of array A.
 *        If row_major=0: lda >= max(1,N). If row_major=1: lda >= max(1,N).
 * @param b [in,out] System input matrix B (N x M). Can be NULL if N=0 or M=0 or not required.
 *        Input if COMUSE='U', output if COMUSE='C'.
 * @param ldb [in] Leading dimension of array B.
 *        If row_major=0: ldb >= max(1,N). If row_major=1: ldb >= max(1,M).
 * @param c [in] System output matrix C (L x N). Can be NULL if N=0 or not required.
 * @param ldc [in] Leading dimension of array C.
 *        If row_major=0: ldc >= L. If row_major=1: ldc >= max(1,N).
 * @param d [in,out] System input-output matrix D (L x M). Can be NULL if M=0 or JOB='B'.
 *        Input if COMUSE='U' and JOB='D', output if COMUSE='C' and JOB='D'.
 * @param ldd [in] Leading dimension of array D.
 *        If row_major=0: ldd >= L. If row_major=1: ldd >= max(1,M).
 * @param u [in,out] Input data matrix U (NSMP x M). Can be NULL if M=0 or NSMP=0.
 *        May be modified if COMUSE='C' and JOB='D'.
 * @param ldu [in] Leading dimension of array U.
 *        If row_major=0: ldu >= max(1,NSMP). If row_major=1: ldu >= max(1,M).
 * @param y [in] Output data matrix Y (NSMP x L). Can be NULL if NSMP=0 or not required.
 * @param ldy [in] Leading dimension of array Y.
 *        If row_major=0: ldy >= max(1,NSMP). If row_major=1: ldy >= max(1,L).
 * @param x0 [out] Estimated initial state vector (length N). Can be NULL if N=0 or JOBX0='N'.
 * @param v [out] Orthogonal matrix V from real Schur factorization (N x N). Can be NULL if N=0.
 * @param ldv [in] Leading dimension of array V.
 *        If row_major=0: ldv >= max(1,N). If row_major=1: ldv >= max(1,N).
 * @param tol [in] Tolerance for rank estimation (TOL <= 1). If TOL <= 0, machine precision is used.
 * @param iwarn [out] Warning indicator (can be NULL if not needed):
 *        = 0: no warning
 *        = 4: rank-deficient coefficient matrix
 *        = 6: matrix A is unstable
 * @param row_major [in] Matrix storage format:
 *        0: column-major (Fortran style)
 *        1: row-major (C style)
 *
 * @return info Error indicator:
 *         = 0: successful exit
 *         < 0: if info = -i, the i-th argument had an illegal value
 *         = 1: QR algorithm failed to compute eigenvalues of matrix A
 *         = 2: SVD algorithm did not converge
 *         = SLICOT_MEMORY_ERROR (-1010): internal memory allocation failed
 */
SLICOT_EXPORT
int slicot_ib01cd(const char* jobx0, const char* comuse, const char* job,
                  int n, int m, int l, int nsmp,
                  const double* a, int lda,
                  double* b, int ldb,
                  const double* c, int ldc,
                  double* d, int ldd,
                  const double* u, int ldu,
                  const double* y, int ldy,
                  double* x0,
                  double* v, int ldv,
                  double tol,
                  int* iwarn,
                  int row_major);

#ifdef __cplusplus
}
#endif

#endif /* SLICOT_WRAPPER_IB01CD_H */
