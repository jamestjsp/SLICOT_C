/**
 * @file ib01rd.h
 * @brief Header for C wrapper of SLICOT routine IB01RD.
 */

#ifndef SLICOT_WRAPPER_IB01RD_H
#define SLICOT_WRAPPER_IB01RD_H

#include "slicot_utils.h" // Provides SLICOT_EXPORT macro

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Estimates the initial state of a linear time-invariant (LTI)
 * discrete-time system, given the system matrices (A,B,C,D) and
 * the input and output trajectories of the system.
 * @details This is a C wrapper for the SLICOT Fortran routine IB01RD.
 * Workspace is allocated internally. Matrix A is assumed to be in a
 * real Schur form.
 *
 * @param job [in] CHARACTER*1. Specifies whether or not the matrix D is zero:
 * = 'Z': the matrix D is zero;
 * = 'N': the matrix D is not zero.
 * @param n [in] INTEGER. The order of the system. N >= 0.
 * @param m [in] INTEGER. The number of system inputs. M >= 0.
 * @param l [in] INTEGER. The number of system outputs. L > 0.
 * @param nsmp [in] INTEGER. The number of rows of matrices U and Y (number of samples). NSMP >= N.
 * @param a [in] const DOUBLE PRECISION array, dimension (LDA, N) or (N, LDA) if row_major.
 * The system state matrix A in a real Schur form.
 * If N=0, this array is not referenced.
 * @param lda [in] INTEGER. Leading dimension of A.
 * If row_major=0, LDA >= MAX(1,N).
 * If row_major=1, LDA >= MAX(1,N) (cols of A).
 * @param b [in] const DOUBLE PRECISION array, dimension (LDB, M) or (N, LDB) if row_major.
 * The system input matrix B.
 * If N=0 or M=0, this array is not referenced.
 * @param ldb [in] INTEGER. Leading dimension of B.
 * If row_major=0, LDB >= N if N>0,M>0; else LDB >= 1.
 * If row_major=1, LDB >= M if N>0,M>0; else LDB >= 1 (cols of B).
 * @param c [in] const DOUBLE PRECISION array, dimension (LDC, N) or (L, LDC) if row_major.
 * The system output matrix C.
 * @param ldc [in] INTEGER. Leading dimension of C.
 * If row_major=0, LDC >= L.
 * If row_major=1, LDC >= N (if N>0), else LDC >= 1 (cols of C).
 * @param d [in] const DOUBLE PRECISION array, dimension (LDD, M) or (L, LDD) if row_major.
 * The system input-output matrix D.
 * If M=0 or JOB='Z', this array is not referenced.
 * @param ldd [in] INTEGER. Leading dimension of D.
 * If row_major=0, LDD >= L if M>0,JOB='N'; else LDD >= 1.
 * If row_major=1, LDD >= M if M>0,JOB='N'; else LDD >= 1 (cols of D).
 * @param u [in] const DOUBLE PRECISION array, dimension (LDU, M) or (NSMP, LDU) if row_major.
 * The NSMP-by-M input-data sequence matrix U.
 * If M=0, this array is not referenced.
 * @param ldu [in] INTEGER. Leading dimension of U.
 * If row_major=0, LDU >= MAX(1,NSMP) if M>0; else LDU >= 1.
 * If row_major=1, LDU >= M if M>0; else LDU >= 1 (cols of U).
 * @param y [in] const DOUBLE PRECISION array, dimension (LDY, L) or (NSMP, LDY) if row_major.
 * The NSMP-by-L output-data sequence matrix Y.
 * @param ldy [in] INTEGER. Leading dimension of Y.
 * If row_major=0, LDY >= MAX(1,NSMP).
 * If row_major=1, LDY >= L (cols of Y).
 * @param x0 [out] DOUBLE PRECISION array, dimension (N).
 * The estimated initial state of the system, x(0).
 * @param tol [in] DOUBLE PRECISION. Tolerance for estimating rank. If TOL <= 0, EPS is used. TOL <= 1.
 * @param iwarn [out] INTEGER pointer. Warning indicator:
 * = 0: no warning;
 * = 4: least squares problem rank-deficient.
 * Can be NULL if not needed.
 * @param dwork_info [out] DOUBLE PRECISION array, dimension (2), optional.
 * If not NULL and INFO=0:
 * dwork_info[0] = optimal LDWORK.
 * dwork_info[1] = reciprocal condition number.
 * If not NULL and INFO=-22:
 * dwork_info[0] = minimum LDWORK.
 * Can be NULL if not needed.
 * @param row_major [in] INTEGER. Specifies matrix storage:
 * 0 for column-major (Fortran default),
 * 1 for row-major (C default).
 *
 * @return info Error indicator:
 * = 0:  successful exit;
 * < 0:  if INFO = -i, the i-th argument had an illegal value;
 * = 2:  SVD algorithm did not converge.
 * = SLICOT_MEMORY_ERROR (-1010): internal memory allocation failed.
 */
SLICOT_EXPORT
int slicot_ib01rd(
    char job, int n, int m, int l, int nsmp,
    const double* a, int lda,
    const double* b, int ldb,
    const double* c, int ldc,
    const double* d, int ldd,
    const double* u, int ldu,
    const double* y, int ldy,
    double* x0,
    double tol,
    int* iwarn,
    double* dwork_info,
    int row_major);

#ifdef __cplusplus
}
#endif

#endif /* SLICOT_WRAPPER_IB01RD_H */
