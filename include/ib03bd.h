/**
 * @file ib03bd.h
 * @brief Header for C wrapper of SLICOT routine IB03BD.
 */

#ifndef SLICOT_WRAPPER_IB03BD_H
#define SLICOT_WRAPPER_IB03BD_H

#include "slicot_utils.h" // Provides SLICOT_EXPORT macro

#ifdef __cplusplus
extern "C" {
#endif

// Used for sizing the out_dwork_summary array in the caller if they want RCOND values.
// This is an estimate; the actual number of RCOND values is returned in out_iwork_summary[2].
#define MAX_EXPECTED_RCONDS_IN_DWORK_HDR_IB03BD 30 

/**
 * @brief Wiener system identification using a MINPACK-like Levenberg-Marquardt algorithm.
 * @details This is a C wrapper for the SLICOT Fortran routine IB03BD.
 * Workspace is allocated internally.
 *
 * @param init_char [in] CHARACTER. Specifies initialization:
 * 'L': initialize linear part only.
 * 'S': initialize static nonlinearity only.
 * 'B': initialize both linear and nonlinear parts.
 * 'N': do not initialize anything.
 * @param nobr_in [in] INTEGER. Number of block rows for Hankel matrices if INIT='L' or 'B'. NOBR > 0.
 * @param m_in [in] INTEGER. Number of system inputs. M >= 0.
 * @param l_in [in] INTEGER. Number of system outputs. L >= 0. (L > 0 if INIT='L' or 'B').
 * @param nsmp_in [in] INTEGER. Number of input/output samples. NSMP >= 0.
 * (NSMP >= 2*(M+L+1)*NOBR - 1 if INIT='L' or 'B').
 * @param n_ptr [in, out] INTEGER pointer. Order of the linear part.
 * If INIT='L'/'B' and *n_ptr < 0 on entry, order is estimated by the routine.
 * If INIT='S'/'N', *n_ptr >= 0 on entry and is used as the fixed order.
 * If INIT='L'/'B' and *n_ptr >=0 on entry, then 0 < *n_ptr < NOBR.
 * On exit, contains the used/estimated order.
 * @param nn_in [in] INTEGER. Number of neurons for nonlinear part. NN >= 0.
 * @param itmax1_in [in] INTEGER. Max iterations for nonlinearity initialization. ITMAX1 >= 0 if INIT='S'/'B'.
 * @param itmax2_in [in] INTEGER. Max iterations for overall optimization. ITMAX2 >= 0.
 * @param nprint_in [in] INTEGER. Controls printing of iterates from Fortran. >0 enables printing.
 * @param u [in] const DOUBLE PRECISION array. Input samples, NSMP-by-M.
 * Stored column-wise if row_major=0, row-wise if row_major=1.
 * Not referenced if M=0.
 * @param ldu [in] INTEGER. Leading dimension of U.
 * If row_major=0, LDU >= MAX(1,NSMP) (if M>0).
 * If row_major=1, LDU >= M (if M>0). Else LDU >= 1.
 * @param y [in] const DOUBLE PRECISION array. Output samples, NSMP-by-L.
 * Stored column-wise if row_major=0, row-wise if row_major=1.
 * Not referenced if L=0 (only possible if INIT='S' or 'N').
 * @param ldy [in] INTEGER. Leading dimension of Y.
 * If row_major=0, LDY >= MAX(1,NSMP) (if L>0).
 * If row_major=1, LDY >= L (if L>0). Else LDY >= 1.
 * @param x [in, out] DOUBLE PRECISION array, dimension (*lx_ptr). System parameters.
 * Structure depends on INIT, NN, L, N, M. See SLICOT docs for IB03BD.
 * If INIT='B', contents on input are not used.
 * If INIT='L', 'S', or 'N', appropriate parts must be pre-filled.
 * On exit, contains the optimized parameters.
 * @param lx_ptr [in, out] INTEGER pointer. Length of array X.
 * On entry, must be large enough for parameters.
 * On exit, if N was estimated (*n_ptr < 0 on entry) and input *lx_ptr was too small,
 * *lx_ptr is updated to the required length (and INFO=-21 is returned by Fortran).
 * @param tol1 [in] DOUBLE PRECISION. Tolerance for nonlinearity initialization. If <0, SQRT(EPS) is used.
 * @param tol2 [in] DOUBLE PRECISION. Tolerance for overall optimization. If <0, SQRT(EPS) is used.
 * @param out_iwork_summary [out] INTEGER array, dimension (3), optional (can be NULL).
 * If not NULL and call is successful/warning:
 * out_iwork_summary[0] = total function evaluations.
 * out_iwork_summary[1] = total Jacobian evaluations.
 * out_iwork_summary[2] = number of RCOND estimates in DWORK (relevant if INIT='L'/'B').
 * @param out_dwork_summary [out] DOUBLE PRECISION array, dimension (at least 8, or 8 + number of RCONDs), optional (can be NULL).
 * Sized as 8 + MAX_EXPECTED_RCONDS_IN_DWORK_HDR_IB03BD for safety by caller.
 * If not NULL and call is successful/warning:
 * out_dwork_summary[0] = optimal/minimum LDWORK.
 * out_dwork_summary[1] = residual error norm (main optimization).
 * out_dwork_summary[2] = iterations performed (main).
 * out_dwork_summary[3] = final Levenberg factor (main).
 * If INIT='S'/'B':
 * out_dwork_summary[4-7] = corresponding values for initialization step. (Note: Fortran DWORK(5-8))
 * If INIT='L'/'B' and out_iwork_summary[2] > 0, elements from index 8 onwards
 * would contain RCOND estimates.
 * @param iwarn_ptr [out] INTEGER pointer. Warning indicator from Fortran. Can be NULL.
 * @param row_major [in] INTEGER. Specifies matrix storage for U, Y:
 * 0 for column-major (Fortran default),
 * 1 for row-major (C default).
 *
 * @return info Error indicator:
 * = 0:  successful exit;
 * < 0:  if INFO = -i, the i-th argument had an illegal value (C wrapper or Fortran);
 * > 0:  Fortran routine specific error (see SLICOT IB03BD documentation).
 * = SLICOT_MEMORY_ERROR (-1010): internal memory allocation failed in C wrapper.
 */
SLICOT_EXPORT
int slicot_ib03bd(
    char init_char, int nobr_in, int m_in, int l_in, int nsmp_in,
    int* n_ptr,
    int nn_in,
    int itmax1_in, int itmax2_in, int nprint_in,
    const double* u, int ldu,
    const double* y, int ldy,
    double* x, int* lx_ptr,
    double tol1, double tol2,
    int* out_iwork_summary,
    double* out_dwork_summary,
    int* iwarn_ptr,
    int row_major);

#ifdef __cplusplus
}
#endif

#endif /* SLICOT_WRAPPER_IB03BD_H */
