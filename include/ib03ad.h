/**
 * @file ib03ad.h
 * @brief Header for C wrapper of SLICOT routine IB03AD.
 */

#ifndef SLICOT_WRAPPER_IB03AD_H
#define SLICOT_WRAPPER_IB03AD_H

#include "slicot_utils.h" // Provides SLICOT_EXPORT macro

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Estimates parameters of a Wiener system using Levenberg-Marquardt.
 * @details This is a C wrapper for the SLICOT Fortran routine IB03AD.
 * Workspace is allocated internally.
 *
 * @param init_char [in] CHARACTER. Specifies initialization:
 * 'L': initialize linear part only.
 * 'S': initialize static nonlinearity only.
 * 'B': initialize both linear and nonlinear parts.
 * 'N': do not initialize anything.
 * @param alg_char [in] CHARACTER. Algorithm for linear systems:
 * 'D': direct Cholesky based.
 * 'I': iterative Conjugate Gradients.
 * @param stor_char [in] CHARACTER. Storage for J'*J if ALG='D':
 * 'F': full storage.
 * 'P': packed storage.
 * @param nobr [in] INTEGER. Number of block rows for Hankel matrices if INIT='L' or 'B'. NOBR > 0.
 * @param m [in] INTEGER. Number of system inputs. M >= 0.
 * @param l [in] INTEGER. Number of system outputs. L >= 0. (L > 0 if INIT='L' or 'B').
 * @param nsmp [in] INTEGER. Number of input/output samples. NSMP >= 0.
 * (NSMP >= 2*(M+L+1)*NOBR - 1 if INIT='L' or 'B').
 * @param n_ptr [in, out] INTEGER pointer. Order of the linear part.
 * If INIT='L'/'B' and *n_ptr < 0 on entry, order is estimated.
 * If INIT='S'/'N', *n_ptr >= 0 on entry.
 * If INIT='L'/'B' and *n_ptr >=0 on entry, then *n_ptr < NOBR and *n_ptr != 0.
 * On exit, contains the used/estimated order.
 * @param nn [in] INTEGER. Number of neurons for nonlinear part. NN >= 0.
 * @param itmax1 [in] INTEGER. Max iterations for nonlinearity initialization. ITMAX1 >= 0 if INIT='S'/'B'.
 * @param itmax2 [in] INTEGER. Max iterations for overall optimization. ITMAX2 >= 0.
 * @param nprint [in] INTEGER. Controls printing of iterates.
 * @param u [in] const DOUBLE PRECISION array. Input samples, NSMP-by-M.
 * Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param ldu [in] INTEGER. Leading dimension of U.
 * If row_major=0, LDU >= MAX(1,NSMP) (if M>0).
 * If row_major=1, LDU >= M (if M>0).
 * @param y [in] const DOUBLE PRECISION array. Output samples, NSMP-by-L.
 * Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param ldy [in] INTEGER. Leading dimension of Y.
 * If row_major=0, LDY >= MAX(1,NSMP) (if L>0).
 * If row_major=1, LDY >= L (if L>0).
 * @param x [in, out] DOUBLE PRECISION array, dimension (*lx_ptr). System parameters.
 * Structure depends on INIT, NN, L, N, M. See SLICOT docs.
 * @param lx_ptr [in, out] INTEGER pointer. Length of array X.
 * On entry, must be large enough for parameters.
 * On exit, if N was estimated, *lx_ptr might be updated if initial was too small.
 * @param tol1 [in] DOUBLE PRECISION. Tolerance for nonlinearity initialization sum of squares.
 * @param tol2 [in] DOUBLE PRECISION. Tolerance for overall optimization sum of squares.
 * @param out_iwork_summary [out] INTEGER array, dimension (3), optional (can be NULL).
 * If not NULL and call is successful/warning:
 * out_iwork_summary[0] = total function evaluations.
 * out_iwork_summary[1] = total Jacobian evaluations.
 * out_iwork_summary[2] = number of RCOND estimates in DWORK.
 * @param out_dwork_summary [out] DOUBLE PRECISION array, dimension (at least 10), optional (can be NULL).
 * If not NULL and call is successful/warning:
 * out_dwork_summary[0] = optimal/minimum LDWORK.
 * out_dwork_summary[1] = residual error norm.
 * out_dwork_summary[2] = iterations performed (main).
 * out_dwork_summary[3] = CG iterations performed (main).
 * out_dwork_summary[4] = final Levenberg factor (main).
 * If INIT='S'/'B':
 * out_dwork_summary[5-9] = corresponding values for initialization step.
 * @param iwarn_ptr [out] INTEGER pointer. Warning indicator. Can be NULL.
 * @param row_major [in] INTEGER. Specifies matrix storage for U, Y:
 * 0 for column-major (Fortran default),
 * 1 for row-major (C default).
 *
 * @return info Error indicator:
 * = 0:  successful exit;
 * < 0:  if INFO = -i, the i-th argument had an illegal value;
 * > 0:  Fortran routine specific error (see SLICOT IB03AD documentation).
 * = SLICOT_MEMORY_ERROR (-1010): internal memory allocation failed.
 */
SLICOT_EXPORT
int slicot_ib03ad(
    char init_char, char alg_char, char stor_char,
    int nobr, int m, int l, int nsmp,
    int* n_ptr,
    int nn,
    int itmax1, int itmax2, int nprint,
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

#endif /* SLICOT_WRAPPER_IB03AD_H */
