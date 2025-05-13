/**
 * @file ab09bd.h
 * @brief C wrapper for SLICOT routine AB09BD
 *
 * This file provides a C interface to the SLICOT routine AB09BD,
 * which performs Singular Perturbation Approximation based model
 * reduction for stable systems.
 */

#ifndef AB09BD_H
#define AB09BD_H

#include "slicot_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Singular Perturbation Approximation based model reduction for stable systems.
 *
 * Computes a reduced order model (Ar,Br,Cr,Dr) for a stable original
 * state-space representation (A,B,C,D) by using either the square-root
 * or the balancing-free square-root Singular Perturbation Approximation (SPA)
 * model reduction method.
 *
 * @param[in] dico      Specifies the type of the original system:
 *                      = 'C': Continuous-time system.
 *                      = 'D': Discrete-time system.
 * @param[in] job       Specifies the model reduction approach:
 *                      = 'B': Use square-root SPA method.
 *                      = 'N': Use balancing-free square-root SPA method.
 * @param[in] equil     Specifies whether to equilibrate the system:
 *                      = 'S': Perform equilibration (scaling).
 *                      = 'N': Do not perform equilibration.
 * @param[in] ordsel    Specifies the order selection method:
 *                      = 'F': The resulting order NR is fixed.
 *                      = 'A': The resulting order NR is automatically determined.
 * @param[in] n         The order of the original state-space representation, n >= 0.
 * @param[in] m         The number of system inputs, m >= 0.
 * @param[in] p         The number of system outputs, p >= 0.
 * @param[in,out] nr    On entry with ORDSEL = 'F', desired order of reduced system, 0<=nr<=n.
 *                      On exit, the order of the resulting reduced order model.
 * @param[in,out] a     Double array, dimension (lda, n) or (n, lda).
 *                      On entry, the state dynamics matrix A.
 *                      On exit, the leading nr-by-nr part contains the state matrix Ar.
 * @param[in] lda       The leading dimension of array A.
 * @param[in,out] b     Double array, dimension (ldb, m) or (n, ldb).
 *                      On entry, the input/state matrix B.
 *                      On exit, the leading nr-by-m part contains the matrix Br.
 * @param[in] ldb       The leading dimension of array B.
 * @param[in,out] c     Double array, dimension (ldc, n) or (p, ldc).
 *                      On entry, the state/output matrix C.
 *                      On exit, the leading p-by-nr part contains the matrix Cr.
 * @param[in] ldc       The leading dimension of array C.
 * @param[in,out] d     Double array, dimension (ldd, m) or (p, ldd).
 *                      On entry and exit, the direct transmission matrix D.
 * @param[in] ldd       The leading dimension of array D.
 * @param[out] hsv      Double array, dimension (n). Contains the Hankel singular values.
 * @param[in] tol1      Tolerance for determining the order of reduced system when ORDSEL = 'A'.
 *                      Recommended: tol1 = c*HNORM(A,B,C), where c is in [0.00001,0.001].
 *                      Default if tol1 <= 0: N*EPS*HNORM(A,B,C).
 * @param[in] tol2      Tolerance for determining order of minimal realization.
 *                      Recommended: tol2 = N*EPS*HNORM(A,B,C).
 *                      Default if tol2 <= 0: N*EPS*HNORM(A,B,C).
 *                      If tol2 > 0, then tol2 <= tol1.
 * @param[out] iwarn    Warning indicator: 
 *                      0 = no warning, 
 *                      1 = selected order nr is greater than minimal realization order.
 * @param[in] row_major Integer flag:
 *                      = 0: Arrays a, b, c, d are column-major (Fortran style).
 *                      = 1: Arrays a, b, c, d are row-major (C style).
 *
 * @return info         Error indicator:
 *                      = 0: successful exit
 *                      < 0: if info = -i, the i-th argument had an illegal value.
 *                      = 1: the reduction of A to real Schur form failed.
 *                      = 2: the state matrix A is not stable.
 *                      = 3: the computation of Hankel singular values failed.
 */
SLICOT_EXPORT
int slicot_ab09bd(
    char dico, char job, char equil, char ordsel,
    int n, int m, int p, int* nr,
    double* a, int lda, double* b, int ldb,
    double* c, int ldc, double* d, int ldd,
    double* hsv, double tol1, double tol2, int* iwarn,
    int row_major
);

#ifdef __cplusplus
}
#endif

#endif /* AB09BD_H */
