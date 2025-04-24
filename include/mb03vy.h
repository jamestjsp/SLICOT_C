/**
 * @file mb03vy.h
 * @brief C wrapper for SLICOT routine MB03VY
 *
 * This file provides a C interface to the SLICOT routine MB03VY,
 * which generates the orthogonal matrices Q_j from the elementary
 * reflectors computed by MB03VD.
 */

#ifndef MB03VY_H
#define MB03VY_H

 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Generates orthogonal matrices from reflectors computed by MB03VD.
 *
 * Generates the real orthogonal matrices Q_1, Q_2, ..., Q_p,
 * which are defined as the product of ihi-ilo elementary reflectors
 * of order n, as returned by SLICOT Library routine MB03VD.
 *
 * @param[in] n         The order of the matrices Q_j, n >= 0.
 * @param[in] p         The number of transformation matrices, p >= 1.
 * @param[in] ilo       Lower row/column index used in MB03VD.
 * @param[in] ihi       Upper row/column index used in MB03VD.
 * 1 <= ilo <= max(1,n); min(ilo,n) <= ihi <= n.
 * @param[in,out] a     Double array, dimension (lda1, lda2, p) stored contiguously
 * in **column-major (Fortran-style)** order.
 * On entry, the strictly lower triangular parts contain reflector
 * data from MB03VD.
 * On exit, contains the computed orthogonal matrices Q_1 to Q_p.
 * @param[in] lda1      The first leading dimension of A. >= max(1,n).
 * @param[in] lda2      The second leading dimension of A. >= max(1,n).
 * @param[in] tau       Double array, dimension (ldtau, p). Contains scalar factors
 * of the elementary reflectors from MB03VD.
 * @param[in] ldtau     The leading dimension of TAU. >= max(1,n-1).
 * @param[in] row_major This flag is ignored for the 3D array 'a'. The data for 'a'
 * **must** be provided in column-major (Fortran-style) layout.
 *
 * @return info         Error indicator:
 * = 0: successful exit
 * < 0: if info = -i, the i-th argument had an illegal value.
 * Memory allocation errors may also be returned by the wrapper.
 */
SLICOT_EXPORT
int slicot_mb03vy(int n, int p, int ilo, int ihi,
                  double* a, int lda1, int lda2,
                  const double* tau, int ldtau,
                  int row_major); // row_major ignored for 'a'

#ifdef __cplusplus
}
#endif

#endif /* MB03VY_H */
