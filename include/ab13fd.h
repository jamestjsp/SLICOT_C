/**
 * @file ab13fd.h
 * @brief C wrapper for SLICOT routine AB13FD
 *
 * This file provides a C interface to the SLICOT routine AB13FD,
 * which computes the distance from a real matrix A to the nearest
 * complex matrix with an eigenvalue on the imaginary axis, using SVD.
 */

#ifndef AB13FD_H
#define AB13FD_H

 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Computes distance from matrix A to nearest matrix with eigenvalue on imaginary axis (SVD method).
 *
 * Computes beta(A), the 2-norm distance from a real matrix A to the
 * nearest complex matrix with an eigenvalue on the imaginary axis.
 * If A is stable, beta(A) is the complex stability radius.
 * Also returns the frequency OMEGA where the minimum occurs.
 *
 * @param[in] n         The order of the matrix A, n >= 0.
 * @param[in] a         Double array for matrix A. Dimension (lda, n) or (n, lda).
 * @param[in] lda       The leading dimension of array A. >= max(1,n).
 * @param[out] beta     The computed distance beta(A) (an upper bound).
 * @param[out] omega    The frequency w where the minimum singular value occurs.
 * @param[in] tol       Tolerance for accuracy. If tol < eps, eps is used.
 * @param[in] row_major Integer flag:
 * = 0: Array a is column-major (Fortran style).
 * = 1: Array a is row-major (C style).
 *
 * @return info         Error indicator:
 * = 0: successful exit
 * < 0: if info = -i, the i-th argument had an illegal value.
 * = 1: Algorithm failed to converge within tolerance. BETA is upper bound.
 * = 2: QR or SVD algorithm failed to converge.
 * Memory allocation errors may also be returned.
 */
SLICOT_EXPORT
int slicot_ab13fd(int n, const double* a, int lda,
                  double* beta, double* omega, double tol,
                  int row_major);

#ifdef __cplusplus
}
#endif

#endif /* AB13FD_H */
