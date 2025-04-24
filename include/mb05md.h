/**
 * @file mb05md.h
 * @brief C wrapper for SLICOT routine MB05MD
 *
 * This file provides a C interface to the SLICOT routine MB05MD,
 * which computes the matrix exponential exp(A*delta) for a real
 * non-defective matrix A.
 */

#ifndef MB05MD_H
#define MB05MD_H

 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Computes the matrix exponential exp(A*delta) for a non-defective matrix.
 *
 * Computes exp(A*delta) using eigenvalue decomposition. Also returns
 * eigenvalues and eigenvectors of A. Optionally performs balancing.
 * The input matrix A is overwritten with the result exp(A*delta).
 *
 * @param[in] balanc    Specifies balancing option:
 * = 'N': Do not diagonally scale.
 * = 'S': Perform diagonal scaling on A.
 * @param[in] n         The order of the matrix A, n >= 0.
 * @param[in] delta     The scalar delta.
 * @param[in,out] a     Double array, dimension (lda, n) or (n, lda).
 * On entry, the matrix A. On exit, exp(A*delta).
 * @param[in] lda       The leading dimension of array A. >= max(1,n).
 * @param[out] v        Double array, dimension (ldv, n) or (n, ldv). Eigenvector matrix V.
 * @param[in] ldv       The leading dimension of array V. >= max(1,n).
 * @param[out] y        Double array, dimension (ldy, n) or (n, ldy). Intermediate result.
 * If eigenvalues are real, Y = exp(Lambda*delta)*inv(V). Generally, exp(A*delta) = V*Y.
 * @param[in] ldy       The leading dimension of array Y. >= max(1,n).
 * @param[out] valr     Double array, dimension (n). Real parts of eigenvalues.
 * @param[out] vali     Double array, dimension (n). Imaginary parts of eigenvalues.
 * @param[in] row_major Integer flag:
 * = 0: Arrays a, v, y are column-major (Fortran style).
 * = 1: Arrays a, v, y are row-major (C style).
 * (Arrays valr, vali are 1D).
 *
 * @return info         Error indicator:
 * = 0: successful exit
 * < 0: if info = -i, the i-th argument had an illegal value.
 * = i: (1 <= i <= N) QR algorithm failed to compute all eigenvalues.
 * = N+1: Eigenvector matrix is singular; inverse could not be formed.
 * = N+2: Matrix A is defective (possibly due to rounding). Consider MB05ND/OD.
 * Memory allocation errors may also be returned.
 */
SLICOT_EXPORT
int slicot_mb05md(char balanc, int n, double delta,
                  double* a, int lda,
                  double* v, int ldv,
                  double* y, int ldy,
                  double* valr, double* vali,
                  int row_major);

#ifdef __cplusplus
}
#endif

#endif /* MB05MD_H */
