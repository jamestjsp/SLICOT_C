/**
 * @file mb03rd.h
 * @brief C wrapper for SLICOT routine MB03RD
 *
 * This file provides a C interface to the SLICOT routine MB03RD,
 * which reduces a matrix A in real Schur form to block-diagonal form
 * using well-conditioned non-orthogonal similarity transformations.
 */

#ifndef MB03RD_H
#define MB03RD_H

 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Reduces a real Schur form matrix to block-diagonal form.
 *
 * Reduces a matrix A (already in real Schur form) to block-diagonal form
 * using well-conditioned non-orthogonal similarity transformations X, such
 * that inv(X)*A*X is block diagonal. The condition number of elementary
 * transformations is bounded by PMAX. Optionally accumulates transformations
 * into a given matrix. Optionally reorders eigenvalues first.
 *
 * @param[in] jobx      Specifies whether to accumulate transformations in X:
 * = 'N': Do not accumulate transformations (X is not referenced).
 * = 'U': Accumulate transformations in X (X must be provided on entry).
 * @param[in] sort      Specifies eigenvalue sorting/strategy:
 * = 'N': Do not reorder; use "closest to mean" strategy.
 * = 'S': Reorder clusters; use "closest to mean" strategy.
 * = 'C': Do not reorder; use "closest-neighbour" strategy.
 * = 'B': Reorder clusters; use "closest-neighbour" strategy.
 * @param[in] n         The order of the matrix A, n >= 0.
 * @param[in] pmax      Upper bound on infinity norm of transformation submatrices, pmax >= 1.0.
 * @param[in,out] a     Double array, dimension (lda, n) or (n, lda).
 * On entry, the matrix A in real Schur form.
 * On exit, the block-diagonalized matrix (non-diagonal blocks zeroed).
 * @param[in] lda       The leading dimension of array A. >= max(1,n).
 * @param[in,out] x     Double array, dimension (ldx, n) or (n, ldx).
 * If jobx='U', contains input matrix on entry, updated by transformations on exit (X_new = X_old * T).
 * If jobx='N', not referenced (can be NULL, ldx ignored).
 * @param[in] ldx       The leading dimension of array X. >= 1 if jobx='N', >= max(1,n) if jobx='U'.
 * @param[out] nblcks   The number of diagonal blocks in the resulting A.
 * @param[out] blsize   Integer array, dimension (n). The first nblcks elements contain block sizes.
 * @param[out] wr       Double array, dimension (n). Real parts of eigenvalues of A.
 * @param[out] wi       Double array, dimension (n). Imaginary parts of eigenvalues of A.
 * @param[in] tol       Tolerance for eigenvalue clustering (if sort='S' or 'B').
 * If tol > 0, absolute tolerance. If tol < 0, relative tolerance.
 * If tol = 0, default relative tolerance sqrt(sqrt(eps)) is used.
 * Ignored if sort='N' or 'C'.
 * @param[in] row_major Integer flag:
 * = 0: Arrays a, x are column-major (Fortran style).
 * = 1: Arrays a, x are row-major (C style).
 * (Arrays blsize, wr, wi are 1D).
 *
 * @return info         Error indicator:
 * = 0: successful exit
 * < 0: if info = -i, the i-th argument had an illegal value.
 * Memory allocation errors may also be returned.
 */
SLICOT_EXPORT
int slicot_mb03rd(char jobx, char sort, int n, double pmax,
                  double* a, int lda, double* x, int ldx,
                  int* nblcks, int* blsize, double* wr, double* wi,
                  double tol, int row_major);

#ifdef __cplusplus
}
#endif

#endif /* MB03RD_H */
