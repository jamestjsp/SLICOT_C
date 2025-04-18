/**
 * @file mb05nd.h
 * @brief C wrapper for SLICOT routine MB05ND
 *
 * This file provides a C interface to the SLICOT routine MB05ND,
 * which computes the matrix exponential exp(A*delta) and its integral.
 */

#ifndef MB05ND_H
#define MB05ND_H

 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Computes the matrix exponential exp(A*delta) and its integral.
 *
 * Calculates F(delta) = exp(A*delta) and H(delta) = integral of F(s) ds
 * from s=0 to s=delta, using Pade approximations.
 *
 * @param[in] n         The order of the matrix A, n >= 0.
 * @param[in] delta     The scalar delta.
 * @param[in] a         Double array for matrix A. Dimension (lda, n) or (n, lda).
 * (Not needed if delta = 0).
 * @param[in] lda       The leading dimension of array A. >= max(1,n).
 * @param[out] ex       Double array, dimension (ldex, n) or (n, ldex).
 * Contains the computed exp(A*delta).
 * @param[in] ldex      The leading dimension of array EX. >= max(1,n).
 * @param[out] exint    Double array, dimension (ldexin, n) or (n, ldexin).
 * Contains the computed integral H(delta).
 * @param[in] ldexin    The leading dimension of array EXINT. >= max(1,n).
 * @param[in] tol       Tolerance for Pade approximation order. Recommended: sqrt(eps).
 * @param[in] row_major Integer flag:
 * = 0: Arrays a, ex, exint are column-major (Fortran style).
 * = 1: Arrays a, ex, exint are row-major (C style).
 *
 * @return info         Error indicator:
 * = 0: successful exit
 * < 0: if info = -i, the i-th argument had an illegal value.
 * = i: (1 <= i <= N) Pade denominator is singular at step i.
 * = N+1: delta * norm(A) may be too large (potential overflow).
 * Memory allocation errors may also be returned.
 */
SLICOT_C_WRAPPER_API
int slicot_mb05nd(int n, double delta, const double* a, int lda,
                  double* ex, int ldex, double* exint, int ldexin,
                  double tol, int row_major);

#ifdef __cplusplus
}
#endif

#endif /* MB05ND_H */
