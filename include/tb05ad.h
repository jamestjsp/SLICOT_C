/**
 * @file tb05ad.h
 * @brief C wrapper for SLICOT routine TB05AD
 *
 * This file provides a C interface to the SLICOT routine TB05AD,
 * which computes the frequency response matrix G(freq) = C*(freq*I - A)^-1 * B
 * for a given state-space representation (A,B,C) at a specified complex frequency.
 */

 #ifndef TB05AD_H
 #define TB05AD_H
 
 #include <stddef.h> // For size_t
 #include "slicot_utils.h" // For slicot_complex_double definition
 
 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Computes the frequency response matrix G(freq) = C*(freq*I - A)^-1 * B.
  *
  * Calculates the complex frequency response matrix G at a given complex frequency freq.
  * Can optionally balance A, compute eigenvalues, and estimate condition number.
  *
  * @param[in] baleig    Specifies balancing/eigenvalue/condition number options:
  * 'N': None.
  * 'C': Condition number estimate only.
  * 'B'/'E': Balance A and compute eigenvalues (requires inita='G').
  * 'A': Balance A, compute eigenvalues and condition number (requires inita='G').
  * @param[in] inita     Specifies if A is already in Hessenberg form:
  * 'G': General matrix (first call).
  * 'H': Upper Hessenberg (subsequent calls if A not modified).
  * @param[in] n         Order of matrix A, n >= 0.
  * @param[in] m         Number of inputs (columns of B), m >= 0.
  * @param[in] p         Number of outputs (rows of C), p >= 0.
  * @param[in] freq      Complex frequency at which to evaluate G(freq).
  * @param[in,out] a     Double array, dimension (lda, n) or (n, lda).
  * On entry, the state matrix A.
  * On exit, if inita='G', contains the upper Hessenberg form of A.
  * @param[in] lda       Leading dimension of A. >= max(1,n).
  * @param[in,out] b     Double array, dimension (ldb, m) or (n, ldb).
  * On entry, the input matrix B.
  * On exit, if inita='G', contains the transformed B.
  * @param[in] ldb       Leading dimension of B. >= max(1,n).
  * @param[in,out] c     Double array, dimension (ldc, n) or (p, ldc).
  * On entry, the output matrix C.
  * On exit, if inita='G', contains the transformed C.
  * @param[in] ldc       Leading dimension of C. >= max(1,p).
  * @param[out] rcond    Estimate of the reciprocal condition number of (freq*I - A) if baleig='C' or 'A'.
  * @param[out] g        Complex array, dimension (ldg, m) or (p, ldg). The frequency response matrix G(freq).
  * @param[in] ldg       Leading dimension of G. >= max(1,p).
  * @param[out] evre     Double array, dimension (n). Real parts of eigenvalues of A (if requested by baleig).
  * @param[out] evim     Double array, dimension (n). Imaginary parts of eigenvalues of A (if requested by baleig).
  * @param[out] hinvb    Complex array, dimension (ldhinv, m) or (n, ldhinv). Contains (freq*I - A)^-1 * B.
  * @param[in] ldhinv    Leading dimension of HINVB. >= max(1,n).
  * @param[in] row_major Integer flag: 0 for column-major, 1 for row-major.
  * Affects interpretation of lda, ldb, ldc, ldg, ldhinv.
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value
  * = 1: eigenvalue computation did not converge (results might still be useful)
  * = 2: freq is too near an eigenvalue, or matrix is ill-conditioned (RCOND too small)
  * Memory allocation errors may also be returned.
  */
 int slicot_tb05ad(char baleig, char inita, int n, int m, int p,
                   slicot_complex_double freq,
                   double* a, int lda, double* b, int ldb,
                   double* c, int ldc, double* rcond,
                   slicot_complex_double* g, int ldg,
                   double* evre, double* evim,
                   slicot_complex_double* hinvb, int ldhinv,
                   int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* TB05AD_H */
 