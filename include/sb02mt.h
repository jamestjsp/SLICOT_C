/**
 * @file sb02mt.h
 * @brief C wrapper for SLICOT routine SB02MT
 *
 * This file provides a C interface to the SLICOT routine SB02MT,
 * which converts linear-quadratic optimal control problems with
 * coupling weighting terms (matrix L) to standard problems by
 * computing modified matrices A, Q and optionally G = B*inv(R)*B'.
 */

 #ifndef SB02MT_H
 #define SB02MT_H
 
 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Converts optimal control problems with coupling terms to standard form.
  *
  * Computes G = B*inv(R)*B', A_new = A - B*inv(R)*L', Q_new = Q - L*inv(R)*L'.
  * Handles factored or unfactored R (symmetric, possibly indefinite).
  * Overwrites input matrices A, Q, and potentially B, R, L, IPIV depending on flags.
  *
  * @param[in] jobg      Specifies whether to compute G:
  * = 'G': Compute G.
  * = 'N': Do not compute G.
  * @param[in] jobl      Specifies if L is zero:
  * = 'Z': L is zero (A and Q are not modified).
  * = 'N': L is nonzero.
  * @param[in] fact      Specifies how R is provided:
  * = 'N': R contains the matrix R.
  * = 'C': R contains the Cholesky factor of R (R must be positive definite).
  * = 'U': R contains factors U/L and D from symmetric indefinite factorization. IPIV must be provided.
  * @param[in] uplo      Specifies which triangle of R, Q, G is stored:
  * = 'U': Upper triangle.
  * = 'L': Lower triangle.
  * @param[in] n         The order of matrices A, Q, G; rows of B, L. n >= 0.
  * @param[in] m         The order of matrix R; columns of B, L. m >= 0.
  * @param[in,out] a     Double array, dimension (lda, n) or (n, lda).
  * On entry (if jobl='N'), the matrix A.
  * On exit (if jobl='N'), the modified matrix A_new. Not referenced if jobl='Z'.
  * @param[in] lda       Leading dimension of A. >=max(1,n) if jobl='N', else >=1.
  * @param[in,out] b     Double array, dimension (ldb, m) or (n, ldb).
  * On entry, the matrix B.
  * On exit (if oufact=1), contains B * inv(chol(R)').
  * @param[in] ldb       Leading dimension of B. >= max(1,n).
  * @param[in,out] q     Double array, dimension (ldq, n) or (n, ldq).
  * On entry (if jobl='N'), the symmetric matrix Q (triangle specified by uplo).
  * On exit (if jobl='N'), the modified symmetric matrix Q_new (triangle specified by uplo).
  * Not referenced if jobl='Z'.
  * @param[in] ldq       Leading dimension of Q. >=max(1,n) if jobl='N', else >=1.
  * @param[in,out] r     Double array, dimension (ldr, m) or (m, ldr).
  * On entry, contains R or its factors based on 'fact'.
  * On exit (if oufact=1 or 2), contains computed factors (Cholesky or UdU'/LdL').
  * @param[in] ldr       Leading dimension of R. >= max(1,m).
  * @param[in,out] l     Double array, dimension (ldl, m) or (n, ldl).
  * On entry (if jobl='N'), the matrix L.
  * On exit (if jobl='N' and oufact=1), contains L * inv(chol(R)').
  * Not referenced if jobl='Z'.
  * @param[in] ldl       Leading dimension of L. >=max(1,n) if jobl='N', else >=1.
  * @param[in,out] ipiv  Integer array, dimension (m).
  * On entry (if fact='U'), pivot information from DSYTRF.
  * On exit (if oufact=2), pivot information from DSYTRF.
  * Not referenced if fact='C'. Used as workspace if fact='N'.
  * @param[out] oufact   Information about factorization used/computed:
  * = 0: M=0, no factorization used.
  * = 1: Cholesky factorization used/computed.
  * = 2: UdU'/LdL' factorization used/computed.
  * @param[out] g        Double array, dimension (ldg, n) or (n, ldg).
  * If jobg='G', contains the computed symmetric matrix G (triangle specified by uplo).
  * Not referenced if jobg='N'.
  * @param[in] ldg       Leading dimension of G. >=max(1,n) if jobg='G', else >=1.
  * @param[in] row_major Integer flag:
  * = 0: Arrays a, b, q, r, l, g are column-major (Fortran style).
  * = 1: Arrays a, b, q, r, l, g are row-major (C style).
  * (Array ipiv is 1D).
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * = i: (1 <= i <= M) if D(i,i) is exactly zero in UdU'/LdL' factorization.
  * = M+1: if R is numerically singular.
  * Memory allocation errors may also be returned.
  */
 SLICOT_EXPORT
 int slicot_sb02mt(char jobg, char jobl, char fact, char uplo,
                   int n, int m,
                   double* a, int lda, double* b, int ldb,
                   double* q, int ldq, double* r, int ldr,
                   double* l, int ldl, int* ipiv, int* oufact,
                   double* g, int ldg, int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* SB02MT_H */
 