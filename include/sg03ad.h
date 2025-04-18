/**
 * @file sg03ad.h
 * @brief C wrapper for SLICOT routine SG03AD
 *
 * This file provides a C interface to the SLICOT routine SG03AD,
 * which solves continuous- or discrete-time generalized Lyapunov
 * equations and optionally estimates the separation condition number.
 */

 #ifndef SG03AD_H
 #define SG03AD_H
 
 #include <stddef.h> // For size_t
 
 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Solves generalized Lyapunov equations and/or estimates separation.
  *
  * Solves op(A)'*X*op(E) + op(E)'*X*op(A) = scale*Y (continuous) or
  * op(A)'*X*op(A) - op(E)'*X*op(E) = scale*Y (discrete), where op(M) = M or M'.
  * Optionally estimates separation and forward error bound FERR.
  * Uses generalized Schur factorization of (A, E).
  *
  * @param[in] dico      Specifies the type of Lyapunov equation: 'C' or 'D'.
  * @param[in] job       Specifies the computation to perform: 'X', 'S', or 'B'.
  * @param[in] fact      Specifies if generalized Schur factorization is provided: 'F' or 'N'.
  * @param[in] trans     Specifies the form of op(A), op(E): 'N' or 'T'/'C'.
  * @param[in] uplo      Specifies which triangle of input Y is used: 'U' or 'L'.
  * @param[in] n         The order of matrices A, E, X, Y, n >= 0.
  * @param[in,out] a     Double array, dimension (lda, n) or (n, lda).
  * On entry, matrix A or Schur factor A_s. On exit (if fact='N'), Schur factor A_s.
  * @param[in] lda       Leading dimension of A. >= max(1,n).
  * @param[in,out] e     Double array, dimension (lde, n) or (n, lde).
  * On entry, matrix E or Schur factor E_s. On exit (if fact='N'), Schur factor E_s.
  * @param[in] lde       Leading dimension of E. >= max(1,n).
  * @param[in,out] q     Double array, dimension (ldq, n) or (n, ldq).
  * If fact='F', orthogonal matrix Q on entry. If fact='N', computed Q on exit.
  * @param[in] ldq       Leading dimension of Q. >= max(1,n).
  * @param[in,out] z     Double array, dimension (ldz, n) or (n, ldz).
  * If fact='F', orthogonal matrix Z on entry. If fact='N', computed Z on exit.
  * @param[in] ldz       Leading dimension of Z. >= max(1,n).
  * @param[in,out] x     Double array, dimension (ldx, n) or (n, ldx).
  * On entry (if job='X' or 'B'), the symmetric matrix Y (relevant triangle).
  * On exit (if job='X' or 'B'), the symmetric solution matrix X.
  * Not referenced if job='S'.
  * @param[in] ldx       Leading dimension of X. >=max(1,n) if job != 'S', else >=1.
  * @param[out] scale    Scale factor (<= 1) applied to Y to avoid overflow in X.
  * @param[out] sep      Estimated separation (if job='S' or 'B'). Not referenced if job='X'.
  * @param[out] ferr     Estimated forward error bound for X (if job='B'). Not referenced if job='X' or 'S'.
  * @param[out] alphar   Double array, dimension (n). Real parts of generalized eigenvalues (if fact='N').
  * @param[out] alphai   Double array, dimension (n). Imaginary parts of generalized eigenvalues (if fact='N').
  * @param[out] beta     Double array, dimension (n). Scaling factors of generalized eigenvalues (if fact='N').
  * @param[in] row_major Integer flag:
  * = 0: Arrays a, e, q, z, x are column-major.
  * = 1: Arrays a, e, q, z, x are row-major.
  * (Arrays alphar, alphai, beta are 1D).
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * = 1: (fact='F') Input A_s is not upper quasi-triangular.
  * = 2: (fact='N') QZ algorithm failed.
  * = 3: (dico='D') Pencil has reciprocal eigenvalues (singular equation). Perturbed solve.
  * = 4: (dico='C') Pencil has degenerate eigenvalues (singular equation). Perturbed solve.
  * Memory allocation errors may also be returned.
  */
 int slicot_sg03ad(char dico, char job, char fact, char trans, char uplo, int n,
                   double* a, int lda, double* e, int lde,
                   double* q, int ldq, double* z, int ldz,
                   double* x, int ldx, double* scale, double* sep, double* ferr,
                   double* alphar, double* alphai, double* beta, int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* SG03AD_H */
 