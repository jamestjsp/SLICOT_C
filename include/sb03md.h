/**
 * @file sb03md.h
 * @brief C wrapper for SLICOT routine SB03MD
 *
 * This file provides a C interface to the SLICOT routine SB03MD,
 * which solves continuous- or discrete-time Lyapunov equations and
 * optionally estimates the separation condition number.
 */

 #ifndef SB03MD_H
 #define SB03MD_H
 
 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Solves real Lyapunov equations and/or estimates separation.
  *
  * Solves op(A)'*X + X*op(A) = scale*C (continuous) or
  * op(A)'*X*op(A) - X = scale*C (discrete), where op(A) = A or A'.
  * Optionally estimates separation sep(op(A), -op(A)') or sep(op(A), op(A)').
  * Optionally computes forward error bound FERR.
  * Uses Schur factorization of A.
  *
  * @param[in] dico      Specifies the type of Lyapunov equation:
  * = 'C': Continuous-time.
  * = 'D': Discrete-time.
  * @param[in] job       Specifies the computation to perform:
  * = 'X': Compute solution X only.
  * = 'S': Compute separation SEP only.
  * = 'B': Compute both solution X and separation SEP (and FERR).
  * @param[in] fact      Specifies if Schur factorization of A is provided:
  * = 'F': A and U contain Schur factors on entry.
  * = 'N': Compute Schur factorization; A and U are overwritten.
  * @param[in] trana     Specifies the form of op(A):
  * = 'N': op(A) = A.
  * = 'T': op(A) = A'.
  * = 'C': op(A) = A'.
  * @param[in] n         The order of the matrices A, X, C, n >= 0.
  * @param[in,out] a     Double array, dimension (lda, n) or (n, lda).
  * On entry, the matrix A. If fact='F', contains upper quasi-triangular Schur factor S.
  * On exit (if fact='N'), contains the computed Schur factor S.
  * @param[in] lda       The leading dimension of array A. >= max(1,n).
  * @param[in,out] u     Double array, dimension (ldu, n) or (n, ldu).
  * If fact='F', contains the orthogonal Schur matrix U on entry.
  * If fact='N', contains the computed orthogonal Schur matrix U on exit.
  * @param[in] ldu       The leading dimension of array U. >= max(1,n).
  * @param[in,out] c     Double array, dimension (ldc, n) or (n, ldc).
  * On entry (if job='X' or 'B'), the symmetric matrix C.
  * On exit (if job='X' or 'B'), the symmetric solution matrix X.
  * Not referenced if job='S'.
  * @param[in] ldc       Leading dimension of C. >=1 if job='S', >=max(1,n) otherwise.
  * @param[out] scale    Scale factor (<= 1) applied to C to avoid overflow in X.
  * @param[out] sep      Estimated separation (if job='S' or 'B'). Not referenced if job='X'.
  * @param[out] ferr     Estimated forward error bound for X (if job='B'). Not referenced if job='X' or 'S'.
  * @param[out] wr       Double array, dimension (n). Real parts of eigenvalues of A (if fact='N'). Not referenced if fact='F'.
  * @param[out] wi       Double array, dimension (n). Imaginary parts of eigenvalues of A (if fact='N'). Not referenced if fact='F'.
  * @param[in] row_major Integer flag:
  * = 0: Arrays a, u, c are column-major (Fortran style).
  * = 1: Arrays a, u, c are row-major (C style).
  * (Arrays wr, wi are 1D).
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * > 0: (if fact='N') if info = i, QR algorithm failed to compute eigenvalues.
  * = N+1: Equation is singular or nearly singular; perturbed values used.
  * Memory allocation errors may also be returned.
  */
 SLICOT_C_WRAPPER_API
 int slicot_sb03md(char dico, char job, char fact, char trana, int n,
                   double* a, int lda, double* u, int ldu,
                   double* c, int ldc, double* scale, double* sep, double* ferr,
                   double* wr, double* wi, int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* SB03MD_H */
 