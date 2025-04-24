/**
 * @file sb03od.h
 * @brief C wrapper for SLICOT routine SB03OD
 *
 * This file provides a C interface to the SLICOT routine SB03OD,
 * which solves stable continuous- or discrete-time Lyapunov equations
 * for the Cholesky factor U of the solution X = op(U)'*op(U).
 */

 #ifndef SB03OD_H
 #define SB03OD_H
 
 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Solves stable Lyapunov equations for the Cholesky factor of the solution.
  *
  * Solves op(A)'*X + X*op(A) = -scale^2 * op(B)'*op(B) (continuous) or
  * op(A)'*X*op(A) - X = -scale^2 * op(B)'*op(B) (discrete),
  * where X = op(U)'*op(U) and U is upper triangular.
  * A must be stable (continuous) or convergent (discrete).
  *
  * @param[in] dico      Specifies the type of Lyapunov equation:
  * = 'C': Continuous-time.
  * = 'D': Discrete-time.
  * @param[in] fact      Specifies if Schur factorization of A is provided:
  * = 'F': A and Q contain Schur factors on entry.
  * = 'N': Compute Schur factorization; A and Q are overwritten.
  * @param[in] trans     Specifies the form of op(K):
  * = 'N': op(K) = K.
  * = 'T': op(K) = K'.
  * @param[in] n         The order of matrix A; columns of op(B). n >= 0.
  * @param[in] m         The number of rows of op(B). m >= 0.
  * @param[in,out] a     Double array, dimension (lda, n) or (n, lda).
  * On entry, the matrix A. If fact='F', contains upper quasi-triangular Schur factor S.
  * On exit (if fact='N'), contains the computed Schur factor S.
  * @param[in] lda       Leading dimension of A. >= max(1,n).
  * @param[in,out] q     Double array, dimension (ldq, n) or (n, ldq).
  * If fact='F', contains the orthogonal Schur matrix Q on entry.
  * If fact='N', contains the computed orthogonal Schur matrix Q on exit.
  * @param[in] ldq       Leading dimension of Q. >= max(1,n).
  * @param[in,out] b     Double array. If trans='N', dimension (ldb, n) or (m, ldb).
  * If trans='T', dimension (ldb, m) or (n, ldb).
  * On entry, contains the matrix B.
  * On exit, the leading N-by-N part contains the upper triangular Cholesky factor U.
  * @param[in] ldb       Leading dimension of B.
  * If trans='N', ldb >= max(1,m,n). If trans='T', ldb >= max(1,n).
  * @param[out] scale    Scale factor (<= 1) applied to B to avoid overflow in X.
  * @param[out] wr       Double array, dimension (n). Real parts of eigenvalues of A.
  * @param[out] wi       Double array, dimension (n). Imaginary parts of eigenvalues of A.
  * @param[in] row_major Integer flag:
  * = 0: Arrays a, q, b are column-major (Fortran style).
  * = 1: Arrays a, q, b are row-major (C style).
  * (Arrays wr, wi are 1D).
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * = 1: Equation is nearly singular (warning). Perturbed values used.
  * = 2: (fact='N') A is not stable/convergent.
  * = 3: (fact='F') Schur factor S is not stable/convergent.
  * = 4: (fact='F') Schur factor S has blocks > 2x2.
  * = 5: (fact='F') Schur factor S has 2x2 block with real eigenvalues.
  * = 6: (fact='N') DGEES failed to converge.
  * Memory allocation errors may also be returned.
  */
 SLICOT_EXPORT
 int slicot_sb03od(char dico, char fact, char trans, int n, int m,
                   double* a, int lda, double* q, int ldq,
                   double* b, int ldb, double* scale,
                   double* wr, double* wi, int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* SB03OD_H */
 