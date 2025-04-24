/**
 * @file sg03bd.h
 * @brief C wrapper for SLICOT routine SG03BD
 *
 * This file provides a C interface to the SLICOT routine SG03BD,
 * which solves stable continuous- or discrete-time generalized
 * Lyapunov equations for the Cholesky factor U of the solution
 * X = op(U)'*op(U).
 */

 #ifndef SG03BD_H
 #define SG03BD_H
 
 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Solves stable generalized Lyapunov equations for the Cholesky factor.
  *
  * Solves op(A)'*X*op(E) + op(E)'*X*op(A) = -scale^2*op(B)'*op(B) (continuous) or
  * op(A)'*X*op(A) - op(E)'*X*op(E) = -scale^2*op(B)'*op(B) (discrete),
  * for the Cholesky factor U, where X = op(U)'*op(U) and U is upper triangular.
  * The pencil A - lambda*E must be stable (continuous) or convergent (discrete).
  *
  * @param[in] dico      Specifies the type of Lyapunov equation: 'C' or 'D'.
  * @param[in] fact      Specifies if generalized Schur factorization is provided: 'F' or 'N'.
  * @param[in] trans     Specifies the form of op(K): 'N' or 'T'.
  * @param[in] n         Order of matrices A, E; columns of op(B). n >= 0.
  * @param[in] m         Number of rows of op(B). m >= 0.
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
  * @param[in,out] b     Double array. If trans='N', dimension (ldb, n) or (m, ldb).
  * If trans='T', dimension (ldb, m) or (n, ldb).
  * On entry, contains the matrix B.
  * On exit, the leading N-by-N part contains the upper triangular Cholesky factor U.
  * @param[in] ldb       Leading dimension of B.
  * If trans='N', ldb >= max(1,m,n). If trans='T', ldb >= max(1,n).
  * @param[out] scale    Scale factor (<= 1) applied to B to avoid overflow in U.
  * @param[out] alphar   Double array, dimension (n). Real parts of generalized eigenvalues.
  * @param[out] alphai   Double array, dimension (n). Imaginary parts of generalized eigenvalues.
  * @param[out] beta     Double array, dimension (n). Scaling factors of generalized eigenvalues.
  * @param[in] row_major Integer flag:
  * = 0: Arrays a, e, q, z, b are column-major.
  * = 1: Arrays a, e, q, z, b are row-major.
  * (Arrays alphar, alphai, beta are 1D).
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * = 1: Pencil A - lambda*E is nearly singular. Perturbed solve.
  * = 2: (fact='F') Input A_s is not upper quasi-triangular.
  * = 3: (fact='F') A 2x2 block in pencil has non-complex eigenvalues.
  * = 4: (fact='N') QZ algorithm failed.
  * = 5: (dico='C') Pencil is not c-stable.
  * = 6: (dico='D') Pencil is not d-stable.
  * = 7: DSYEVX failed (discrete-time case only).
  * Memory allocation errors may also be returned.
  */
 SLICOT_EXPORT
 int slicot_sg03bd(char dico, char fact, char trans, int n, int m,
                   double* a, int lda, double* e, int lde,
                   double* q, int ldq, double* z, int ldz,
                   double* b, int ldb, double* scale,
                   double* alphar, double* alphai, double* beta,
                   int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* SG03BD_H */
 