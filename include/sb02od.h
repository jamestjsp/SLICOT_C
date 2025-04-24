/**
 * @file sb02od.h
 * @brief C wrapper for SLICOT routine SB02OD
 *
 * This file provides a C interface to the SLICOT routine SB02OD,
 * which solves continuous- or discrete-time algebraic Riccati equations
 * using the generalized Schur vectors method.
 */

 #ifndef SB02OD_H
 #define SB02OD_H
 
 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Solves continuous or discrete algebraic Riccati equations (Generalized Schur method).
  *
  * Solves Q + A'X + XA - (L+XB)*inv(R)*(L+XB)' = 0 (continuous) or
  * X = A'XA - (L+A'XB)*inv(R + B'XB)*(L+A'XB)' + Q (discrete).
  * Handles factored Q=C'C, R=D'D, and zero L. Can optionally use G=B*inv(R)*B'.
  *
  * @param[in] dico      Specifies the type of Riccati equation:
  * = 'C': Continuous-time.
  * = 'D': Discrete-time.
  * @param[in] jobb      Specifies if G = B*inv(R)*B' is provided instead of B and R:
  * = 'B': B and R (or factors D) are provided.
  * = 'G': G is provided in array B. R and L are not used. JOBL='Z' assumed.
  * @param[in] fact      Specifies if Q and/or R are factored:
  * = 'N': Q and R are provided.
  * = 'C': Factor C (for Q=C'C) and matrix R are provided.
  * = 'D': Matrix Q and factor D (for R=D'D) are provided.
  * = 'B': Factors C and D are provided.
  * @param[in] uplo      Specifies which triangle of Q, R, G is stored:
  * = 'U': Upper triangle.
  * = 'L': Lower triangle.
  * Not used if FACT='C' or 'B' for Q, or FACT='D' or 'B' for R.
  * @param[in] jobl      Specifies if L is zero (not used if JOBB='G'):
  * = 'Z': L is zero.
  * = 'N': L is nonzero and provided.
  * @param[in] sort      Specifies eigenvalue ordering in generalized Schur form:
  * = 'S': Stable eigenvalues first.
  * = 'U': Unstable eigenvalues first.
  * @param[in] n         The order of matrix A, n >= 0.
  * @param[in] m         Number of inputs (columns of B, L; order of R), m >= 0. (Not used if JOBB='G').
  * @param[in] p         Number of outputs (rows of C, D), p >= 0. (Used only if FACT involves C or D).
  * @param[in] a         Double array, dimension (lda, n) or (n, lda). The matrix A.
  * @param[in] lda       Leading dimension of A. >= max(1,n).
  * @param[in] b         Double array. If JOBB='B', dimension (ldb, m) or (n, ldb), contains B.
  * If JOBB='G', dimension (ldb, n) or (n, ldb), contains symmetric G (triangle specified by uplo).
  * @param[in] ldb       Leading dimension of B/G array. >= max(1,n).
  * @param[in] q         Double array. If FACT='N' or 'D', dimension (ldq, n) or (n, ldq), contains symmetric Q.
  * If FACT='C' or 'B', dimension (ldq, n) or (p, ldq), contains C.
  * @param[in] ldq       Leading dimension of Q/C array. >=max(1,n) if FACT='N'/'D', >=max(1,p) if FACT='C'/'B'.
  * @param[in] r         Double array. If JOBB='B' and FACT='N' or 'C', dimension (ldr, m) or (m, ldr), contains symmetric R.
  * If JOBB='B' and FACT='D' or 'B', dimension (ldr, m) or (p, ldr), contains D.
  * Not referenced if JOBB='G'.
  * @param[in] ldr       Leading dimension of R/D array. >=max(1,m) if FACT='N'/'C', >=max(1,p) if FACT='D'/'B'. Ignored if JOBB='G'.
  * @param[in] l         Double array, dimension (ldl, m) or (n, ldl). Matrix L. Not referenced if JOBL='Z' or JOBB='G'.
  * @param[in] ldl       Leading dimension of L. >=max(1,n) if JOBL='N' and JOBB='B', else >=1.
  * @param[out] rcond    Estimated reciprocal condition number of the system solved for X.
  * @param[out] x        Double array, dimension (ldx, n) or (n, ldx). The solution matrix X.
  * @param[in] ldx       Leading dimension of X. >= max(1,n).
  * @param[out] alfar    Double array, dimension (2*n). Real parts of generalized eigenvalues.
  * @param[out] alfai    Double array, dimension (2*n). Imaginary parts of generalized eigenvalues.
  * @param[out] beta     Double array, dimension (2*n). Scaling factors for generalized eigenvalues.
  * The first n eigenvalues correspond to the closed-loop spectrum.
  * @param[out] s        Double array, dimension (lds, *) or (2*n+m or 2*n, lds).
  * Ordered real Schur form S of the pencil matrix (or Hamiltonian if dico='C', jobb='G').
  * @param[in] lds       Leading dimension of S. >=max(1, 2*n+m) if jobb='B', >=max(1, 2*n) if jobb='G'.
  * @param[out] t        Double array, dimension (ldt, 2*n) or (2*n+m or 2*n, ldt).
  * Ordered upper triangular form T of the second matrix in the pencil.
  * Not referenced if dico='C' and jobb='G'.
  * @param[in] ldt       Leading dimension of T. >=max(1, 2*n+m) if jobb='B', >=max(1, 2*n) if jobb='G'/'D', >=1 if jobb='G'/'C'.
  * @param[out] u        Double array, dimension (ldu, 2*n) or (2*n, ldu).
  * Orthogonal transformation matrix U.
  * @param[in] ldu       Leading dimension of U. >= max(1, 2*n).
  * @param[in] tol       Tolerance for rank determination during pencil reduction. If tol<=0, default used.
  * @param[in] row_major Integer flag:
  * = 0: Arrays a, b/g, q/c, r/d, l, x, s, t, u are column-major.
  * = 1: Arrays a, b/g, q/c, r/d, l, x, s, t, u are row-major.
  * (Arrays alfar, alfai, beta are 1D).
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * = 1: Extended pencil is singular.
  * = 2: QZ/QR algorithm failed.
  * = 3: Reordering of eigenvalues failed.
  * = 4: Reordered eigenvalues do not satisfy stability condition.
  * = 5: Computed dimension of solution does not equal N.
  * = 6: Singular matrix encountered during solution for X.
  * Memory allocation errors may also be returned.
  */
 SLICOT_EXPORT
 int slicot_sb02od(char dico, char jobb, char fact, char uplo, char jobl, char sort,
                   int n, int m, int p,
                   const double* a, int lda, const double* b, int ldb,
                   const double* q, int ldq, const double* r, int ldr,
                   const double* l, int ldl, double* rcond,
                   double* x, int ldx, double* alfar, double* alfai, double* beta,
                   double* s, int lds, double* t, int ldt,
                   double* u, int ldu, double tol, int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* SB02OD_H */
 