/**
 * @file sg02ad.h
 * @brief C wrapper for SLICOT routine SG02AD
 *
 * This file provides a C interface to the SLICOT routine SG02AD,
 * which solves continuous- or discrete-time algebraic Riccati equations
 * for descriptor systems using the generalized Schur vectors method.
 */

 #ifndef SG02AD_H
 #define SG02AD_H
 
 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Solves continuous or discrete generalized algebraic Riccati equations (Descriptor).
  *
  * Solves Q + A'XE + E'XA - (L+E'XB)*inv(R)*(L+E'XB)' = 0 (continuous) or
  * E'XE = A'XA - (L+A'XB)*inv(R + B'XB)*(L+A'XB)' + Q (discrete).
  * Handles factored Q=C'C, R=D'D, and zero L. Can optionally use G=B*inv(R)*B'.
  * Assumes E is nonsingular.
  *
  * @param[in] dico      Specifies the type of Riccati equation: 'C' or 'D'.
  * @param[in] jobb      Specifies if G = B*inv(R)*B' is provided: 'B' or 'G'.
  * @param[in] fact      Specifies if Q and/or R are factored: 'N', 'C', 'D', or 'B'.
  * @param[in] uplo      Specifies which triangle of symmetric matrices is stored: 'U' or 'L'.
  * @param[in] jobl      Specifies if L is zero (not used if JOBB='G'): 'Z' or 'N'.
  * @param[in] scal      Specifies scaling strategy (not used if JOBB='G'): 'G' or 'N'.
  * @param[in] sort      Specifies eigenvalue ordering: 'S' (stable first) or 'U' (unstable first).
  * @param[in] acc       Specifies if iterative refinement is used for solving for X: 'R' or 'N'.
  * @param[in] n         Order of matrices A, E, Q, X, n >= 0.
  * @param[in] m         Number of inputs (columns of B, L; order of R), m >= 0. (Not used if JOBB='G').
  * @param[in] p         Number of outputs (rows of C, D), p >= 0. (Used only if FACT involves C or D).
  * @param[in] a         Double array, dimension (lda, n) or (n, lda). Matrix A.
  * @param[in] lda       Leading dimension of A. >= max(1,n).
  * @param[in] e         Double array, dimension (lde, n) or (n, lde). Matrix E.
  * @param[in] lde       Leading dimension of E. >= max(1,n).
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
  * @param[out] rcondu   Estimated reciprocal condition number of the system solved for X.
  * @param[out] x        Double array, dimension (ldx, n) or (n, ldx). The solution matrix X.
  * @param[in] ldx       Leading dimension of X. >= max(1,n).
  * @param[out] alfar    Double array, dimension (2*n). Real parts of generalized eigenvalues.
  * @param[out] alfai    Double array, dimension (2*n). Imaginary parts of generalized eigenvalues.
  * @param[out] beta     Double array, dimension (2*n). Scaling factors for generalized eigenvalues.
  * @param[out] s        Double array, dimension (lds, *) or (2*n+m or 2*n, lds). Ordered real Schur form S.
  * @param[in] lds       Leading dimension of S. >=max(1, 2*n+m) if jobb='B', >=max(1, 2*n) if jobb='G'.
  * @param[out] t        Double array, dimension (ldt, 2*n) or (2*n+m or 2*n, ldt). Ordered upper triangular form T.
  * @param[in] ldt       Leading dimension of T. >=max(1, 2*n+m) if jobb='B', >=max(1, 2*n) if jobb='G'.
  * @param[out] u        Double array, dimension (ldu, 2*n) or (2*n, ldu). Orthogonal transformation matrix U.
  * @param[in] ldu       Leading dimension of U. >= max(1, 2*n).
  * @param[in] tol       Tolerance for rank determination. If tol<=0, default used. Not used if JOBB='G'.
  * @param[out] iwarn    Warning indicator (0=no warning, 1=potential inaccuracy).
  * @param[in] row_major Integer flag:
  * = 0: All 2D arrays are column-major.
  * = 1: All 2D arrays are row-major.
  * (Arrays alfar, alfai, beta are 1D).
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * = 1: Extended pencil is singular.
  * = 2: QZ algorithm failed.
  * = 3: Reordering of eigenvalues failed.
  * = 4: Reordered eigenvalues do not satisfy stability condition.
  * = 5: Computed dimension of solution does not equal N.
  * = 6: Spectrum too close to stability boundary.
  * = 7: Singular matrix encountered during solution for X.
  * Memory allocation errors may also be returned.
  */
 SLICOT_C_WRAPPER_API
 int slicot_sg02ad(char dico, char jobb, char fact, char uplo, char jobl, char scal, char sort, char acc,
                   int n, int m, int p,
                   const double* a, int lda, const double* e, int lde,
                   const double* b, int ldb, const double* q, int ldq,
                   const double* r, int ldr, const double* l, int ldl,
                   double* rcondu, double* x, int ldx,
                   double* alfar, double* alfai, double* beta,
                   double* s, int lds, double* t, int ldt,
                   double* u, int ldu, double tol, int* iwarn, int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* SG02AD_H */
 