/**
 * @file sb02md.h
 * @brief C wrapper for SLICOT routine SB02MD
 *
 * This file provides a C interface to the SLICOT routine SB02MD,
 * which solves continuous- or discrete-time algebraic Riccati
 * equations using the Schur vectors method.
 */

 #ifndef SB02MD_H
 #define SB02MD_H
 
 #include <stddef.h> // For size_t
 
 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Solves continuous or discrete algebraic Riccati equations (Schur method).
  *
  * Solves Q + A'X + XA - XGX = 0 (continuous) or
  * X = A'XA - A'XB*(R + B'XB)^-1*B'XA + Q (discrete, where G = B*inv(R)*B').
  * The routine actually uses G directly. Use SB02MT to compute G from B and R.
  *
  * @param[in] dico      Specifies the type of Riccati equation:
  * = 'C': Continuous-time.
  * = 'D': Discrete-time.
  * @param[in] hinv      If dico='D', specifies symplectic matrix construction:
  * = 'D': Construct H.
  * = 'I': Construct inv(H).
  * Not used if dico='C'.
  * @param[in] uplo      Specifies which triangle of G and Q is stored:
  * = 'U': Upper triangle.
  * = 'L': Lower triangle.
  * @param[in] scal      Specifies scaling strategy:
  * = 'G': General scaling.
  * = 'N': No scaling.
  * @param[in] sort      Specifies eigenvalue ordering in Schur form:
  * = 'S': Stable eigenvalues first.
  * = 'U': Unstable eigenvalues first.
  * @param[in] n         The order of the matrices A, G, Q, X, n >= 0.
  * @param[in,out] a     Double array, dimension (lda, n) or (n, lda).
  * On entry, the matrix A.
  * On exit, if dico='D' and info=0 or >1, contains inv(A). Otherwise unchanged.
  * @param[in] lda       The leading dimension of array A. >= max(1,n).
  * @param[in] g         Double array, dimension (ldg, n) or (n, ldg).
  * Symmetric matrix G = B*inv(R)*B'. Only triangle specified by uplo is used.
  * @param[in] ldg       The leading dimension of array G. >= max(1,n).
  * @param[in,out] q     Double array, dimension (ldq, n) or (n, ldq).
  * On entry, symmetric matrix Q. Only triangle specified by uplo is used.
  * On exit, the solution matrix X.
  * @param[in] ldq       The leading dimension of array Q. >= max(1,n).
  * @param[out] rcond    Estimated reciprocal condition number of the system solved for X.
  * @param[out] wr       Double array, dimension (2*n). Real parts of eigenvalues of H.
  * @param[out] wi       Double array, dimension (2*n). Imaginary parts of eigenvalues of H.
  * The first n elements correspond to the closed-loop spectrum.
  * @param[out] s        Double array, dimension (lds, 2*n) or (2*n, lds).
  * The ordered real Schur form S of the Hamiltonian/symplectic matrix H.
  * @param[in] lds       The leading dimension of array S. >= max(1,2*n).
  * @param[out] u        Double array, dimension (ldu, 2*n) or (2*n, ldu).
  * The orthogonal transformation matrix U such that U'*H*U = S.
  * @param[in] ldu       The leading dimension of array U. >= max(1,2*n).
  * @param[in] row_major Integer flag:
  * = 0: Arrays a, g, q, s, u are column-major (Fortran style).
  * = 1: Arrays a, g, q, s, u are row-major (C style).
  * (Arrays wr, wi are 1D).
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * = 1: A is singular (discrete-time only).
  * = 2: QR/QZ algorithm failed to compute Schur form.
  * = 3: Reordering of eigenvalues failed.
  * = 4: Less than N stable eigenvalues found after reordering.
  * = 5: System for X is singular to working precision.
  * Memory allocation errors may also be returned.
  */
 int slicot_sb02md(char dico, char hinv, char uplo, char scal, char sort,
                   int n, double* a, int lda, const double* g, int ldg,
                   double* q, int ldq, double* rcond,
                   double* wr, double* wi,
                   double* s, int lds, double* u, int ldu,
                   int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* SB02MD_H */
 