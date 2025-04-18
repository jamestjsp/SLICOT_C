/**
 * @file sb01bd.h
 * @brief C wrapper for SLICOT routine SB01BD
 *
 * This file provides a C interface to the SLICOT routine SB01BD,
 * which determines the state feedback matrix F for a given system (A,B)
 * such that the closed-loop state matrix A+B*F has specified eigenvalues.
 */

 #ifndef SB01BD_H
 #define SB01BD_H
 
 #include <stddef.h> // For size_t
 
 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Computes state feedback F for pole assignment.
  *
  * Determines the state feedback matrix F for a given system (A,B)
  * such that the closed-loop state matrix A+B*F has specified
  * eigenvalues (provided in WR, WI). The routine attempts to assign
  * MIN(NP, N-NFP) eigenvalues, where NFP is the number of eigenvalues
  * already satisfying the stability criterion defined by ALPHA.
  * The matrix A is overwritten with the Schur form of the resulting A+B*F.
  *
  * @param[in] dico      Specifies the type of the system:
  * = 'C': Continuous-time system.
  * = 'D': Discrete-time system.
  * @param[in] n         The order of the matrix A, n >= 0.
  * @param[in] m         The number of inputs (columns of B, rows of F), m >= 0.
  * @param[in] np        The number of desired eigenvalues provided in WR/WI, 0 <= np.
  * @param[in] alpha     Stability boundary value. Real parts (if dico='C') or moduli
  * (if dico='D') of eigenvalues already satisfying this are not modified.
  * @param[in,out] a     Double array, dimension (lda, n) or (n, lda).
  * On entry, the state matrix A.
  * On exit, the real Schur form of Z'*(A+B*F)*Z.
  * @param[in] lda       The leading dimension of array A. >= max(1,n).
  * @param[in] b         Double array, dimension (ldb, m) or (n, ldb). The input matrix B.
  * @param[in] ldb       The leading dimension of array B. >= max(1,n).
  * @param[in,out] wr    Double array, dimension (np). On entry, real parts of desired eigenvalues.
  * On exit, real parts of assigned eigenvalues (first NAP) and unassigned (remaining).
  * @param[in,out] wi    Double array, dimension (np). On entry, imaginary parts of desired eigenvalues.
  * On exit, imaginary parts of assigned eigenvalues (first NAP) and unassigned (remaining).
  * Complex pairs must be consecutive on entry.
  * @param[out] nfp      Number of fixed eigenvalues (already stable w.r.t. ALPHA).
  * @param[out] nap      Number of successfully assigned eigenvalues. nap = n - nfp - nup.
  * @param[out] nup      Number of uncontrollable eigenvalues detected.
  * @param[out] f        Double array, dimension (ldf, n) or (m, ldf). The computed feedback matrix F.
  * @param[in] ldf       The leading dimension of array F. >= max(1,m).
  * @param[out] z        Double array, dimension (ldz, n) or (n, ldz). The orthogonal transformation matrix Z.
  * @param[in] ldz       The leading dimension of array Z. >= max(1,n).
  * @param[in] tol       Tolerance for controllability tests. If tol<=0, a default is used.
  * @param[out] iwarn    Warning indicator (0=no warning, K=K stability violations for F).
  * @param[in] row_major Integer flag:
  * = 0: Arrays a, b, f, z are column-major (Fortran style).
  * = 1: Arrays a, b, f, z are row-major (C style).
  * (Arrays wr, wi are 1D).
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * = 1: Schur reduction failed.
  * = 2: Eigenvalue reordering failed.
  * = 3: NP < number of assignable eigenvalues; some were left unmodified.
  * = 4: Attempt to place complex pair on real eigenvalue location. NAP eigenvalues still assigned.
  * Memory allocation errors may also be returned.
  */
 int slicot_sb01bd(char dico, int n, int m, int np, double alpha,
                   double* a, int lda, const double* b, int ldb,
                   double* wr, double* wi,
                   int* nfp, int* nap, int* nup,
                   double* f, int ldf, double* z, int ldz,
                   double tol, int* iwarn, int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* SB01BD_H */
 