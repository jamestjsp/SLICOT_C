/**
 * @file ab13dd.h
 * @brief C wrapper for SLICOT routine AB13DD
 *
 * This file provides a C interface to the SLICOT routine AB13DD,
 * which computes the L-infinity norm of a standard or descriptor system.
 */

 #ifndef AB13DD_H
 #define AB13DD_H
 
 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Computes the L-infinity norm of a standard or descriptor system.
  *
  * Calculates the L-infinity norm (peak gain) of the transfer function
  * G(lambda) = C*( lambda*E - A )^-1 * B + D.
  * The norm is finite iff the pencil (A, E) has no eigenvalue on the
  * stability boundary (imaginary axis or unit circle). Assumes E is nonsingular.
  *
  * @param[in] dico      Specifies the type of the system:
  * = 'C': Continuous-time system.
  * = 'D': Discrete-time system.
  * @param[in] jobe      Specifies whether E is identity or general:
  * = 'G': E is a general square matrix.
  * = 'I': E is the identity matrix (not referenced).
  * @param[in] equil     Specifies whether to preliminarily equilibrate:
  * = 'S': Perform equilibration (scaling).
  * = 'N': Do not perform equilibration.
  * @param[in] jobd      Specifies whether D is zero or present:
  * = 'D': D is present and provided.
  * = 'Z': D is assumed zero (not referenced).
  * @param[in] n         The order of the system (order of A, E), n >= 0.
  * @param[in] m         The number of inputs (columns of B, D), m >= 0.
  * @param[in] p         The number of outputs (rows of C, D), p >= 0.
  * @param[in,out] fpeak Double array, dimension (2).
  * On entry, estimate of peak frequency [omega_real, omega_imag].
  * Set fpeak[1]=0 for infinite frequency estimate. Use {0, 1} if no estimate.
  * On exit, the frequency where the peak gain GPEAK is achieved.
  * @param[in] a         Double array for matrix A. Dimension (lda, n) or (n, lda).
  * @param[in] lda       Leading dimension of A. >= max(1,n).
  * @param[in] e         Double array for matrix E (if jobe='G'). Dimension (lde, n) or (n, lde).
  * Not referenced if jobe='I'.
  * @param[in] lde       Leading dimension of E. >= max(1,n) if jobe='G', else >=1.
  * @param[in] b         Double array for matrix B. Dimension (ldb, m) or (n, ldb).
  * @param[in] ldb       Leading dimension of B. >= max(1,n).
  * @param[in] c         Double array for matrix C. Dimension (ldc, n) or (p, ldc).
  * @param[in] ldc       Leading dimension of C. >= max(1,p).
  * @param[in] d         Double array for matrix D (if jobd='D'). Dimension (ldd, m) or (p, ldd).
  * Not referenced if jobd='Z'.
  * @param[in] ldd       Leading dimension of D. >= max(1,p) if jobd='D', else >=1.
  * @param[out] gpeak    Double array, dimension (2). The computed L-infinity norm.
  * gpeak[0] = norm value, gpeak[1] = 0.0.
  * @param[in] tol       Relative tolerance for norm accuracy, 0 <= tol < 1.
  * @param[in] row_major Integer flag:
  * = 0: Arrays a, e, b, c, d are column-major (Fortran style).
  * = 1: Arrays a, e, b, c, d are row-major (C style).
  * (Arrays fpeak, gpeak are small and unaffected).
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * = 1: E is singular (if jobe='G').
  * = 2: Eigenvalue computation failed (QR/QZ).
  * = 3: SVD computation failed.
  * = 4: Algorithm did not converge within tolerance/iterations. Increase TOL.
  * Memory allocation errors may also be returned.
  */
 int slicot_ab13dd(char dico, char jobe, char equil, char jobd,
                   int n, int m, int p, double* fpeak,
                   const double* a, int lda, const double* e, int lde,
                   const double* b, int ldb, const double* c, int ldc,
                   const double* d, int ldd, double* gpeak,
                   double tol, int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* AB13DD_H */
 