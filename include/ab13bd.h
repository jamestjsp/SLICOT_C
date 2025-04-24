/**
 * @file ab13bd.h
 * @brief C wrapper for SLICOT function AB13BD
 *
 * This file provides a C interface to the SLICOT function AB13BD,
 * which computes the H2 or L2 norm of a system (A,B,C,D).
 */

 #ifndef AB13BD_H
 #define AB13BD_H
 
 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Computes the H2 or L2 norm of a system.
  *
  * Computes the H2 or L2 norm of the transfer-function matrix G
  * of the system (A,B,C,D). G must not have poles on the stability
  * boundary (imaginary axis or unit circle). For H2 norm, the system
  * must be stable. The matrices A, B, C, D are overwritten by the
  * state-space representation of the numerator factor Q of a right
  * coprime factorization G = Q*R^-1 (Q=G if stable).
  *
  * @param[in] dico      Specifies the type of the system:
  * = 'C': Continuous-time system.
  * = 'D': Discrete-time system.
  * @param[in] jobn      Specifies the norm to compute:
  * = 'H': H2-norm (system must be stable).
  * = 'L': L2-norm.
  * @param[in] n         The order of the matrix A, n >= 0.
  * @param[in] m         The number of inputs (columns of B, D), m >= 0.
  * @param[in] p         The number of outputs (rows of C, D), p >= 0.
  * @param[in,out] a     Double array, dimension (lda, n) or (n, lda).
  * On entry, the state matrix A. On exit, the state matrix AQ of the numerator Q.
  * @param[in] lda       The leading dimension of array A. >= max(1,n).
  * @param[in,out] b     Double array, dimension (ldb, m) or (n, ldb).
  * On entry, the input matrix B. On exit, the input matrix BQ of the numerator Q.
  * @param[in] ldb       The leading dimension of array B. >= max(1,n).
  * @param[in,out] c     Double array, dimension (ldc, n) or (p, ldc).
  * On entry, the output matrix C. On exit, the output matrix CQ of the numerator Q.
  * @param[in] ldc       The leading dimension of array C. >= max(1,p).
  * @param[in,out] d     Double array, dimension (ldd, m) or (p, ldd).
  * On entry, the feedthrough matrix D (must be zero if dico='C').
  * On exit, the feedthrough matrix DQ of the numerator Q.
  * @param[in] ldd       The leading dimension of array D. >= max(1,p).
  * @param[out] nq       The order of the resulting numerator state matrix AQ.
  * @param[in] tol       Tolerance used for controllability tests. If tol<=0, a default is used.
  * @param[out] iwarn    Warning indicator (related to eigenvalue assignment).
  * @param[out] info     Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * = 1: Schur reduction failed.
  * = 2: Eigenvalue reordering failed.
  * = 3: System has pole on stability boundary.
  * = 4: Lyapunov equation solution failed (singular).
  * = 5: D is non-zero for continuous-time H2 norm.
  * = 6: System unstable for H2 norm computation.
  * Memory allocation errors may also be returned by the wrapper.
  * @param[in] row_major Integer flag:
  * = 0: Arrays a, b, c, d are column-major (Fortran style).
  * = 1: Arrays a, b, c, d are row-major (C style).
  *
  * @return              The computed H2 or L2 norm if info = 0, otherwise typically 0.0.
  * Check info for success status.
  */
 SLICOT_EXPORT
 double slicot_ab13bd(char dico, char jobn, int n, int m, int p,
                      double* a, int lda, double* b, int ldb,
                      double* c, int ldc, double* d, int ldd,
                      int* nq, double tol, int* iwarn, int* info,
                      int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* AB13BD_H */
 