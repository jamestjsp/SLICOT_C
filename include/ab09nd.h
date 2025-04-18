/**
 * @file ab09nd.h
 * @brief C wrapper for SLICOT routine AB09ND
 *
 * This file provides a C interface to the SLICOT routine AB09ND,
 * which computes a reduced order model (Ar,Br,Cr,Dr) for the ALPHA-stable
 * part of an original state-space representation (A,B,C,D) using
 * Singular Perturbation Approximation (SPA) methods.
 */

 #ifndef AB09ND_H
 #define AB09ND_H
 
 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Singular Perturbation Approximation model reduction for the ALPHA-stable part.
  *
  * Computes a reduced order model (Ar,Br,Cr,Dr) for an original
  * state-space representation (A,B,C,D) by using either the
  * square-root or the balancing-free square-root Singular
  * Perturbation Approximation (SPA) model reduction method for the
  * ALPHA-stable part of the system. The original A, B, C, D arrays are
  * overwritten with the reduced system.
  *
  * @param[in] dico      Specifies the type of the original system:
  * = 'C': Continuous-time system.
  * = 'D': Discrete-time system.
  * @param[in] job       Specifies the model reduction approach:
  * = 'B': Use square-root SPA method.
  * = 'N': Use balancing-free square-root SPA method.
  * @param[in] equil     Specifies whether to preliminarily equilibrate (A,B,C):
  * = 'S': Perform equilibration (scaling).
  * = 'N': Do not perform equilibration.
  * @param[in] ordsel    Specifies the order selection method:
  * = 'F': The resulting order NR is fixed by input NR.
  * = 'A': The resulting order NR is automatically determined by TOL1.
  * @param[in] n         The order of the original matrix A, n >= 0.
  * @param[in] m         The number of system inputs, m >= 0.
  * @param[in] p         The number of system outputs, p >= 0.
  * @param[in,out] nr    On entry (if ordsel='F'), the desired reduced order, 0<=nr<=n.
  * On exit, the computed order of the reduced model (includes unstable part).
  * @param[in] alpha     ALPHA-stability boundary value.
  * @param[in,out] a     Double array, dimension (lda, n) or (n, lda).
  * On entry, the state matrix A.
  * On exit, the leading nr-by-nr part contains the reduced state matrix Ar.
  * @param[in] lda       The leading dimension of array A. >= max(1,n).
  * @param[in,out] b     Double array, dimension (ldb, m) or (n, ldb).
  * On entry, the input matrix B.
  * On exit, the leading nr-by-m part contains the reduced input matrix Br.
  * @param[in] ldb       The leading dimension of array B. >= max(1,n).
  * @param[in,out] c     Double array, dimension (ldc, n) or (p, ldc).
  * On entry, the output matrix C.
  * On exit, the leading p-by-nr part contains the reduced output matrix Cr.
  * @param[in] ldc       The leading dimension of array C. >= max(1,p).
  * @param[in,out] d     Double array, dimension (ldd, m) or (p, ldd).
  * On entry, the feedthrough matrix D.
  * On exit, the leading p-by-m part contains the reduced feedthrough matrix Dr.
  * @param[in] ldd       The leading dimension of array D. >= max(1,p).
  * @param[out] ns       The dimension of the ALPHA-stable subsystem.
  * @param[out] hsv      Double array, dimension (n). The leading ns elements contain
  * the Hankel singular values of the ALPHA-stable part.
  * @param[in] tol1      Tolerance for order selection (if ordsel='A'). If tol1<=0, a default is used. Ignored if ordsel='F'.
  * @param[in] tol2      Tolerance for determining minimal realization order of stable part. If tol2<=0, a default is used. Must be <= tol1 if both > 0.
  * @param[out] iwarn    Warning indicator (0=no warning, 1=NR>minimal, 2=NR<unstable).
  * @param[in] row_major Integer flag:
  * = 0: Arrays a, b, c, d are column-major (Fortran style).
  * = 1: Arrays a, b, c, d are row-major (C style).
  * (Array hsv is 1D).
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * = 1: Schur decomposition failed.
  * = 2: Eigenvalue reordering failed.
  * = 3: Hankel singular value computation failed.
  * Memory allocation errors may also be returned.
  */
 SLICOT_C_WRAPPER_API
 int slicot_ab09nd(char dico, char job, char equil, char ordsel,
                   int n, int m, int p, int* nr, double alpha,
                   double* a, int lda,
                   double* b, int ldb,
                   double* c, int ldc,
                   double* d, int ldd,
                   int* ns, double* hsv, double tol1, double tol2,
                   int* iwarn, int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* AB09ND_H */
 