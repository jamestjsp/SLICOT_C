/**
 * @file ab09ad.h
 * @brief C wrapper for SLICOT routine AB09AD
 *
 * This file provides a C interface to the SLICOT routine AB09AD,
 * which computes a reduced order model (Ar,Br,Cr) for a stable
 * original state-space representation (A,B,C) using Balance &
 * Truncate methods.
 */

 #ifndef AB09AD_H
 #define AB09AD_H
 
 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Balance & Truncate model reduction for stable systems.
  *
  * Computes a reduced order model (Ar,Br,Cr) for a stable original
  * state-space representation (A,B,C) by using either the square-root
  * or the balancing-free square-root Balance & Truncate (B & T)
  * model reduction method. The original A, B, C arrays are overwritten.
  *
  * @param[in] dico      Specifies the type of the original system:
  * = 'C': Continuous-time system.
  * = 'D': Discrete-time system.
  * @param[in] job       Specifies the model reduction approach:
  * = 'B': Use square-root Balance & Truncate method.
  * = 'N': Use balancing-free square-root B & T method.
  * @param[in] equil     Specifies whether to preliminarily equilibrate (A,B,C):
  * = 'S': Perform equilibration (scaling).
  * = 'N': Do not perform equilibration.
  * @param[in] ordsel    Specifies the order selection method:
  * = 'F': The resulting order NR is fixed by input NR.
  * = 'A': The resulting order NR is automatically determined by TOL.
  * @param[in] n         The order of the original matrix A, n >= 0.
  * @param[in] m         The number of system inputs, m >= 0.
  * @param[in] p         The number of system outputs, p >= 0.
  * @param[in,out] nr    On entry (if ordsel='F'), the desired reduced order, 0<=nr<=n.
  * On exit, the computed order of the reduced model.
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
  * @param[out] hsv      Double array, dimension (n). Contains the Hankel singular values.
  * @param[in] tol       Tolerance for order selection (if ordsel='A'). If tol<=0, a default is used. Ignored if ordsel='F'.
  * @param[out] iwarn    Warning indicator (0 = no warning, 1 = NR was adjusted).
  * @param[in] row_major Integer flag:
  * = 0: Arrays a, b, c are column-major (Fortran style).
  * = 1: Arrays a, b, c are row-major (C style).
  * (Array hsv is 1D).
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * = 1: Schur reduction failed.
  * = 2: System is not stable/convergent.
  * = 3: Hankel singular value computation failed.
  * Memory allocation errors may also be returned.
  */
 SLICOT_C_WRAPPER_API
 int slicot_ab09ad(char dico, char job, char equil, char ordsel,
                   int n, int m, int p, int* nr,
                   double* a, int lda,
                   double* b, int ldb,
                   double* c, int ldc,
                   double* hsv, double tol, int* iwarn,
                   int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* AB09AD_H */
 