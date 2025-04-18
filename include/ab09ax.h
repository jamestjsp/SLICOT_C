/**
 * @file ab09ax.h
 * @brief C wrapper for SLICOT routine AB09AX
 *
 * This file provides a C interface to the SLICOT routine AB09AX,
 * which computes a reduced order model (Ar,Br,Cr) for a stable system
 * (A,B,C) using Balance & Truncate methods, where A is already in
 * real Schur form. It also returns the truncation matrices T and TI.
 */

 #ifndef AB09AX_H
 #define AB09AX_H
 
 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Balance & Truncate model reduction for stable systems (Schur input).
  *
  * Computes a reduced order model (Ar,Br,Cr) for a stable original
  * state-space representation (A,B,C), where A is assumed to be in real
  * Schur form, by using either the square-root or the balancing-free
  * square-root Balance & Truncate method. The original A, B, C arrays
  * are overwritten with the reduced system. The truncation matrices T and TI
  * are also computed.
  *
  * @param[in] dico      Specifies the type of the original system:
  * = 'C': Continuous-time system.
  * = 'D': Discrete-time system.
  * @param[in] job       Specifies the model reduction approach:
  * = 'B': Use square-root Balance & Truncate method.
  * = 'N': Use balancing-free square-root B & T method.
  * @param[in] ordsel    Specifies the order selection method:
  * = 'F': The resulting order NR is fixed by input NR.
  * = 'A': The resulting order NR is automatically determined by TOL.
  * @param[in] n         The order of the original matrix A, n >= 0.
  * @param[in] m         The number of system inputs, m >= 0.
  * @param[in] p         The number of system outputs, p >= 0.
  * @param[in,out] nr    On entry (if ordsel='F'), the desired reduced order, 0<=nr<=n.
  * On exit, the computed order of the reduced model.
  * @param[in,out] a     Double array, dimension (lda, n) or (n, lda).
  * On entry, the state matrix A in real Schur form.
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
  * @param[out] t        Double array, dimension (ldt, n) or (n, ldt).
  * If nr > 0, the leading n-by-nr part contains the right truncation matrix T.
  * @param[in] ldt       The leading dimension of array T. >= max(1,n).
  * @param[out] ti       Double array, dimension (ldti, n) or (nr, ldti).
  * If nr > 0, the leading nr-by-n part contains the left truncation matrix TI.
  * @param[in] ldti      The leading dimension of array TI. >= max(1,nr) if nr>0, else >=1.
  * @param[in] tol       Tolerance for order selection (if ordsel='A'). If tol<=0, a default is used. Ignored if ordsel='F'.
  * @param[out] iwarn    Warning indicator (0 = no warning, 1 = NR was adjusted).
  * @param[in] row_major Integer flag:
  * = 0: Arrays a, b, c, t, ti are column-major (Fortran style).
  * = 1: Arrays a, b, c, t, ti are row-major (C style).
  * (Array hsv is 1D).
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * = 1: System is not stable/convergent.
  * = 2: Hankel singular value computation failed.
  * Memory allocation errors may also be returned.
  */
 SLICOT_C_WRAPPER_API
 int slicot_ab09ax(char dico, char job, char ordsel,
                   int n, int m, int p, int* nr,
                   double* a, int lda,
                   double* b, int ldb,
                   double* c, int ldc,
                   double* hsv,
                   double* t, int ldt,
                   double* ti, int ldti,
                   double tol, int* iwarn,
                   int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* AB09AX_H */
 