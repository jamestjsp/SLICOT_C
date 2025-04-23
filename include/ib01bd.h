/**
 * @file ib01bd.h
 * @brief Header for C wrapper of SLICOT routine IB01BD.
 * @details This header provides the interface for the C wrapper of the
 * SLICOT routine IB01BD, which estimates system matrices (A, B, C, D),
 * noise covariances (Q, Ry, S), and the Kalman gain (K) using processed
 * data from IB01AD.
 * **Workspace (IWORK, DWORK, BWORK) is allocated internally.**
 */

 #ifndef SLICOT_WRAPPER_IB01BD_H
 #define SLICOT_WRAPPER_IB01BD_H
 
 #include <stdbool.h>      // Include for bool type (though not directly used in signature)
 #include "slicot_utils.h" // Provides SLICOT_C_WRAPPER_API macro
 
 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Estimates system matrices, covariances, and Kalman gain from processed Hankel data.
  * @details This routine uses the processed upper triangular factor R, computed by IB01AD,
  * to estimate the system matrices (A, C, B, D), noise covariance matrices (Q, Ry, S),
  * and optionally the Kalman gain matrix (K) of a linear time-invariant state-space model,
  * using MOESP or N4SID subspace identification methods.
  * Workspace is allocated internally by this wrapper.
  *
  * @param meth Specifies the subspace identification method:
  * 'M' = MOESP algorithm.
  * 'N' = N4SID algorithm.
  * 'C' = Combined (MOESP for A/C, N4SID for B/D).
  * @param job Specifies which matrices should be computed:
  * 'A' = Compute A, B, C, D.
  * 'C' = Compute A, C only.
  * 'B' = Compute B only (requires A, C as input if METH='N' or 'C').
  * 'D' = Compute B, D only (requires A, C as input if METH='N' or 'C').
  * @param jobck Specifies computation of covariances and Kalman gain:
  * 'C' = Compute covariance matrices (Q, Ry, S) only.
  * 'K' = Compute covariance matrices and Kalman gain (K).
  * 'N' = Do not compute covariances or Kalman gain.
  * @param nobr (Input) Number of block rows (s) used in IB01AD. nobr > 1.
  * @param n (Input) The estimated system order. nobr > n > 0.
  * @param m (Input) Number of system inputs. m >= 0.
  * @param l (Input) Number of system outputs. l > 0.
  * @param nsmpl (Input) Total number of samples used for covariance calculation (if jobck='C' or 'K').
  * nsmpl >= 2*(m+l)*nobr. Not used if jobck='N'.
  * @param r (Input) Processed upper triangular factor from IB01AD.
  * Dimensions (ldr, 2*(m+l)*nobr). **Must be column-major.**
  * @param ldr (Input) Leading dimension of R. ldr >= 2*(m+l)*nobr.
  * @param a [in,out] System state matrix A (n x n).
  * Input if meth='N'/'C' and job='B'/'D'. Output if job='A'/'C'.
  * Stored column-wise if row_major=0, row-wise if row_major=1.
  * @param lda (Input) Leading dimension of C array A.
  * If row_major=0, lda >= max(1, n).
  * If row_major=1, lda >= max(1, n) (number of columns).
  * Required if A is input or output.
  * @param c [in,out] System output matrix C (l x n).
  * Input if meth='N'/'C' and job='B'/'D'. Output if job='A'/'C'.
  * Stored column-wise if row_major=0, row-wise if row_major=1.
  * @param ldc (Input) Leading dimension of C array C.
  * If row_major=0, ldc >= max(1, l).
  * If row_major=1, ldc >= max(1, n) (number of columns).
  * Required if C is input or output.
  * @param b [out] System input matrix B (n x m). Output if job='A'/'B'/'D' and m > 0.
  * Stored column-wise if row_major=0, row-wise if row_major=1.
  * @param ldb (Input) Leading dimension of C array B.
  * If row_major=0, ldb >= max(1, n).
  * If row_major=1, ldb >= max(1, m) (number of columns).
  * Required if B is output.
  * @param d [out] System input-output matrix D (l x m). Output if job='A'/'D' and m > 0.
  * Stored column-wise if row_major=0, row-wise if row_major=1.
  * @param ldd (Input) Leading dimension of C array D.
  * If row_major=0, ldd >= max(1, l).
  * If row_major=1, ldd >= max(1, m) (number of columns).
  * Required if D is output.
  * @param q [out] State covariance matrix Q (n x n). Output if jobck='C'/'K'.
  * Stored column-wise if row_major=0, row-wise if row_major=1.
  * @param ldq (Input) Leading dimension of C array Q.
  * If row_major=0, ldq >= max(1, n).
  * If row_major=1, ldq >= max(1, n) (number of columns).
  * Required if Q is output.
  * @param ry [out] Output covariance matrix Ry (l x l). Output if jobck='C'/'K'.
  * Stored column-wise if row_major=0, row-wise if row_major=1.
  * @param ldry (Input) Leading dimension of C array Ry.
  * If row_major=0, ldry >= max(1, l).
  * If row_major=1, ldry >= max(1, l) (number of columns).
  * Required if Ry is output.
  * @param s [out] State-output cross-covariance matrix S (n x l). Output if jobck='C'/'K'.
  * Stored column-wise if row_major=0, row-wise if row_major=1.
  * @param lds (Input) Leading dimension of C array S.
  * If row_major=0, lds >= max(1, n).
  * If row_major=1, lds >= max(1, l) (number of columns).
  * Required if S is output.
  * @param k [out] Kalman gain matrix K (n x l). Output if jobck='K'.
  * Stored column-wise if row_major=0, row-wise if row_major=1.
  * @param ldk (Input) Leading dimension of C array K.
  * If row_major=0, ldk >= max(1, n).
  * If row_major=1, ldk >= max(1, l) (number of columns).
  * Required if K is output.
  * @param tol (Input) Tolerance for rank estimation. If tol <= 0, a default is used.
  * @param[out] iwarn Warning indicator (see SLICOT documentation). Can be NULL if warning is not needed.
  * @param row_major Specifies matrix storage for A, C, B, D, Q, Ry, S, K:
  * 0 for column-major (Fortran default),
  * 1 for row-major (C default).
  * **Note: Input matrix R must always be column-major.**
  *
  * @return info Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value (wrapper or Fortran validation)
  * > 0: Fortran routine specific error (see SLICOT documentation for IB01BD)
  * = SLICOT_MEMORY_ERROR (-1010): internal memory allocation failed.
  */
 SLICOT_C_WRAPPER_API
 int slicot_ib01bd(char meth, char job, char jobck, int nobr, int n, int m, int l,
                   int nsmpl, double *r, int ldr,
                   double *a, int lda, double *c, int ldc,
                   double *b, int ldb, double *d, int ldd,
                   double *q, int ldq, double *ry, int ldry, double *s, int lds,
                   double *k, int ldk, double tol, int *iwarn, int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* SLICOT_WRAPPER_IB01BD_H */
 