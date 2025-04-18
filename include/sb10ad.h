/**
 * @file sb10ad.h
 * @brief C wrapper for SLICOT routine SB10AD
 *
 * This file provides a C interface to the SLICOT routine SB10AD,
 * which computes an H-infinity optimal controller for a continuous-time
 * system using modified Glover's and Doyle's formulas.
 */

 #ifndef SB10AD_H
 #define SB10AD_H
 
 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Computes H-infinity optimal controller and closed-loop system.
  *
  * Computes the matrices of an H-infinity optimal n-state controller
  * K = [AK, BK; CK, DK] for the system P = [A, B; C, D] using modified
  * Glover's and Doyle's formulas. It estimates the minimal possible gamma
  * and computes the controller for that gamma. It also computes the
  * matrices of the resulting closed-loop system G = [AC, BC; CC, DC].
  * Assumes standard H-infinity problem setup and conditions (A1)-(A4).
  *
  * @param[in] job       Strategy for reducing gamma: 1=Bisection, 2=Scan, 3=Bisect+Scan, 4=Suboptimal only.
  * @param[in] n         Order of the system (state dimension), n >= 0.
  * @param[in] m         Column size of B and D, m >= 0.
  * @param[in] np        Row size of C and D, np >= 0.
  * @param[in] ncon      Number of control inputs (columns of B2, rows of CK, DK), M >= ncon >= 0.
  * @param[in] nmeas     Number of measurements (rows of C2, columns of BK, DK), NP >= nmeas >= 0.
  * @param[in,out] gamma On entry, initial large gamma. On exit, minimal estimated gamma. gamma >= 0.
  * @param[in] a         Double array, dimension (lda, n) or (n, lda). System matrix A.
  * @param[in] lda       Leading dimension of A. >= max(1,n).
  * @param[in] b         Double array, dimension (ldb, m) or (n, ldb). System matrix B = [B1 B2].
  * @param[in] ldb       Leading dimension of B. >= max(1,n).
  * @param[in] c         Double array, dimension (ldc, n) or (np, ldc). System matrix C = [C1; C2].
  * @param[in] ldc       Leading dimension of C. >= max(1,np).
  * @param[in] d         Double array, dimension (ldd, m) or (np, ldd). System matrix D = [D11 D12; D21 D22].
  * @param[in] ldd       Leading dimension of D. >= max(1,np).
  * @param[out] ak       Double array, dimension (ldak, n) or (n, ldak). Controller state matrix AK.
  * @param[in] ldak      Leading dimension of AK. >= max(1,n).
  * @param[out] bk       Double array, dimension (ldbk, nmeas) or (n, ldbk). Controller input matrix BK.
  * @param[in] ldbk      Leading dimension of BK. >= max(1,n).
  * @param[out] ck       Double array, dimension (ldck, n) or (ncon, ldck). Controller output matrix CK.
  * @param[in] ldck      Leading dimension of CK. >= max(1,ncon).
  * @param[out] dk       Double array, dimension (lddk, nmeas) or (ncon, lddk). Controller matrix DK.
  * @param[in] lddk      Leading dimension of DK. >= max(1,ncon).
  * @param[out] ac       Double array, dimension (ldac, 2*n) or (2*n, ldac). Closed-loop state matrix AC.
  * @param[in] ldac      Leading dimension of AC. >= max(1,2*n).
  * @param[out] bc       Double array, dimension (ldbc, m-ncon) or (2*n, ldbc). Closed-loop input matrix BC.
  * @param[in] ldbc      Leading dimension of BC. >= max(1,2*n).
  * @param[out] cc       Double array, dimension (ldcc, 2*n) or (np-nmeas, ldcc). Closed-loop output matrix CC.
  * @param[in] ldcc      Leading dimension of CC. >= max(1,np-nmeas).
  * @param[out] dc       Double array, dimension (lddc, m-ncon) or (np-nmeas, lddc). Closed-loop matrix DC.
  * @param[in] lddc      Leading dimension of DC. >= max(1,np-nmeas).
  * @param[out] rcond    Double array, dimension (4). Reciprocal condition numbers for transformations and Riccati equations.
  * @param[in] gtol      Tolerance for gamma accuracy. If <= 0, sqrt(eps) is used.
  * @param[in] actol     Tolerance for closed-loop stability check (poles real part < actol). <= 0 for stability.
  * @param[in] row_major Integer flag:
  * = 0: All 2D arrays are column-major (Fortran style).
  * = 1: All 2D arrays are row-major (C style).
  * (Array rcond is 1D).
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * = 1..5: Rank condition failures (A3, A4, D12, D21) or SVD failure.
  * = 6: Controller not admissible (gamma too small initially).
  * = 7: X-Riccati solver failed.
  * = 8: Y-Riccati solver failed.
  * = 9..11: Numerical problems during controller construction/check.
  * = 12: Stabilizing controller could not be found.
  * Memory allocation errors may also be returned.
  */
 SLICOT_C_WRAPPER_API
 int slicot_sb10ad(int job, int n, int m, int np, int ncon, int nmeas,
                   double* gamma, const double* a, int lda,
                   const double* b, int ldb, const double* c, int ldc,
                   const double* d, int ldd, double* ak, int ldak,
                   double* bk, int ldbk, double* ck, int ldck,
                   double* dk, int lddk, double* ac, int ldac,
                   double* bc, int ldbc, double* cc, int ldcc,
                   double* dc, int lddc, double* rcond, double gtol,
                   double actol, int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* SB10AD_H */
 