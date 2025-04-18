/**
 * @file sb10fd.h
 * @brief C wrapper for SLICOT routine SB10FD
 *
 * This file provides a C interface to the SLICOT routine SB10FD,
 * which computes an H-infinity (sub)optimal state controller for a
 * continuous-time system.
 */

 #ifndef SB10FD_H
 #define SB10FD_H
 
 #include <stddef.h> // For size_t
 
 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Computes H-infinity (sub)optimal state controller for a continuous-time system.
  *
  * Computes the matrices of an H-infinity (sub)optimal n-state controller
  * K = [AK, BK; CK, DK] for the continuous-time system P = [A, B; C, D]
  * for a given value of gamma. Uses modified Glover's and Doyle's formulas.
  * Assumes standard H-infinity problem setup and conditions (A1)-(A4).
  *
  * @param[in] n         Order of the system (state dimension), n >= 0.
  * @param[in] m         Column size of B and D, m >= 0.
  * @param[in] np        Row size of C and D, np >= 0.
  * @param[in] ncon      Number of control inputs (columns of B2, rows of CK, DK), M >= ncon >= 0.
  * @param[in] nmeas     Number of measurements (rows of C2, columns of BK, DK), NP >= nmeas >= 0.
  * @param[in] gamma     The H-infinity performance level, gamma >= 0.
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
  * @param[out] rcond    Double array, dimension (4). Reciprocal condition number estimates.
  * @param[in] tol       Tolerance for rank determination. If <= 0, default used.
  * @param[in] row_major Integer flag:
  * = 0: All 2D arrays are column-major (Fortran style).
  * = 1: All 2D arrays are row-major (C style).
  * (Array rcond is 1D).
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * = 1..5: Rank condition failures (A3, A4, D12, D21) or SVD failure.
  * = 6: Controller not admissible (gamma too small).
  * = 7: X-Riccati solver failed.
  * = 8: Y-Riccati solver failed.
  * = 9: Determinant involving Tu, Ty, D11HAT is zero.
  * Memory allocation errors may also be returned.
  */
 int slicot_sb10fd(int n, int m, int np, int ncon, int nmeas,
                   double gamma, const double* a, int lda,
                   const double* b, int ldb, const double* c, int ldc,
                   const double* d, int ldd, double* ak, int ldak,
                   double* bk, int ldbk, double* ck, int ldck,
                   double* dk, int lddk, double* rcond, double tol,
                   int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* SB10FD_H */
 