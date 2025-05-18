/**
 * @file tb03ad.h
 * @brief Header for C wrapper of SLICOT routine TB03AD.
 */

#ifndef SLICOT_WRAPPER_TB03AD_H
#define SLICOT_WRAPPER_TB03AD_H

#include "slicot_utils.h" // Provides SLICOT_EXPORT macro

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Finds a relatively prime left or right polynomial matrix representation
 * and a minimal state-space realization for a given system.
 * @details This is a C wrapper for the SLICOT Fortran routine TB03AD.
 * Matrices A, B, C are modified in place to return Amin, Bmin, Cmin.
 * NR, INDEX, PCOEFF, QCOEFF, VCOEFF are outputs.
 * Workspace is allocated internally.
 * For PCOEFF, QCOEFF, VCOEFF, the C arrays are expected to be flat arrays
 * representing 3D data in row-major order: array[dim1_idx][dim2_idx][dim3_idx]
 * where dim3_idx (coefficient index) is the fastest varying.
 *
 * @param leri_param    (input) Character, 'L' for left, 'R' for right polynomial matrix representation.
 * @param equil_param   (input) Character, 'S' to scale, 'N' not to scale (A,B,C) before minimal realization.
 * @param n_param       (input) Order of the original state matrix A. N >= 0.
 * @param m_param       (input) Number of system inputs. M >= 0.
 * @param p_param       (input) Number of system outputs. P >= 0.
 * @param a_io          (input/output) On entry, the N x N state matrix A.
 * On exit, the leading NR x NR part contains the minimal state matrix Amin.
 * @param lda           (input) Leading dimension of A.
 * If row_major=0, lda >= max(1,N).
 * If row_major=1, lda >= max(1,N) (number of columns).
 * @param b_io          (input/output) On entry, the N x M input matrix B.
 * (Fortran array B is dimensioned (LDB, MAX(M,P)) for workspace).
 * On exit, the leading NR x M part contains Bmin.
 * @param ldb           (input) Leading dimension of B.
 * If row_major=0, ldb >= max(1,N).
 * If row_major=1, ldb >= max(1,M) (actual cols of B data); the wrapper handles Fortran's MAX(M,P) col requirement.
 * @param c_io          (input/output) On entry, the P x N output matrix C.
 * (Fortran array C is dimensioned (LDC, N) but LDC must be >= MAX(1,M,P) for workspace).
 * On exit, the leading P x NR part contains Cmin.
 * @param ldc           (input) Leading dimension of C.
 * If row_major=0, ldc >= max(1,P) (actual rows of C data); wrapper ensures Fortran LDC >= MAX(1,M,P).
 * If row_major=1, ldc >= max(1,N) (number of columns).
 * @param d_in          (input) The P x M direct transmission matrix D.
 * (Fortran array D is dimensioned (LDD, MAX(M,P)) for workspace).
 * @param ldd           (input) Leading dimension of D.
 * If row_major=0, ldd >= max(1,P).
 * If row_major=1, ldd >= max(1,M) (actual cols of D data).
 * @param[out] nr_out   (output) Order of the minimal state-space representation (Amin,Bmin,Cmin).
 * @param[out] index_out (output) Integer array. Dimension P if LERI='L', M if LERI='R'.
 * Contains maximum degrees of polynomials in rows/columns of P(s).
 * @param[out] pcoeff_out (output) Double array for coefficients of P(s).
 * Dimensions (porm x porm x (N+1)). porm = P if LERI='L', M if LERI='R'.
 * Stored flat, row-major: PCOEFF[i][j][k] with k fastest.
 * @param ldpco1_c      (input) First dimension of the C array pcoeff_out (porm).
 * @param ldpco2_c      (input) Second dimension of the C array pcoeff_out (porm).
 * @param[out] qcoeff_out (output) Double array for coefficients of Q(s).
 * Dimensions (porm x porp x (N+1)) if LERI='L', (porp x porm x (N+1)) if LERI='R'.
 * porm = P if LERI='L', M if LERI='R'. porp = M if LERI='L', P if LERI='R'.
 * Stored flat, row-major.
 * @param ldqco1_c      (input) First dimension of the C array qcoeff_out.
 * @param ldqco2_c      (input) Second dimension of the C array qcoeff_out.
 * @param[out] vcoeff_out (output) Double array for coefficients of V(s).
 * Dimensions (porm x NR_out x (N+1)). porm = P if LERI='L', M if LERI='R'.
 * Stored flat, row-major.
 * @param ldvco1_c      (input) First dimension of the C array vcoeff_out (porm).
 * @param ldvco2_c      (input) Second dimension of the C array vcoeff_out (NR_out at exit, pass N for allocation).
 * @param tol_param     (input) Tolerance for rank determination. If <= 0, default is used.
 * @param row_major     (input) Integer, 0 for column-major, 1 for row-major storage of matrices A,B,C,D and output PCOEFF,QCOEFF,VCOEFF.
 *
 * @return info Error indicator:
 * = 0: successful exit
 * < 0: if info = -i, the i-th argument had an illegal value
 * = 1: singular matrix in V(s) computation
 * = 2: singular matrix in P(s) computation
 * = SLICOT_MEMORY_ERROR (-1010): internal memory allocation failed.
 */
SLICOT_EXPORT
int slicot_tb03ad(
    char leri_param, char equil_param,
    int n_param, int m_param, int p_param,
    double* a_io, int lda,
    double* b_io, int ldb,
    double* c_io, int ldc,
    const double* d_in, int ldd,
    int* nr_out,
    int* index_out,
    double* pcoeff_out, int ldpco1_c, int ldpco2_c,
    double* qcoeff_out, int ldqco1_c, int ldqco2_c,
    double* vcoeff_out, int ldvco1_c, int ldvco2_c,
    double tol_param,
    int row_major);

#ifdef __cplusplus
}
#endif

#endif /* SLICOT_WRAPPER_TB03AD_H */
