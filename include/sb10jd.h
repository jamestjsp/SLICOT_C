/**
 * @file sb10jd.h
 * @brief Header for C wrapper of SLICOT routine SB10JD.
 */

#ifndef SLICOT_WRAPPER_SB10JD_H
#define SLICOT_WRAPPER_SB10JD_H

#include "slicot_utils.h" // Provides SLICOT_EXPORT macro

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Converts a descriptor state-space system into regular state-space form.
 * @details This is a C wrapper for the SLICOT Fortran routine SB10JD.
 * Matrices A, B, C, D are modified to contain the matrices Ad, Bd, Cd, Dd
 * of the converted system. Matrix E is used as input and its content is
 * destroyed. Workspace is allocated internally.
 *
 * @param n_param (input) The order of the descriptor system. N >= 0.
 * @param m_param (input) The column size of the matrix B. M >= 0.
 * @param np_param (input) The row size of the matrix C. NP >= 0.
 * @param a      (input/output) On entry, state matrix A (N x N).
 * On exit, state matrix Ad (NSYS x NSYS) of the converted system.
 * @param lda    Leading dimension of A.
 * @param b      (input/output) On entry, input matrix B (N x M).
 * On exit, input matrix Bd (NSYS x M) of the converted system.
 * @param ldb    Leading dimension of B.
 * @param c      (input/output) On entry, output matrix C (NP x N).
 * On exit, output matrix Cd (NP x NSYS) of the converted system.
 * @param ldc    Leading dimension of C.
 * @param d      (input/output) On entry, matrix D (NP x M).
 * On exit, matrix Dd (NP x M) of the converted system.
 * @param ldd    Leading dimension of D.
 * @param e      (input/output) On entry, matrix E (N x N) of the descriptor system.
 * On exit, E contains no useful information.
 * @param lde    Leading dimension of E.
 * @param[out] nsys The order of the converted state-space system.
 * @param row_major Specifies matrix storage: 0 for column-major, 1 for row-major.
 *
 * @return info Error indicator:
 * = 0:  successful exit
 * < 0:  if info = -i, the i-th argument had an illegal value
 * = 1:  SVD algorithm did not converge.
 * = SLICOT_MEMORY_ERROR (-1010): internal memory allocation failed.
 */
SLICOT_EXPORT
int slicot_sb10jd(int n_param, int m_param, int np_param,
                  double* a, int lda,
                  double* b, int ldb,
                  double* c, int ldc,
                  double* d, int ldd,
                  double* e, int lde,
                  int* nsys,
                  int row_major);

#ifdef __cplusplus
}
#endif

#endif /* SLICOT_WRAPPER_SB10JD_H */
