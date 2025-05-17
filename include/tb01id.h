/**
 * @file tb01id.h
 * @brief Header for C wrapper of SLICOT routine TB01ID.
 */

#ifndef SLICOT_WRAPPER_TB01ID_H
#define SLICOT_WRAPPER_TB01ID_H

#include "slicot_utils.h" // Provides SLICOT_EXPORT macro

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Balances a system matrix corresponding to a triplet (A,B,C).
 * @details This is a C wrapper for the SLICOT Fortran routine TB01ID.
 * Matrices A, B, C are modified in place. MAXRED is input/output.
 * SCALE is output. No explicit workspace arrays required.
 *
 * @param job_param (input) Specifies which matrices are involved:
 * 'A': All matrices (A, B, C).
 * 'B': B and A matrices.
 * 'C': C and A matrices.
 * 'N': Only A matrix (B and C not involved in balancing by JOB, but still scaled if part of S).
 * @param n_param   (input) Order of A, rows of B, columns of C. N >= 0.
 * @param m_param   (input) Columns of B. M >= 0.
 * @param p_param   (input) Rows of C. P >= 0.
 * @param maxred_io (input/output) On entry, max allowed reduction in 1-norm.
 * If <= 0.0, 10.0 is used. On exit, ratio of original to balanced norm.
 * @param a_io      (input/output) Matrix A (N x N). On exit, balanced inv(D)*A*D.
 * @param lda       Leading dimension of A.
 * @param b_io      (input/output) Matrix B (N x M). On exit, balanced inv(D)*B.
 * Not referenced if M=0 or if JOB implies B is not used.
 * @param ldb       Leading dimension of B.
 * @param c_io      (input/output) Matrix C (P x N). On exit, balanced C*D.
 * Not referenced if P=0 or if JOB implies C is not used.
 * @param ldc       Leading dimension of C.
 * @param[out] scale_out Scaling factors D(j) applied. Dimension (N).
 * @param row_major Specifies matrix storage (0 for column-major, 1 for row-major).
 *
 * @return info Error indicator: 0 for success, <0 if illegal argument.
 * SLICOT_MEMORY_ERROR (-1010) for internal memory allocation failure.
 */
SLICOT_EXPORT
int slicot_tb01id(
    char job_param,
    int n_param, int m_param, int p_param,
    double* maxred_io,
    double* a_io, int lda,
    double* b_io, int ldb,
    double* c_io, int ldc,
    double* scale_out,
    int row_major);

#ifdef __cplusplus
}
#endif

#endif /* SLICOT_WRAPPER_TB01ID_H */
