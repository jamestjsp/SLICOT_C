/**
 * @file sg03ad.h
 * @brief Header for C wrapper of SLICOT routine SG03AD.
 */

#ifndef SLICOT_WRAPPER_SG03AD_H
#define SLICOT_WRAPPER_SG03AD_H

#include "slicot_utils.h" // Provides SLICOT_EXPORT macro

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Solves generalized continuous- or discrete-time Lyapunov equations
 * and estimates separation.
 * @details This is a C wrapper for the SLICOT Fortran routine SG03AD.
 * Workspace is allocated internally.
 *
 * @param dico_param Specifies equation type: 'C' (continuous), 'D' (discrete).
 * @param job_param  Specifies task: 'X' (solve only), 'S' (separation only), 'B' (both).
 * @param fact_param Specifies if Schur factorization of (A,E) is supplied: 'N' (no), 'F' (yes).
 * @param trans_param Specifies if transposed equation is solved: 'N' (no), 'T' (yes).
 * @param uplo_param Specifies which triangle of input Y (in x_out) is used: 'L' (lower), 'U' (upper).
 * @param n_param    (input) Order of matrices A, E, Y, X. N >= 0.
 * @param a_io       (input/output) On entry, matrix A (N x N). If FACT='F', A_s.
 * On exit, Schur factor A_s.
 * @param lda        Leading dimension of A.
 * @param e_io       (input/output) On entry, matrix E (N x N). If FACT='F', E_s.
 * On exit, Schur factor E_s.
 * @param lde        Leading dimension of E.
 * @param q_io       (input/output) Orthogonal matrix Q from Schur factorization (N x N).
 * Input if FACT='F', output if FACT='N'.
 * @param ldq        Leading dimension of Q.
 * @param z_io       (input/output) Orthogonal matrix Z from Schur factorization (N x N).
 * Input if FACT='F', output if FACT='N'.
 * @param ldz        Leading dimension of Z.
 * @param y_rhs_in_x_out (input/output) On entry (if JOB='X' or 'B'), right-hand side Y (N x N).
 * On exit (if JOB='X' or 'B'), solution X (N x N).
 * @param ldx        Leading dimension of X/Y.
 * @param[out] scale_out Scale factor SCALE.
 * @param[out] sep_out   Estimate of separation SEP (if JOB='S' or 'B').
 * @param[out] ferr_out  Estimated forward error FERR (if JOB='B').
 * @param[out] alphar_out Real parts of eigenvalues of A - lambda*E (if FACT='N'). Dim (N).
 * @param[out] alfai_out Imaginary parts of eigenvalues (if FACT='N'). Dim (N).
 * @param[out] beta_out Denominators of eigenvalues (if FACT='N'). Dim (N).
 * @param row_major  Specifies matrix storage (0 for column-major, 1 for row-major).
 *
 * @return info Error indicator (see SLICOT documentation for SG03AD).
 * SLICOT_MEMORY_ERROR (-1010) for internal memory allocation failure.
 */
SLICOT_EXPORT
int slicot_sg03ad(
    char dico_param, char job_param, char fact_param, char trans_param, char uplo_param,
    int n_param,
    double* a_io, int lda,
    double* e_io, int lde,
    double* q_io, int ldq,
    double* z_io, int ldz,
    double* y_rhs_in_x_out, int ldx,
    double* scale_out, double* sep_out, double* ferr_out,
    double* alphar_out, double* alfai_out, double* beta_out,
    int row_major);

#ifdef __cplusplus
}
#endif

#endif /* SLICOT_WRAPPER_SG03AD_H */
