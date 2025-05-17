/**
 * @file sg02ad.h
 * @brief Header for C wrapper of SLICOT routine SG02AD.
 */

#ifndef SLICOT_WRAPPER_SG02AD_H
#define SLICOT_WRAPPER_SG02AD_H

#include "slicot_utils.h" // Provides SLICOT_EXPORT macro

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Solves continuous- or discrete-time algebraic Riccati equations for descriptor systems.
 * @details This is a C wrapper for the SLICOT Fortran routine SG02AD.
 * Workspace is allocated internally.
 *
 * @param dico_param Specifies continuous ('C') or discrete-time ('D').
 * @param jobb_param Specifies if B, R ('B') or G ('G') are given.
 * @param fact_param Specifies if Q, R are factored ('N', 'C', 'D', 'B').
 * @param uplo_param Specifies if upper ('U') or lower ('L') triangle is stored for symmetric Q, R, G.
 * @param jobl_param Specifies if L is zero ('Z') or nonzero ('N').
 * @param scal_param Specifies scaling strategy ('G', 'N'). Not used if JOBB='G'.
 * @param sort_param Specifies eigenvalue sorting ('S' for stable, 'U' for unstable).
 * @param acc_param  Specifies iterative refinement for X ('R', 'N').
 * @param n_param    (input) Order of matrices A, E, Q, X. N >= 0.
 * @param m_param    (input) Number of system inputs. M >= 0. Not used if JOBB='G'.
 * @param p_param    (input) Number of system outputs (used if FACT is 'C','D','B'). P >= 0.
 * @param a_io       (input) State matrix A (N x N).
 * @param lda        Leading dimension of A.
 * @param e_io       (input) Descriptor matrix E (N x N).
 * @param lde        Leading dimension of E.
 * @param b_io       (input) If JOBB='B', input matrix B (N x M). If JOBB='G', matrix G (N x N, triangular).
 * @param ldb        Leading dimension of B/G.
 * @param q_io       (input) If FACT='N'/'D', matrix Q (N x N, triangular). If FACT='C'/'B', factor C (P x N).
 * @param ldq        Leading dimension of Q/C_factor.
 * @param r_io       (input) If JOBB='B': if FACT='N'/'C', matrix R (M x M, triangular). If FACT='D'/'B', factor D (P x M). Not used if JOBB='G'.
 * @param ldr        Leading dimension of R/D_factor.
 * @param l_io       (input) If JOBB='B' and JOBL='N', matrix L (N x M). Not used otherwise.
 * @param ldl        Leading dimension of L.
 * @param[out] rcondu Estimate of reciprocal condition number for system solving for X.
 * @param[out] x_out  Solution matrix X (N x N).
 * @param ldx        Leading dimension of X.
 * @param[out] alfar_out Real parts of generalized eigenvalues (2*N).
 * @param[out] alfai_out Imaginary parts of generalized eigenvalues (2*N).
 * @param[out] beta_out Denominators of generalized eigenvalues (2*N).
 * @param[out] s_out  Schur form S. Dim: LDS x (2N+M if JOBB='B', else 2N).
 * @param lds        Leading dimension of S.
 * @param[out] t_out  Upper triangular T. Dim: LDT x 2N.
 * @param ldt        Leading dimension of T.
 * @param[out] u_out  Transformation matrix U. Dim: LDU x 2N.
 * @param ldu        Leading dimension of U.
 * @param tol_param  (input) Tolerance for rank determination.
 * @param[out] iwarn_out Warning indicator.
 * @param row_major  Specifies matrix storage (0 for column-major, 1 for row-major).
 *
 * @return info Error indicator (see SLICOT documentation for SG02AD).
 * SLICOT_MEMORY_ERROR (-1010) for internal memory allocation failure.
 */
SLICOT_EXPORT
int slicot_sg02ad(
    char dico_param, char jobb_param, char fact_param, char uplo_param,
    char jobl_param, char scal_param, char sort_param, char acc_param,
    int n_param, int m_param, int p_param,
    double* a_io, int lda,
    double* e_io, int lde,
    double* b_io, int ldb,
    double* q_io, int ldq,
    double* r_io, int ldr,
    double* l_io, int ldl,
    double* rcondu,
    double* x_out, int ldx,
    double* alfar_out, double* alfai_out, double* beta_out,
    double* s_out, int lds,
    double* t_out, int ldt,
    double* u_out, int ldu,
    double tol_param,
    int* iwarn_out,
    int row_major);

#ifdef __cplusplus
}
#endif

#endif /* SLICOT_WRAPPER_SG02AD_H */
