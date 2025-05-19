/**
 * @file sg03bd.h
 * @brief Header for C wrapper of SLICOT routine SG03BD.
 */

#ifndef SG03BD_H
#define SG03BD_H

#include "slicot_utils.h" // Provides SLICOT_EXPORT macro

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Solves for Cholesky factor U of the solution X to generalized stable
 * continuous- or discrete-time Lyapunov equations. X = op(U)' * op(U).
 * @details This is a C wrapper for the SLICOT Fortran routine SG03BD.
 * Workspace is allocated internally.
 *
 * @param dico_param Specifies equation type: 'C' (continuous), 'D' (discrete).
 * @param fact_param Specifies if Schur factorization of (A,E) is supplied: 'N' (no), 'F' (yes).
 * @param trans_param Specifies if transposed equation is solved: 'N' (op(M)=M), 'T' (op(M)=M').
 * @param n_param    (input) Order of matrix A. N >= 0.
 * @param m_param    (input) Number of rows in op(B) (matrix from RHS). M >= 0.
 * @param a_io       (input/output) On entry, matrix A (N x N). If FACT='F', A_s.
 * On exit (if FACT='N'), Schur factor A_s.
 * @param lda        Leading dimension of A.
 * @param e_io       (input/output) On entry, matrix E (N x N). If FACT='F', E_s.
 * On exit (if FACT='N'), Schur factor E_s.
 * @param lde        Leading dimension of E.
 * @param q_io       (input/output) Orthogonal matrix Q from Schur factorization (N x N).
 * Input if FACT='F', output if FACT='N'.
 * @param ldq        Leading dimension of Q.
 * @param z_io       (input/output) Orthogonal matrix Z from Schur factorization (N x N).
 * Input if FACT='F', output if FACT='N'.
 * @param ldz        Leading dimension of Z.
 * @param b_in_u_out (input/output) On entry, contains matrix B for RHS.
 * Dimensions of input B: if TRANS='T', N x M; if TRANS='N', M x N.
 * On exit, the leading N x N part contains Cholesky factor U.
 * The C array must be dimensioned to hold an N x N matrix for output U.
 * @param ldb        Leading dimension of the C array b_in_u_out.
 * For row_major=1, ldb must be >= N (cols of output U).
 * For row_major=0, ldb must be >= MAX(1, N, (TRANS='N' ? M : N_if_TRANS_T_and_N>M)).
 * @param[out] scale_out Scale factor SCALE.
 * @param[out] alphar_out Real parts of eigenvalues of A - lambda*E. Dim (N).
 * @param[out] alfai_out Imaginary parts of eigenvalues. Dim (N).
 * @param[out] beta_out Denominators of eigenvalues. Dim (N).
 * @param row_major  Specifies matrix storage (0 for column-major, 1 for row-major).
 *
 * @return info Error indicator (see SLICOT documentation for SG03BD).
 * SLICOT_MEMORY_ERROR (-1010) for internal memory allocation failure.
 */
SLICOT_EXPORT
int slicot_sg03bd(
    char dico_param, char fact_param, char trans_param,
    int n_param, int m_param,
    double* a_io, int lda,
    double* e_io, int lde,
    double* q_io, int ldq,
    double* z_io, int ldz,
    double* b_in_u_out, int ldb,
    double* scale_out,
    double* alphar_out, double* alfai_out, double* beta_out,
    int row_major);

#ifdef __cplusplus
}
#endif

#endif /* SG03BD_H */
