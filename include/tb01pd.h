/**
 * @file tb01pd.h
 * @brief Header for C wrapper of SLICOT routine TB01PD.
 */

#ifndef TB01PD_H
#define TB01PD_H

#include "slicot_utils.h" // Provides SLICOT_EXPORT macro

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Finds a reduced (controllable, observable, or minimal) state-space
 * representation (Ar,Br,Cr) for a given (A,B,C). Ar is upper block Hessenberg.
 * @details This is a C wrapper for the SLICOT Fortran routine TB01PD.
 * Matrices A, B, C are modified in place. NR is output.
 * Workspace is allocated internally.
 *
 * @param job_param   (input) Task: 'M' (minimal), 'C' (controllable), 'O' (observable).
 * @param equil_param (input) Balancing: 'S' (scale), 'N' (no scaling).
 * @param n_param     (input) Order of A. N >= 0.
 * @param m_param     (input) Columns of B (number of inputs). M >= 0.
 * @param p_param     (input) Rows of C (number of outputs). P >= 0.
 * @param a_io        (input/output) On entry, matrix A (N x N).
 * On exit, leading NR x NR part contains Ar.
 * @param lda         Leading dimension of A.
 * @param b_io        (input/output) On entry, matrix B.
 * For C row-major, array should be N rows by (M if JOB='C' else MAX(M,P)) columns.
 * On exit, leading NR x M part contains Br.
 * @param ldb         Leading dimension of B (columns for C row-major).
 * @param c_io        (input/output) On entry, matrix C (P x N).
 * For C row-major, array should be P rows by N columns.
 * (Fortran LDC is MAX(M,P) for workspace if JOB='M' or 'O').
 * On exit, leading P x NR part contains Cr.
 * @param ldc         Leading dimension of C (columns for C row-major, i.e., N).
 * @param[out] nr_out Order of the reduced system (Ar,Br,Cr).
 * @param tol_param   (input) Tolerance for rank determination. If <= 0, default is used.
 * @param row_major   Specifies matrix storage (0 for column-major, 1 for row-major).
 *
 * @return info Error indicator: 0 for success, <0 if illegal argument.
 * SLICOT_MEMORY_ERROR (-1010) for internal memory allocation failure.
 */
SLICOT_EXPORT
int slicot_tb01pd(
    char job_param, char equil_param,
    int n_param, int m_param, int p_param,
    double* a_io, int lda,
    double* b_io, int ldb,
    double* c_io, int ldc,
    int* nr_out,
    double tol_param,
    int row_major);

#ifdef __cplusplus
}
#endif

#endif /* TB01PD_H */
