/**
 * @file sb02mt.h
 * @brief C wrapper for SLICOT routine SB02MT
 *
 * This file provides a C interface to the SLICOT routine SB02MT,
 * which converts optimal problems with coupling weighting terms to
 * standard problems.
 */

#ifndef SB02MT_H
#define SB02MT_H

#include <stddef.h> // For size_t
#include "slicot_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Converts optimal problems with coupling weighting terms to standard problems.
 *
 * Computes the matrices:
 * G_out = B * R^(-1) * B'  (if jobg == 'G')
 * A_io  = A_io - B_io * R_io^(-1) * L_io' (if jobl == 'N')
 * Q_io  = Q_io - L_io * R_io^(-1) * L_io' (if jobl == 'N')
 *
 * Workspace (IWORK, DWORK) is allocated internally by this wrapper.
 * Input/output matrix format is handled via the row_major parameter.
 *
 * @param[in] jobg      Specifies whether or not the matrix G is to be computed:
 * = 'G': Compute G_out.
 * = 'N': Do not compute G_out.
 * @param[in] jobl      Specifies whether or not the matrix L_io is zero:
 * = 'Z': L_io is zero (A_io and Q_io are not modified by L_io terms).
 * = 'N': L_io is non-zero.
 * @param[in] fact      Specifies how the matrix R_io is given:
 * = 'N': R_io contains the matrix R.
 * = 'C': R_io contains the Cholesky factor of R.
 * = 'U': R_io contains the UdU' or LdL' factorization of R.
 * @param[in] uplo      Specifies which triangle of symmetric matrices R_io, Q_io, G_out is stored:
 * = 'U': Upper triangle.
 * = 'L': Lower triangle.
 * @param[in] n         Order of matrices A_io, Q_io, G_out; number of rows of B_io, L_io. N >= 0.
 * @param[in] m         Order of matrix R_io; number of columns of B_io, L_io. M >= 0.
 * @param[in,out] a_io  On entry (if jobl='N'), the N-by-N matrix A.
 * On exit (if jobl='N', info=0), the updated N-by-N matrix A_io.
 * Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param[in] lda       Leading dimension of a_io.
 * If row_major=0: lda >= max(1,N) if jobl='N'; lda >= 1 if jobl='Z'.
 * If row_major=1: lda >= max(1,N) if jobl='N' (cols); lda >= 1 if jobl='Z'.
 * @param[in,out] b_io  On entry, the N-by-M matrix B.
 * On exit (if oufact_out=1, info=0), B_io contains B*chol(R)^(-1). Unchanged if oufact_out != 1.
 * Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param[in] ldb       Leading dimension of b_io.
 * If row_major=0: ldb >= max(1,N).
 * If row_major=1: ldb >= max(1,M) (cols).
 * @param[in,out] q_io  On entry (if jobl='N'), the N-by-N symmetric matrix Q.
 * On exit (if jobl='N', info=0), the updated N-by-N symmetric matrix Q_io.
 * Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param[in] ldq       Leading dimension of q_io.
 * If row_major=0: ldq >= max(1,N) if jobl='N'; ldq >= 1 if jobl='Z'.
 * If row_major=1: ldq >= max(1,N) if jobl='N' (cols); ldq >= 1 if jobl='Z'.
 * @param[in,out] r_io  On entry, the M-by-M symmetric matrix R or its factors (see fact).
 * On exit (if oufact_out=1 or 2, info=0 or M+1), R_io contains computed factors. Unchanged if fact='C' or 'U'.
 * Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param[in] ldr       Leading dimension of r_io.
 * If row_major=0: ldr >= max(1,M).
 * If row_major=1: ldr >= max(1,M) (cols).
 * @param[in,out] l_io  On entry (if jobl='N'), the N-by-M matrix L.
 * On exit (if jobl='N', oufact_out=1, info=0), L_io contains L*chol(R)^(-1). Unchanged if oufact_out != 1.
 * Not referenced if jobl='Z'.
 * Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param[in] ldl       Leading dimension of l_io.
 * If row_major=0: ldl >= max(1,N) if jobl='N'; ldl >= 1 if jobl='Z'.
 * If row_major=1: ldl >= max(1,M) if jobl='N' (cols); ldl >= 1 if jobl='Z'.
 * @param[in,out] ipiv_io Integer array, dimension (M).
 * On entry (if fact='U'), details of interchanges for R.
 * On exit (if oufact_out=2), details of interchanges for R.
 * Not referenced if fact='C'.
 * @param[out] g_out    If jobg='G' and info=0, the N-by-N symmetric matrix G.
 * Not referenced if jobg='N'.
 * Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param[in] ldg       Leading dimension of g_out.
 * If row_major=0: ldg >= max(1,N) if jobg='G'; ldg >= 1 if jobg='N'.
 * If row_major=1: ldg >= max(1,N) if jobg='G' (cols); ldg >= 1 if jobg='N'.
 * @param[out] oufact_out Information about the factorization of R used:
 * =0: M=0, no factorization.
 * =1: Cholesky factorization.
 * =2: UdU' or LdL' factorization.
 * @param[in] row_major Specifies matrix storage: 0 for column-major, 1 for row-major.
 *
 * @return info         Error indicator:
 * = 0: successful exit.
 * < 0: if info = -i, the i-th argument had an illegal value.
 * = i (1<=i<=M): i-th element of d factor in R's factorization is zero.
 * = M+1: R is numerically singular.
 * = SLICOT_MEMORY_ERROR (-1010): internal memory allocation failed.
 */
SLICOT_EXPORT
int slicot_sb02mt(
    char jobg, char jobl, char fact, char uplo,
    int n, int m,
    double* a_io, int lda,
    double* b_io, int ldb,
    double* q_io, int ldq,
    double* r_io, int ldr,
    double* l_io, int ldl,
    int* ipiv_io,
    double* g_out, int ldg,
    int* oufact_out,
    int row_major
);

#ifdef __cplusplus
}
#endif

#endif /* SB02MT_H */