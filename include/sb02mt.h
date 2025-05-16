/**
 * @file sb02mt.h
 * @brief C wrapper for SLICOT routine SB02MT
 *
 * This file provides a C interface to the SLICOT routine SB02MT,
 * which converts optimal problems with coupling weighting terms to standard problems.
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
 *   G = B*R^(-1)*B'
 *   A_bar = A - B*R^(-1)*L'
 *   Q_bar = Q - L*R^(-1)*L'
 *
 * @param[in] dico      Character specifying the type of the system:
 *                      = 'C': Continuous-time system
 *                      = 'D': Discrete-time system
 * @param[in] jobb      Character specifying whether or not B is non-zero:
 *                      = 'B': B is non-zero
 *                      = 'Z': B is zero
 * @param[in] fact      Character specifying whether R is factored:
 *                      = 'N': R is not factored
 *                      = 'C': R contains Cholesky factor
 *                      = 'U': R contains UdU'/LdL' factorization
 * @param[in] uplo      Character specifying which triangle is stored:
 *                      = 'U': Upper triangle is stored
 *                      = 'L': Lower triangle is stored
 * @param[in] jobl      Character specifying whether L is zero:
 *                      = 'Z': L is zero
 *                      = 'N': L is non-zero
 * @param[in] sort      Character specifying ordering of eigenvalues:
 *                      = 'S': Stable eigenvalues first
 *                      = 'U': Unstable eigenvalues first
 * @param[in] n         Order of the state matrix A, n >= 0
 * @param[in] m         Number of inputs, m >= 0
 * @param[in] p         Number of outputs, p >= 0
 * @param[in,out] a     State matrix A, dimension (lda,n) or (n,lda)
 *                      On output if jobl='N': A_bar = A - B*R^(-1)*L'
 * @param[in] lda       Leading dimension of a, >= max(1,n) if jobl='N', >= 1 otherwise
 * @param[in,out] b     Input matrix B, dimension (ldb,m) or (n,ldb)
 * @param[in] ldb       Leading dimension of b, >= max(1,n)
 * @param[in,out] q     Input symmetric matrix Q, dimension (ldq,n) or (n,ldq)
 *                      On output if jobl='N': Q_bar = Q - L*R^(-1)*L'
 * @param[in] ldq       Leading dimension of q, >= max(1,n) if jobl='N', >= 1 otherwise
 * @param[in,out] r     Input matrix R, dimension (ldr,m) or (m,ldr)
 * @param[in] ldr       Leading dimension of r, >= max(1,m)
 * @param[in,out] l     Cross-term weighting matrix L, dimension (ldl,m) or (n,ldl)
 *                      Only used if jobl='N'
 * @param[in] ldl       Leading dimension of l, >= max(1,n) if jobl='N', >= 1 otherwise
 * @param[out] g        Output matrix G = B*R^(-1)*B', dimension (ldg,n) or (n,ldg)
 *                      Only computed if jobg='G'
 * @param[in] ldg       Leading dimension of g, >= max(1,n) if jobg='G', >= 1 otherwise
 * @param[out] rcond    Array containing reciprocal condition numbers:
 *                      rcond[0] = Reciprocal condition number of R
 *                      rcond[1] = Reciprocal condition number used in rank determination
 * @param[out] rank     Integer giving the effective rank of matrix R
 * @param[in] tol       Tolerance used for numerical rank determination
 * @param[in] row_major Flag indicating storage format:
 *                      = 0: Column-major (Fortran) storage
 *                      = 1: Row-major (C) storage
 *
 * @return info         Status code:
 *                      = 0: successful exit
 *                      < 0: if info = -i, the i-th argument had an illegal value
 *                      = i: i-th diagonal element of d factor is zero (1<=i<=m)
 *                      = m+1: R is numerically singular
 */
SLICOT_EXPORT
int slicot_sb02mt(
    char dico, char jobb, char fact, char uplo, char jobl, char sort,
    int n, int m, int p,
    double* a, int lda,     // Input/Output A
    double* b, int ldb,     // Input B
    double* q, int ldq,     // Input/Output Q
    double* r, int ldr,     // Input R
    double* l, int ldl,     // Input L (if JOBL='N')
    double* x, int ldx,     // Output X (solution)
    double* g, int ldg,     // Output G (gain)
    double* rcond,          // Output RCOND (array of 2)
    int* rank,              // Output RANK
    double* s, int lds,     // Output S (Schur matrix)
    double* u, int ldu,     // Output U (transformation matrix)
    double* wr,             // Output WR (real part of eigenvalues)
    double* wi,             // Output WI (imaginary part of eigenvalues)
    double tol,             // Input TOL
    int row_major           // Input row_major
);

#ifdef __cplusplus
}
#endif

#endif /* SB02MT_H */
