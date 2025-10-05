#ifndef SB04ND_H
#define SB04ND_H

/**
 * @file sb04nd.h
 * @brief C wrapper for SLICOT routine SB04ND.
 *
 * Provides a C interface to the SLICOT routine SB04ND, which solves the
 * continuous-time Sylvester equation \f$ A X + X B = C \f$ assuming at least
 * one of the coefficient matrices is already in real Schur form and the other
 * is in Hessenberg or Schur form. The wrapper manages workspace allocation and
 * optional row-major inputs for C callers.
 */

#include "slicot_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Solve the Sylvester equation \f$ A X + X B = C \f$ with structured inputs.
 *
 * The routine assumes that matrix \f$A\f$ and/or \f$B\f$ is available in real
 * Schur form while the other is in Schur or Hessenberg form. Depending on the
 * forms specified by @p abschu, @p ula and @p ulb, the algorithm performs a block
 * back-substitution using the Hessenberg-Schur scheme described by Golub, Nash
 * and Van Loan. On exit, @p c contains the solution matrix \f$X\f$.
 *
 * @param[in] abschu  Mode flag indicating the structural form of \f$A\f$ and \f$B\f$:
 *                    - 'A': \f$A\f$ in Schur form, \f$B\f$ in (upper/lower) Hessenberg form.
 *                    - 'B': \f$B\f$ in Schur form, \f$A\f$ in (upper/lower) Hessenberg form.
 *                    - 'S': Both \f$A\f$ and \f$B\f$ in Schur form.
 * @param[in] ula     Indicates whether \f$A\f$ is upper ('U') or lower ('L') Schur/Hessenberg.
 * @param[in] ulb     Indicates whether \f$B\f$ is upper ('U') or lower ('L') Schur/Hessenberg.
 * @param[in] n       Order of \f$A\f$ and the number of rows in \f$C\f$ and \f$X\f$; @p n >= 0.
 * @param[in] m       Order of \f$B\f$ and the number of columns in \f$C\f$ and \f$X\f$; @p m >= 0.
 * @param[in,out] a   Double array holding \f$A\f\; dimension (lda, n) for column-major or (n, lda) for row-major.
 *                    Contents are preserved or updated according to the Fortran routine.
 * @param[in] lda     Leading dimension of @p a. For column-major storage, @p lda >= max(1, n);
 *                    for row-major, @p lda must be >= n when @p n > 0.
 * @param[in,out] b   Double array holding \f$B\f\; dimension (ldb, m) or (m, ldb) accordingly.
 *                    Contents may be modified by the routine.
 * @param[in] ldb     Leading dimension of @p b. For column-major storage, @p ldb >= max(1, m);
 *                    for row-major, @p ldb must be >= m when @p m > 0.
 * @param[in,out] c   Double array containing \f$C\f$ on entry; overwritten with the solution \f$X\f$ on exit.
 *                    Dimension (ldc, m) or (n, ldc) depending on @p row_major.
 * @param[in] ldc     Leading dimension of @p c. For column-major storage, @p ldc >= max(1, n);
 *                    for row-major, @p ldc must be >= m when both dimensions are non-zero.
 * @param[in] tol     Tolerance used to detect near-singularity when solving triangular systems.
 *                    If @p tol <= 0, the routine uses machine precision. This argument is
 *                    ignored when @p abschu = 'S' and @p ula = @p ulb = 'U'.
 * @param[in] row_major
 *                    Storage flag: 0 for column-major (Fortran-style) arrays, 1 for row-major
 *                    (C-style) arrays. The wrapper transposes data as needed when set to 1.
 *
 * @return int        Status code mirroring the Fortran routine:
 *                    - 0: successful exit;
 *                    - <0: if @c info = -i, the i-th argument had an illegal value;
 *                    - 1: a numerically singular block was encountered (reciprocal condition number <= tol).
 *                    The wrapper may return @c SLICOT_MEMORY_ERROR if workspace allocation fails.
 */
SLICOT_EXPORT
int slicot_sb04nd(char abschu, char ula, char ulb,
                  int n, int m,
                  double* a, int lda,
                  double* b, int ldb,
                  double* c, int ldc,
                  double tol,
                  int row_major);

#ifdef __cplusplus
}
#endif

#endif
