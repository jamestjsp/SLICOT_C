/**
 * @file mb02ed.h
 * @brief C wrapper for SLICOT routine MB02ED
 *
 * This file provides a C interface to the SLICOT routine MB02ED,
 * which solves T*X = B or X*T = B for a symmetric positive definite
 * block Toeplitz matrix T.
 */

#ifndef MB02ED_H
#define MB02ED_H

 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Solves T*X = B or X*T = B for a block Toeplitz matrix T.
 *
 * Solves a system of linear equations where T is a symmetric positive
 * definite (s.p.d.) block Toeplitz matrix, defined by its first
 * block row or column.
 *
 * @param[in] typet     Specifies the form of T and the system to solve:
 * = 'R': T array contains first block row, solves X*T = B.
 * = 'C': T array contains first block column, solves T*X = B.
 * @param[in] k         The block size (rows/columns of blocks in T), k >= 0.
 * @param[in] n         The number of blocks in T, n >= 0.
 * @param[in] nrhs      The number of right hand sides (columns of B if typet='C', rows of B if typet='R'), nrhs >= 0.
 * @param[in,out] t     Double array.
 * If typet='R', dimension (ldt, n*k) or (k, ldt). Contains first block row on entry.
 * If typet='C', dimension (ldt, k) or (n*k, ldt). Contains first block column on entry.
 * On exit, contains information about the Cholesky factor of inv(T).
 * @param[in] ldt       Leading dimension of T.
 * If typet='R', ldt >= max(1,k).
 * If typet='C', ldt >= max(1,n*k).
 * @param[in,out] b     Double array containing the right hand side B on entry
 * and the solution X on exit.
 * If typet='R', dimension (ldb, n*k) or (nrhs, ldb).
 * If typet='C', dimension (ldb, nrhs) or (n*k, ldb).
 * @param[in] ldb       Leading dimension of B.
 * If typet='R', ldb >= max(1,nrhs).
 * If typet='C', ldb >= max(1,n*k).
 * @param[in] row_major Integer flag:
 * = 0: Arrays t, b are column-major (Fortran style).
 * = 1: Arrays t, b are row-major (C style).
 *
 * @return info         Error indicator:
 * = 0: successful exit
 * < 0: if info = -i, the i-th argument had an illegal value.
 * = 1: Reduction algorithm failed; T is not numerically positive definite.
 * Memory allocation errors may also be returned.
 */
SLICOT_EXPORT
int slicot_mb02ed(char typet, int k, int n, int nrhs,
                  double* t, int ldt, double* b, int ldb,
                  int row_major);

#ifdef __cplusplus
}
#endif

#endif /* MB02ED_H */
