/**
 * @file tb01id.h
 * @brief C wrapper for SLICOT routine TB01ID
 *
 * This file provides a C interface to the SLICOT routine TB01ID,
 * which balances a system matrix S = [A B; C 0] corresponding to a
 * state-space triplet (A,B,C) using diagonal similarity transformations.
 */

 #ifndef TB01ID_H
 #define TB01ID_H
 
 #include <stddef.h> // For size_t
 
 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Balances a system matrix S = [A B; C 0].
  *
  * Reduces the 1-norm of the system matrix S by applying a diagonal
  * similarity transformation diag(inv(D), I) * S * diag(D, I).
  * This aims to make the 1-norms of corresponding rows and columns
  * (for the first N rows/columns related to A) nearly equal.
  * Optionally balances only A, (A,B), or (A,C).
  *
  * @param[in] job       Specifies which matrices are involved:
  * = 'A': Balance S = [A B; C 0].
  * = 'B': Balance S = [A B].
  * = 'C': Balance S = [A; C].
  * = 'N': Balance S = A.
  * @param[in] n         Order of A (rows of B, columns of C), n >= 0.
  * @param[in] m         Number of columns of B, m >= 0.
  * @param[in] p         Number of rows of C, p >= 0.
  * @param[in,out] maxred On entry, max allowed reduction ratio (> 1.0) if zeros found. If <=0, default 10.0 used.
  * On exit, ratio of 1-norm before vs after balancing.
  * @param[in,out] a     Double array, dimension (lda, n) or (n, lda). Matrix A. Overwritten with inv(D)*A*D.
  * @param[in] lda       Leading dimension of A. >= max(1,n).
  * @param[in,out] b     Double array, dimension (ldb, m) or (n, ldb). Matrix B. Overwritten with inv(D)*B. Not referenced if m=0.
  * @param[in] ldb       Leading dimension of B. >=max(1,n) if m>0, else >=1.
  * @param[in,out] c     Double array, dimension (ldc, n) or (p, ldc). Matrix C. Overwritten with C*D. Not referenced if p=0.
  * @param[in] ldc       Leading dimension of C. >= max(1,p).
  * @param[out] scale    Double array, dimension (n). Contains the diagonal elements of D.
  * @param[in] row_major Integer flag:
  * = 0: Arrays a, b, c are column-major (Fortran style).
  * = 1: Arrays a, b, c are row-major (C style).
  * (Array scale is 1D).
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * Memory allocation errors may also be returned by the wrapper.
  */
 int slicot_tb01id(char job, int n, int m, int p, double* maxred,
                   double* a, int lda, double* b, int ldb,
                   double* c, int ldc, double* scale,
                   int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* TB01ID_H */
 