/**
 * @file tf01md.h
 * @brief Header for C wrapper of SLICOT routine TF01MD.
 */

 #ifndef TF01MD_H // Use unique guard
 #define TF01MD_H
 
 #include "slicot_utils.h" // Provides SLICOT_EXPORT macro, MAX, MIN etc.
 
 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Computes the output response sequence of a linear time-invariant discrete-time system.
  * @details Calculates y(k) = C*x(k) + D*u(k) and x(k+1) = A*x(k) + B*u(k)
  * for k = 1, ..., NY, given the state-space model (A,B,C,D), the initial
  * state vector x(1), and the input sequence u(1), ..., u(NY).
  * This is a C wrapper for the SLICOT Fortran routine TF01MD.
  * **Workspace is allocated internally.**
  *
  * @param[in] n         Order of the state matrix A. n >= 0.
  * @param[in] m         Number of system inputs. m >= 0.
  * @param[in] p         Number of system outputs. p >= 0.
  * @param[in] ny        Number of output vectors y(k) to be computed. ny >= 0.
  * @param[in] a         Pointer to the state matrix A. Dimensions (n x n).
  * Stored column-wise if row_major=0, row-wise if row_major=1.
  * @param[in] lda       Leading dimension of the C array storing A.
  * If row_major=0 (column-major), lda >= max(1, n) (number of rows).
  * If row_major=1 (row-major), lda >= max(1, n) (number of columns).
  * @param[in] b         Pointer to the input matrix B. Dimensions (n x m).
  * Stored column-wise if row_major=0, row-wise if row_major=1.
  * @param[in] ldb       Leading dimension of the C array storing B.
  * If row_major=0 (column-major), ldb >= max(1, n) (number of rows).
  * If row_major=1 (row-major), ldb >= max(1, m) (number of columns).
  * @param[in] c         Pointer to the output matrix C. Dimensions (p x n).
  * Stored column-wise if row_major=0, row-wise if row_major=1.
  * @param[in] ldc       Leading dimension of the C array storing C.
  * If row_major=0 (column-major), ldc >= max(1, p) (number of rows).
  * If row_major=1 (row-major), ldc >= max(1, n) (number of columns).
  * @param[in] d         Pointer to the direct link matrix D. Dimensions (p x m).
  * Stored column-wise if row_major=0, row-wise if row_major=1.
  * @param[in] ldd       Leading dimension of the C array storing D.
  * If row_major=0 (column-major), ldd >= max(1, p) (number of rows).
  * If row_major=1 (row-major), ldd >= max(1, m) (number of columns).
  * @param[in] u         Pointer to the input sequence matrix U. Dimensions (m x ny).
  * Stored column-wise if row_major=0, row-wise if row_major=1.
  * The k-th column (or row if row_major=1) contains u(k).
  * @param[in] ldu       Leading dimension of the C array storing U.
  * If row_major=0 (column-major), ldu >= max(1, m) (number of rows).
  * If row_major=1 (row-major), ldu >= max(1, ny) (number of columns).
  * @param[in,out] x     Pointer to the state vector. Dimension (n).
  * On entry, this array must contain the initial state vector x(1).
  * On exit, this array contains the final state vector x(NY+1).
  * @param[out] y        Pointer to the output sequence matrix Y. Dimensions (p x ny).
  * Stored column-wise if row_major=0, row-wise if row_major=1.
  * The k-th column (or row if row_major=1) will contain y(k).
  * @param[in] ldy       Leading dimension of the C array storing Y.
  * If row_major=0 (column-major), ldy >= max(1, p) (number of rows).
  * If row_major=1 (row-major), ldy >= max(1, ny) (number of columns).
  * @param[in] row_major Specifies matrix storage for A, B, C, D, U, Y:
  * 0 for column-major (Fortran default),
  * 1 for row-major (C default).
  *
  * @return info Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value (wrapper or Fortran validation)
  * -1: n < 0
  * -2: m < 0
  * -3: p < 0
  * -4: ny < 0
  * -5: a == NULL and n > 0
  * -6: lda illegal
  * -7: b == NULL and n > 0 and m > 0
  * -8: ldb illegal
  * -9: c == NULL and p > 0 and n > 0
  * -10: ldc illegal
  * -11: d == NULL and p > 0 and m > 0
  * -12: ldd illegal
  * -13: u == NULL and m > 0 and ny > 0
  * -14: ldu illegal
  * -15: x == NULL and n > 0
  * -16: y == NULL and p > 0 and ny > 0
  * -17: ldy illegal
  * = SLICOT_MEMORY_ERROR (-1010): internal memory allocation failed.
  */
 SLICOT_EXPORT
 int slicot_tf01md(int n, int m, int p, int ny,
                   const double* a, int lda, const double* b, int ldb,
                   const double* c, int ldc, const double* d, int ldd,
                   const double* u, int ldu,
                   double* x, double* y, int ldy, int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* TF01MD_H */
 