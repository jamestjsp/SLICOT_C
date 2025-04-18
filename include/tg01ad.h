/**
 * @file tg01ad.h
 * @brief C wrapper for SLICOT routine TG01AD
 *
 * This file provides a C interface to the SLICOT routine TG01AD,
 * which balances the matrices of the system pencil corresponding to a
 * descriptor triple (A-lambda E,B,C) using diagonal similarity
 * transformations.
 */

 #ifndef TG01AD_H
 #define TG01AD_H
 
 #include <stddef.h> // For size_t
 
 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Balances the system pencil S = [A B; C 0] - lambda*[E 0; 0 0].
  *
  * Applies diagonal similarity transformations diag(Dl,I)*S*diag(Dr,I)
  * to make row and column 1-norms of the pencil matrices nearly equal.
  * Can operate on subsets of the pencil (e.g., A-lambda*E only).
  *
  * @param[in] job       Specifies which matrices are involved:
  * 'A': All matrices (A, E, B, C).
  * 'B': Matrices A, E, B.
  * 'C': Matrices A, E, C.
  * 'N': Matrices A, E only.
  * @param[in] l         Number of rows of A, E, B, l >= 0.
  * @param[in] n         Number of columns of A, E, C, n >= 0.
  * @param[in] m         Number of columns of B, m >= 0.
  * @param[in] p         Number of rows of C, p >= 0.
  * @param[in] thresh    Threshold for ignoring small elements, thresh >= 0.
  * @param[in,out] a     Double array, dimension (lda, n) or (l, lda). State matrix A. Overwritten by Dl*A*Dr.
  * @param[in] lda       Leading dimension of A. >= max(1,l).
  * @param[in,out] e     Double array, dimension (lde, n) or (l, lde). Descriptor matrix E. Overwritten by Dl*E*Dr.
  * @param[in] lde       Leading dimension of E. >= max(1,l).
  * @param[in,out] b     Double array, dimension (ldb, m) or (l, ldb). Input matrix B. Overwritten by Dl*B. Not referenced if m=0.
  * @param[in] ldb       Leading dimension of B. >= max(1,l) if m>0, else >=1.
  * @param[in,out] c     Double array, dimension (ldc, n) or (p, ldc). Output matrix C. Overwritten by C*Dr. Not referenced if p=0.
  * @param[in] ldc       Leading dimension of C. >= max(1,p).
  * @param[out] lscale   Double array, dimension (l). Left scaling factors Dl.
  * @param[out] rscale   Double array, dimension (n). Right scaling factors Dr.
  * @param[in] row_major Integer flag: 0 for column-major, 1 for row-major.
  * Affects interpretation of lda, lde, ldb, ldc.
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value
  * Memory allocation errors may also be returned.
  */
 int slicot_tg01ad(char job, int l, int n, int m, int p, double thresh,
                   double* a, int lda, double* e, int lde,
                   double* b, int ldb, double* c, int ldc,
                   double* lscale, double* rscale, int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* TG01AD_H */
 