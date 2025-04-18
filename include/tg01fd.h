/**
 * @file tg01fd.h
 * @brief C wrapper for SLICOT routine TG01FD
 *
 * This file provides a C interface to the SLICOT routine TG01FD,
 * which computes an orthogonal reduction of a descriptor system
 * (A-lambda E,B,C) to a SVD-like coordinate form.
 */

 #ifndef TG01FD_H
 #define TG01FD_H
 
 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Reduces a descriptor system (A-lambda*E, B, C) to SVD-like coordinate form.
  *
  * Computes orthogonal Q, Z such that Q'*E*Z = [Er 0; 0 0] and Q'*A*Z, Q'*B, C*Z
  * are transformed accordingly. Er is upper triangular and invertible.
  * Optionally reduces the (2,2) block of transformed A further.
  *
  * @param[in] compq     Compute Q: 'N' (no), 'I' (initialize), 'U' (update Q1 -> Q1*Q).
  * @param[in] compz     Compute Z: 'N' (no), 'I' (initialize), 'U' (update Z1 -> Z1*Z).
  * @param[in] joba      Reduce A22 block: 'N' (no), 'R' (SVD-like triangular), 'T' (upper trapezoidal).
  * @param[in] l         Number of rows of A, E, B, l >= 0.
  * @param[in] n         Number of columns of A, E, C, n >= 0.
  * @param[in] m         Number of columns of B, m >= 0.
  * @param[in] p         Number of rows of C, p >= 0.
  * @param[in,out] a     Double array, dimension (lda, n) or (l, lda). State matrix A. Overwritten by Q'*A*Z.
  * @param[in] lda       Leading dimension of A. >= max(1,l).
  * @param[in,out] e     Double array, dimension (lde, n) or (l, lde). Descriptor matrix E. Overwritten by Q'*E*Z.
  * @param[in] lde       Leading dimension of E. >= max(1,l).
  * @param[in,out] b     Double array, dimension (ldb, m) or (l, ldb). Input matrix B. Overwritten by Q'*B. Not referenced if m=0.
  * @param[in] ldb       Leading dimension of B. >= max(1,l) if m>0, else >=1.
  * @param[in,out] c     Double array, dimension (ldc, n) or (p, ldc). Output matrix C. Overwritten by C*Z. Not referenced if p=0.
  * @param[in] ldc       Leading dimension of C. >= max(1,p).
  * @param[in,out] q     Double array, dimension (ldq, l) or (l, ldq). Orthogonal matrix Q. See compq.
  * @param[in] ldq       Leading dimension of Q. >=1 if compq='N', >=max(1,l) otherwise.
  * @param[in,out] z     Double array, dimension (ldz, n) or (n, ldz). Orthogonal matrix Z. See compz.
  * @param[in] ldz       Leading dimension of Z. >=1 if compz='N', >=max(1,n) otherwise.
  * @param[out] ranke    Estimated rank of E.
  * @param[out] rnka22   Estimated rank of A22 block (if joba='R' or 'T').
  * @param[in] tol       Tolerance for rank determination. <= 0 for default. tol < 1.
  * @param[in] row_major Integer flag: 0 for column-major, 1 for row-major.
  * Affects interpretation of lda, lde, ldb, ldc, ldq, ldz.
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value
  * Memory allocation errors may also be returned.
  */
 SLICOT_C_WRAPPER_API
 int slicot_tg01fd(char compq, char compz, char joba, int l, int n, int m, int p,
                   double* a, int lda, double* e, int lde,
                   double* b, int ldb, double* c, int ldc,
                   double* q, int ldq, double* z, int ldz,
                   int* ranke, int* rnka22, double tol, int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* TG01FD_H */
 