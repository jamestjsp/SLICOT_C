/**
 * @file tf01rd.h
 * @brief C wrapper for SLICOT routine TF01RD
 *
 * This file provides a C interface to the SLICOT routine TF01RD,
 * which computes N Markov parameters M(k) = C*A^(k-1)*B for a
 * linear time-invariant system (A,B,C).
 */

 #ifndef TF01RD_H
 #define TF01RD_H
 
 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Computes N Markov parameters M(k) = C*A^(k-1)*B for k=1..N.
  *
  * Calculates the first N Markov parameters for the system (A,B,C).
  * The parameters are stored concatenated in the output array H.
  *
  * @param[in] na        Order of matrix A, na >= 0.
  * @param[in] nb        Number of system inputs (columns of B), nb >= 0.
  * @param[in] nc        Number of system outputs (rows of C), nc >= 0.
  * @param[in] n         Number of Markov parameters M(k) to compute, n >= 0.
  * @param[in] a         Double array, dimension (lda, na) or (na, lda). State matrix A.
  * @param[in] lda       Leading dimension of A. >= max(1,na).
  * @param[in] b         Double array, dimension (ldb, nb) or (na, ldb). Input matrix B.
  * @param[in] ldb       Leading dimension of B. >= max(1,na).
  * @param[in] c         Double array, dimension (ldc, na) or (nc, ldc). Output matrix C.
  * @param[in] ldc       Leading dimension of C. >= max(1,nc).
  * @param[out] h        Double array, dimension (ldh, n*nb) or (nc, ldh).
  * Contains the concatenated Markov parameters [M(1) M(2) ... M(N)].
  * Each M(k) is NC x NB.
  * @param[in] ldh       Leading dimension of H. >= max(1,nc).
  * @param[in] row_major Integer flag: 0 for column-major, 1 for row-major.
  * Affects interpretation of lda, ldb, ldc, ldh.
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value
  * Memory allocation errors may also be returned.
  */
 SLICOT_EXPORT
 int slicot_tf01rd(int na, int nb, int nc, int n,
                   const double* a, int lda, const double* b, int ldb,
                   const double* c, int ldc, double* h, int ldh,
                   int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* TF01RD_H */
 