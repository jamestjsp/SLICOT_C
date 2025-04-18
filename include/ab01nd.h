/**
 * @file ab01nd.h
 * @brief C wrapper for SLICOT routine AB01ND
 *
 * This file provides a C interface to the SLICOT routine AB01ND,
 * which finds a controllable realization for a linear time-invariant
 * multi-input system using orthogonal transformations.
 */

 #ifndef AB01ND_H
 #define AB01ND_H
 
 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Finds a controllable realization for a multi-input system (A,B).
  *
  * Reduces the system (A,B) using orthogonal similarity transformations Z
  * such that Ac = Z'*A*Z is upper block Hessenberg and Bc = Z'*B has
  * non-zero elements only in its first block B1.
  * Determines the order NCONT and controllability index INDCON of the
  * controllable subsystem and the block sizes NBLK.
  *
  * @param[in] jobz      Specifies whether to compute the transformation matrix Z:
  * = 'N': Do not compute Z.
  * = 'F': Compute Z in factored form (reflectors stored in Z and TAU).
  * = 'I': Compute Z explicitly (initialize Z to identity before call).
  * @param[in] n         The order of the matrix A, n >= 0.
  * @param[in] m         The number of system inputs (columns of B), m >= 0.
  * @param[in,out] a     Double array, dimension (lda, n) (col-major) or (n, lda) (row-major).
  * On entry, the state dynamics matrix A.
  * On exit, the leading NCONT-by-NCONT part contains the
  * transformed upper block Hessenberg matrix Ac.
  * @param[in] lda       The leading dimension of array A.
  * If row_major=0, lda >= max(1,n).
  * If row_major=1, lda >= max(1,n) (number of columns).
  * @param[in,out] b     Double array, dimension (ldb, m) (col-major) or (n, ldb) (row-major).
  * On entry, the input matrix B.
  * On exit, the leading NCONT-by-M part contains the
  * transformed input matrix Bc.
  * @param[in] ldb       The leading dimension of array B.
  * If row_major=0, ldb >= max(1,n).
  * If row_major=1, ldb >= max(1,m) (number of columns).
  * @param[out] ncont    The order of the controllable state-space representation.
  * @param[out] indcon   The controllability index of the controllable part.
  * @param[out] nblk     Integer array, dimension (n).
  * The leading INDCON elements contain the orders of the
  * diagonal blocks of the transformed matrix Ac.
  * @param[out] z        Double array, dimension (ldz, n) (col-major) or (n, ldz) (row-major).
  * If jobz = 'I', contains the N-by-N orthogonal transformation matrix Z.
  * If jobz = 'F', contains details of the elementary reflectors.
  * If jobz = 'N', it is not referenced (can be NULL, ldz ignored).
  * @param[in] ldz       The leading dimension of array Z.
  * If jobz = 'I' or 'F':
  * If row_major=0, ldz >= max(1,n).
  * If row_major=1, ldz >= max(1,n) (number of columns).
  * If jobz = 'N', ldz >= 1 (but ignored).
  * @param[out] tau      Double array, dimension (n).
  * Contains the scalar factors of the elementary reflectors.
  * @param[in] tol       Tolerance for rank determination.
  * If tol > 0, used as lower bound for reciprocal condition number.
  * If tol <= 0, a default tolerance is computed internally.
  * @param[in] row_major Integer flag:
  * = 0: Arrays a, b, z are column-major (Fortran style).
  * = 1: Arrays a, b, z are row-major (C style).
  * (Arrays nblk, tau, iwork are 1D and unaffected).
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * Memory allocation errors may also be returned.
  */
 SLICOT_C_WRAPPER_API
 int slicot_ab01nd(char jobz, int n, int m,
                   double* a, int lda,
                   double* b, int ldb,
                   int* ncont, int* indcon, int* nblk,
                   double* z, int ldz,
                   double* tau, double tol,
                   int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* AB01ND_H */
 