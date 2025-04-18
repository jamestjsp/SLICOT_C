/**
 * @file ab01md.h
 * @brief C wrapper for SLICOT routine AB01MD
 *
 * This file provides a C interface to the SLICOT routine AB01MD,
 * which finds a controllable realization for a linear time-invariant
 * single-input system using orthogonal transformations.
 */

 #ifndef AB01MD_H
 #define AB01MD_H
 
 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Finds a controllable realization for a single-input system (A,B).
  *
  * Reduces the system (A,B) using orthogonal similarity transformations Z
  * such that Z'*A*Z is upper Hessenberg and Z'*B = [beta, 0, ..., 0]^T.
  * Determines the order NCONT of the controllable subsystem.
  *
  * @param[in] jobz      Specifies whether to compute the transformation matrix Z:
  * = 'N': Do not compute Z.
  * = 'F': Compute Z in factored form (reflectors stored in Z and TAU).
  * = 'I': Compute Z explicitly (initialize Z to identity before call).
  * @param[in] n         The order of the matrix A, n >= 0.
  * @param[in,out] a     Double array, dimension (lda, n) (col-major) or (n, lda) (row-major).
  * On entry, the state dynamics matrix A.
  * On exit, the leading NCONT-by-NCONT upper Hessenberg part
  * contains the transformed matrix Z'*A*Z. Elements below the
  * first subdiagonal are zeroed.
  * @param[in] lda       The leading dimension of array A.
  * If row_major=0, lda >= max(1,n).
  * If row_major=1, lda >= max(1,n) (number of columns).
  * @param[in,out] b     Double array, dimension (n).
  * On entry, the input/state vector B.
  * On exit, the leading NCONT elements contain the transformed
  * vector Z'*B, with only the first element potentially non-zero.
  * @param[out] ncont    The order of the controllable state-space representation.
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
  * @param[in] tol       Tolerance for determining controllability.
  * If tol > 0, used as an absolute tolerance.
  * If tol <= 0, a default tolerance is computed internally.
  * @param[in] row_major Integer flag:
  * = 0: Arrays a, z are column-major (Fortran style).
  * = 1: Arrays a, z are row-major (C style).
  * (Arrays b and tau are 1D and unaffected).
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * Memory allocation errors may also be returned.
  */
 int slicot_ab01md(char jobz, int n,
                   double* a, int lda,
                   double* b,
                   int* ncont,
                   double* z, int ldz,
                   double* tau, double tol, int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* AB01MD_H */
 