/**
 * @file ab08nd.h
 * @brief C wrapper for SLICOT routine AB08ND
 *
 * This file provides a C interface to the SLICOT routine AB08ND,
 * which constructs a regular pencil for a given system such that its
 * generalized eigenvalues are invariant zeros of the system.
 */

 #ifndef AB08ND_H
 #define AB08ND_H
 
 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Constructs a regular pencil for a system to compute invariant zeros
  *
  * This function constructs a regular pencil (A_f - lambda*B_f) for a linear
  * multivariable system described by a state-space model (A,B,C,D) which has
  * the invariant zeros of the system as generalized eigenvalues. It also computes
  * the orders of the infinite zeros and the right and left Kronecker indices.
  * Input matrices might be scaled in place if balancing is requested.
  *
  * @param[in] equil     Specifies whether to balance the compound matrix:
  * - 'S': Perform balancing (scaling). Input matrices a,b,c,d might be modified.
  * - 'N': Do not perform balancing. Input matrices are read-only.
  * @param[in] n         The order of matrix A (number of state variables), n >= 0.
  * @param[in] m         The number of system inputs, m >= 0.
  * @param[in] p         The number of system outputs, p >= 0.
  * @param[in,out] a     Double array, dimension (lda, n) (col-major) or (n, lda) (row-major).
  * On entry, the state dynamics matrix A.
  * On exit, if equil='S', it may be overwritten by the scaled matrix.
  * @param[in] lda       The leading dimension of array A.
  * If row_major=0, lda >= max(1,n).
  * If row_major=1, lda >= max(1,n) (number of columns).
  * @param[in,out] b     Double array, dimension (ldb, m) (col-major) or (n, ldb) (row-major).
  * On entry, the input/state matrix B.
  * On exit, if equil='S', it may be overwritten by the scaled matrix.
  * @param[in] ldb       The leading dimension of array B.
  * If row_major=0, ldb >= max(1,n).
  * If row_major=1, ldb >= max(1,m) (number of columns).
  * @param[in,out] c     Double array, dimension (ldc, n) (col-major) or (p, ldc) (row-major).
  * On entry, the state/output matrix C.
  * On exit, if equil='S', it may be overwritten by the scaled matrix.
  * @param[in] ldc       The leading dimension of array C.
  * If row_major=0, ldc >= max(1,p).
  * If row_major=1, ldc >= max(1,n) (number of columns).
  * @param[in,out] d     Double array, dimension (ldd, m) (col-major) or (p, ldd) (row-major).
  * On entry, the direct transmission matrix D.
  * On exit, if equil='S', it may be overwritten by the scaled matrix.
  * @param[in] ldd       The leading dimension of array D.
  * If row_major=0, ldd >= max(1,p).
  * If row_major=1, ldd >= max(1,m) (number of columns).
  * @param[out] nu       The number of (finite) invariant zeros.
  * @param[out] rank     The normal rank of the transfer function matrix.
  * @param[out] dinfz    The maximum degree of infinite elementary divisors.
  * @param[out] nkror    The number of right Kronecker indices.
  * @param[out] nkrol    The number of left Kronecker indices.
  * @param[out] infz     Integer array, dimension (n).
  * The leading dinfz elements contain information on the
  * infinite elementary divisors: infz[i-1] contains the number
  * of infinite elementary divisors of degree i for i=1,2,...,dinfz.
  * @param[out] kronr    Integer array, dimension (max(n,m)+1).
  * The leading nkror elements contain the right Kronecker indices.
  * @param[out] kronl    Integer array, dimension (max(n,p)+1).
  * The leading nkrol elements contain the left Kronecker indices.
  * @param[out] af       Double array, dimension (ldaf, *) (col-major) or (*, ldaf) (row-major).
  * Fortran dimension is (LDAF, N+MIN(P,M)).
  * The leading nu-by-nu part contains the coefficient matrix Af of
  * the reduced pencil.
  * @param[in] ldaf      The leading dimension of array AF.
  * If row_major=0, ldaf >= max(1,n+m).
  * If row_major=1, ldaf >= max(1,nu) (number of columns).
  * @param[out] bf       Double array, dimension (ldbf, *) (col-major) or (*, ldbf) (row-major).
  * Fortran dimension is (LDBF, N+M).
  * The leading nu-by-nu part contains the coefficient matrix Bf of
  * the reduced pencil.
  * @param[in] ldbf      The leading dimension of array BF.
  * If row_major=0, ldbf >= max(1,n+p).
  * If row_major=1, ldbf >= max(1,nu) (number of columns).
  * @param[in] tol       Tolerance used for rank decisions. If tol is less than
  * sqrt((n+p)*(n+m))*eps, a default value based on machine
  * precision is used internally by the Fortran routine.
  * @param[in] row_major Integer flag:
  * - 0: Arrays a, b, c, d, af, bf are column-major (Fortran style).
  * - 1: Arrays a, b, c, d, af, bf are row-major (C style).
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * Memory allocation errors may also be returned.
  * Special codes -20, -22 may be returned if ldaf/ldbf are too small for computed nu.
  */
 SLICOT_C_WRAPPER_API
 int slicot_ab08nd(char equil, int n, int m, int p,
                   double* a, int lda,
                   double* b, int ldb,
                   double* c, int ldc,
                   double* d, int ldd,
                   int* nu, int* rank,
                   int* dinfz, int* nkror, int* nkrol,
                   int* infz, int* kronr, int* kronl,
                   double* af, int ldaf,
                   double* bf, int ldbf,
                   double tol, int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* AB08ND_H */
 