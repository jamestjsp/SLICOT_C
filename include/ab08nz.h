/**
 * @file ab08nz.h
 * @brief C wrapper for SLICOT routine AB08NZ
 *
 * This file provides a C interface to the SLICOT routine AB08NZ,
 * which constructs a regular pencil for a given complex system such that its
 * generalized eigenvalues are invariant zeros of the system.
 */

 #ifndef AB08NZ_H
 #define AB08NZ_H
 
 #include <stddef.h> // For size_t
 #include "slicot_utils.h" // For slicot_complex_double definition
 
 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Constructs a regular pencil for a complex system to compute invariant zeros.
  *
  * This function constructs a regular pencil (Af - lambda*Bf) for a linear
  * multivariable complex system described by a state-space model (A,B,C,D) which
  * has the invariant zeros of the system as generalized eigenvalues. It also computes
  * the orders of the infinite zeros and the right and left Kronecker indices.
  * Input matrices might be scaled in place if balancing is requested.
  *
  * @param[in] equil     Specifies whether to balance the compound matrix:
  * = 'S': Perform balancing (scaling). Input matrices a,b,c,d might be modified.
  * = 'N': Do not perform balancing. Input matrices are read-only.
  * @param[in] n         The order of matrix A (number of state variables), n >= 0.
  * @param[in] m         The number of system inputs (columns of B, D), m >= 0.
  * @param[in] p         The number of system outputs (rows of C, D), p >= 0.
  * @param[in,out] a     Complex array for matrix A. Dimension (lda, n) or (n, lda).
  * On entry, the state matrix A. On exit, if equil='S', potentially scaled A.
  * @param[in] lda       The leading dimension of array A. >= max(1,n).
  * @param[in,out] b     Complex array for matrix B. Dimension (ldb, m) or (n, ldb).
  * On entry, the input matrix B. On exit, if equil='S', potentially scaled B.
  * @param[in] ldb       The leading dimension of array B. >= max(1,n).
  * @param[in,out] c     Complex array for matrix C. Dimension (ldc, n) or (p, ldc).
  * On entry, the output matrix C. On exit, if equil='S', potentially scaled C.
  * @param[in] ldc       The leading dimension of array C. >= max(1,p).
  * @param[in,out] d     Complex array for matrix D. Dimension (ldd, m) or (p, ldd).
  * On entry, the feedthrough matrix D. On exit, if equil='S', potentially scaled D.
  * @param[in] ldd       The leading dimension of array D. >= max(1,p).
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
  * @param[out] af       Complex array, dimension (ldaf, *) or (*, ldaf).
  * Fortran dimension is (LDAF, N+MIN(P,M)).
  * The leading nu-by-nu part contains the coefficient matrix Af of the reduced pencil.
  * @param[in] ldaf      The leading dimension of array AF. >= max(1,n+m).
  * @param[out] bf       Complex array, dimension (ldbf, *) or (*, ldbf).
  * Fortran dimension is (LDBF, N+M).
  * The leading nu-by-nu part contains the coefficient matrix Bf of the reduced pencil.
  * @param[in] ldbf      The leading dimension of array BF. >= max(1,n+p).
  * @param[in] tol       Tolerance used for rank decisions. If tol <= 0, a default is used.
  * @param[in] row_major Integer flag:
  * = 0: Arrays a, b, c, d, af, bf are column-major (Fortran style).
  * = 1: Arrays a, b, c, d, af, bf are row-major (C style).
  * (Arrays infz, kronr, kronl are 1D integer arrays).
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * Memory allocation errors may also be returned.
  * Special codes -20, -22 may be returned if ldaf/ldbf are too small for computed nu.
  */
 SLICOT_C_WRAPPER_API
 int slicot_ab08nz(char equil, int n, int m, int p,
                   slicot_complex_double* a, int lda,
                   slicot_complex_double* b, int ldb,
                   slicot_complex_double* c, int ldc,
                   slicot_complex_double* d, int ldd,
                   int* nu, int* rank,
                   int* dinfz, int* nkror, int* nkrol,
                   int* infz, int* kronr, int* kronl,
                   slicot_complex_double* af, int ldaf,
                   slicot_complex_double* bf, int ldbf,
                   double tol, int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* AB08NZ_H */
 