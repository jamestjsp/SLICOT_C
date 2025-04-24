/**
 * @file ag08bd.h
 * @brief C wrapper for SLICOT routine AG08BD
 *
 * This file provides a C interface to the SLICOT routine AG08BD,
 * which computes the zeros and Kronecker structure of a descriptor
 * system pencil S(lambda) = [A-lambda*E, B; C, D].
 */

#ifndef AG08BD_H
#define AG08BD_H

 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Computes zeros and Kronecker structure of a descriptor system pencil.
 *
 * Extracts from the system pencil S(lambda) = [A-lambda*E, B; C, D]
 * a regular pencil Af-lambda*Ef which has the finite Smith zeros of
 * S(lambda) as generalized eigenvalues. Also computes infinite zero/pole
 * structure and Kronecker indices.
 * The matrices A and E are overwritten with the reduced pencil Af, Ef.
 * The matrices B and C are overwritten during computation.
 *
 * @param[in] equil     Specifies whether to balance the system matrix:
 * = 'S': Perform balancing (scaling).
 * = 'N': Do not perform balancing.
 * @param[in] l         Number of rows of A, E, B, l >= 0.
 * @param[in] n         Number of columns of A, E, C, n >= 0.
 * @param[in] m         Number of columns of B, D, m >= 0.
 * @param[in] p         Number of rows of C, D, p >= 0.
 * @param[in,out] a     Double array, dimension (lda, n) or (l, lda).
 * On entry, the matrix A.
 * On exit, the leading nfz-by-nfz part contains the matrix Af of the reduced pencil.
 * @param[in] lda       The leading dimension of array A. >= max(1,l).
 * @param[in,out] e     Double array, dimension (lde, n) or (l, lde).
 * On entry, the matrix E.
 * On exit, the leading nfz-by-nfz part contains the matrix Ef of the reduced pencil.
 * @param[in] lde       The leading dimension of array E. >= max(1,l).
 * @param[in,out] b     Double array, dimension (ldb, m) or (l, ldb).
 * On entry, the matrix B. On exit, overwritten.
 * @param[in] ldb       The leading dimension of array B. >= max(1,l) if m>0, else >=1.
 * @param[in,out] c     Double array, dimension (ldc, n) or (p, ldc).
 * On entry, the matrix C. On exit, overwritten.
 * @param[in] ldc       The leading dimension of array C. >= max(1,p).
 * @param[in] d         Double array, dimension (ldd, m) or (p, ldd). Matrix D.
 * @param[in] ldd       The leading dimension of array D. >= max(1,p).
 * @param[out] nfz      The number of finite zeros (order of Af, Ef).
 * @param[out] nrank    The normal rank of the system pencil S(lambda).
 * @param[out] niz      The number of infinite zeros.
 * @param[out] dinfz    The maximal multiplicity of infinite Smith zeros.
 * @param[out] nkror    The number of right Kronecker indices.
 * @param[out] ninfe    The number of elementary infinite blocks (poles at infinity).
 * @param[out] nkrol    The number of left Kronecker indices.
 * @param[out] infz     Integer array, dimension (n+1). Infinite zero structure.
 * @param[out] kronr    Integer array, dimension (n+m+1). Right Kronecker indices.
 * @param[out] infe     Integer array, dimension (1+min(l+p,n+m)). Infinite pole structure.
 * @param[out] kronl    Integer array, dimension (l+p+1). Left Kronecker indices.
 * @param[in] tol       Tolerance for rank decisions. If tol<=0, defaults are used. tol < 1.
 * @param[in] row_major Integer flag:
 * = 0: Arrays a, e, b, c, d are column-major (Fortran style).
 * = 1: Arrays a, e, b, c, d are row-major (C style).
 * (Arrays infz, kronr, infe, kronl are 1D).
 *
 * @return info         Error indicator:
 * = 0: successful exit
 * < 0: if info = -i, the i-th argument had an illegal value.
 * Memory allocation errors may also be returned.
 */
SLICOT_EXPORT
int slicot_ag08bd(char equil, int l, int n, int m, int p,
                  double* a, int lda, double* e, int lde,
                  double* b, int ldb, double* c, int ldc,
                  const double* d, int ldd,
                  int* nfz, int* nrank, int* niz, int* dinfz,
                  int* nkror, int* ninfe, int* nkrol,
                  int* infz, int* kronr, int* infe, int* kronl,
                  double tol, int row_major);

#ifdef __cplusplus
}
#endif

#endif /* AG08BD_H */
