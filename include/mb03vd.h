/**
 * @file mb03vd.h
 * @brief C wrapper for SLICOT routine MB03VD
 *
 * This file provides a C interface to the SLICOT routine MB03VD,
 * which reduces a product of p real general matrices A = A_1*...*A_p
 * to periodic Hessenberg form H = H_1*...*H_p using orthogonal
 * similarity transformations.
 */

#ifndef MB03VD_H
#define MB03VD_H

 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Reduces a product of matrices to periodic Hessenberg form.
 *
 * Reduces A = A_1*...*A_p to H = H_1*...*H_p, where H_1 is upper
 * Hessenberg and H_2..H_p are upper triangular, using orthogonal
 * similarity transformations Q_j such that Q_j'*A_j*Q_{j+1} = H_j
 * (with Q_{p+1}=Q_1). The elementary reflectors defining Q_j are
 * stored implicitly in the output A and TAU arrays.
 *
 * @param[in] n         The order of the square matrices A_j, n >= 0.
 * @param[in] p         The number of matrices in the product, p >= 1.
 * @param[in] ilo       Lower row/column index for the submatrices to be processed.
 * @param[in] ihi       Upper row/column index for the submatrices to be processed.
 * 1 <= ilo <= max(1,n); min(ilo,n) <= ihi <= n.
 * @param[in,out] a     Double array, dimension (lda1, lda2, p) stored contiguously
 * in **column-major (Fortran-style)** order.
 * On entry, contains the input matrices A_1 to A_p.
 * On exit, A(*,*,1) contains H_1 in its upper triangle and
 * first subdiagonal, and reflector data below. A(*,*,j) for j>1
 * contains H_j in its upper triangle and reflector data below.
 * @param[in] lda1      The first leading dimension of A (stride between rows). >= max(1,n).
 * @param[in] lda2      The second leading dimension of A (stride between columns). >= max(1,n).
 * @param[out] tau      Double array, dimension (ldtau, p). Contains scalar factors
 * of the elementary reflectors for each Q_j.
 * @param[in] ldtau     The leading dimension of TAU. >= max(1,n-1).
 * @param[in] row_major This flag is ignored for the 3D array 'a'. The data for 'a'
 * **must** be provided in column-major (Fortran-style) layout.
 *
 * @return info         Error indicator:
 * = 0: successful exit
 * < 0: if info = -i, the i-th argument had an illegal value.
 * Memory allocation errors may also be returned by the wrapper (though unlikely here).
 */
SLICOT_C_WRAPPER_API
int slicot_mb03vd(int n, int p, int ilo, int ihi,
                  double* a, int lda1, int lda2,
                  double* tau, int ldtau,
                  int row_major); // row_major ignored for 'a'

#ifdef __cplusplus
}
#endif

#endif /* MB03VD_H */
