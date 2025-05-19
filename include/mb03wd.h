/**
 * @file mb03wd.h
 * @brief C wrapper for SLICOT routine MB03WD
 *
 * This file provides a C interface to the SLICOT routine MB03WD,
 * which computes the Schur decomposition and eigenvalues of a product
 * of matrices H = H_1*...*H_p in periodic Hessenberg form.
 */

#ifndef MB03WD_H
#define MB03WD_H

 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Computes Schur decomposition of a product of matrices in periodic Hessenberg form.
 *
 * Computes the Schur decomposition and eigenvalues of H = H_1*...*H_p,
 * where H_1 is upper Hessenberg and H_2..H_p are upper triangular.
 * Computes orthogonal Z_i such that Z_1'*H_1*Z_2=T_1, ..., Z_p'*H_p*Z_1=T_p,
 * where T_1 is real Schur form and T_2..T_p are upper triangular.
 *
 * @param[in] job       Specifies computation level:
 * = 'E': Compute eigenvalues only.
 * = 'S': Compute Schur factors T_1, ..., T_p (returned in H).
 * @param[in] compz     Specifies computation of Schur vectors Z_i:
 * = 'N': Do not compute Z_i.
 * = 'I': Initialize Z_i to identity and return computed Z_i.
 * = 'V': Update Z_i by post-multiplying with computed transformations (Z_i must be input).
 * @param[in] n         The order of the matrices H_j, n >= 0.
 * @param[in] p         The number of matrices in the product, p >= 1.
 * @param[in] ilo       Lower row/column index defining the subproblem.
 * @param[in] ihi       Upper row/column index defining the subproblem.
 * 1 <= ilo <= max(1,n); min(ilo,n) <= ihi <= n.
 * @param[in] iloz      Lower row index for applying transformations to Z.
 * @param[in] ihiz      Upper row index for applying transformations to Z.
 * 1 <= iloz <= ilo; ihi <= ihiz <= n.
 * @param[in,out] h     Double array, dimension (ldh1, ldh2, p) stored contiguously.
 * On entry, contains H_1 (Hessenberg) and H_2..H_p (triangular).
 * On exit (if job='S'), contains Schur factors T_1..T_p.
 * If job='E', contents are unspecified on exit.
 * @param[in] ldh1      The first dimension of H:
 *                      If row_major=0, this is the leading dimension (stride between rows). >= max(1,n).
 *                      If row_major=1, this is the number of rows in each 2D slice. >= n.
 * @param[in] ldh2      The second dimension of H:
 *                      If row_major=0, this is the second dimension (stride between columns). >= max(1,n).
 *                      If row_major=1, this is the leading dimension (stride between columns). >= n.
 * @param[in,out] z     Double array, dimension (ldz1, ldz2, p) stored contiguously.
 * If compz='V', contains input orthogonal matrices on entry.
 * On exit (if compz='I' or 'V'), contains the computed/updated
 * Schur vectors Z_1 to Z_p. Not referenced if compz='N'.
 * @param[in] ldz1      The first dimension of Z:
 *                      If row_major=0, this is the leading dimension (stride between rows).
 *                      >=1 if compz='N', >=max(1,n) otherwise.
 *                      If row_major=1, this is the number of rows in each 2D slice.
 *                      >=1 if compz='N', >=n otherwise.
 * @param[in] ldz2      The second dimension of Z:
 *                      If row_major=0, this is the second dimension (stride between columns).
 *                      >=1 if compz='N', >=max(1,n) otherwise.
 *                      If row_major=1, this is the leading dimension (stride between columns).
 *                      >=1 if compz='N', >=n otherwise.
 * @param[out] wr       Double array, dimension (n). Real parts of eigenvalues.
 * @param[out] wi       Double array, dimension (n). Imaginary parts of eigenvalues.
 * @param[in] row_major Specifies matrix storage format:
 *                      0 for column-major (Fortran-style) layout.
 *                      1 for row-major (C-style) layout.
 *
 * @return info         Error indicator:
 * = 0: successful exit
 * < 0: if info = -i, the i-th argument had an illegal value.
 * > 0: if info = i, QR algorithm failed to compute eigenvalues i+1:ihi.
 * = SLICOT_MEMORY_ERROR (-1010): internal memory allocation failed.
 */
SLICOT_EXPORT
int slicot_mb03wd(char job, char compz, int n, int p, int ilo, int ihi,
                  int iloz, int ihiz,
                  double* h, int ldh1, int ldh2,
                  double* z, int ldz1, int ldz2,
                  double* wr, double* wi,
                  int row_major);

#ifdef __cplusplus
}
#endif

#endif /* MB03WD_H */
