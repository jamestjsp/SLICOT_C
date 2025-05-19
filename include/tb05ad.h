/**  
 * @file tb05ad.h  
 * @brief Header for C wrapper of SLICOT routine TB05AD.  
 */

#ifndef TB05AD_H
#define TB05AD_H

#include "slicot_utils.h" // Provides SLICOT_EXPORT macro

#ifdef __cplusplus  
extern "C" {  
#endif

/**  
 * @brief Frequency response matrix of a given state-space representation (A,B,C).
 * @details Finds the complex frequency response matrix (transfer matrix)
 * G(freq) of the state-space representation (A,B,C) given by:
 *                              -1
 *     G(freq) = C * ((freq*I - A)  ) * B
 *
 * where A, B and C are real N-by-N, N-by-M and P-by-N matrices
 * respectively and freq is a complex scalar.
 *
 * This is a C wrapper for the SLICOT Fortran routine TB05AD.  
 * **Workspace is allocated internally.**  
 *  
 * @param baleig Determines whether to balance matrix A or compute eigenvalues:
 *               'N': No balancing or eigenvalue calculation
 *               'C': No balancing but calculate condition number
 *               'B','E': Balance A and calculate eigenvalues
 *               'A': Balance A, calculate eigenvalues and condition number
 * @param inita  Specifies if A is already in upper Hessenberg form:
 *               'G': General matrix
 *               'H': Upper Hessenberg form (no balancing/eigenvalues needed)
 * @param n      Number of states, i.e., order of matrix A. n >= 0. 
 * @param m      Number of inputs, i.e., columns in matrix B. m >= 0.
 * @param p      Number of outputs, i.e., rows in matrix C. p >= 0.
 * @param freq_re Real part of the frequency freq at which to evaluate response.
 * @param freq_im Imaginary part of the frequency freq at which to evaluate response.
 * @param a      [in,out] State transition matrix A. Dimensions (n x n).
 *               On exit, contains upper Hessenberg matrix (if inita = 'G').
 *               Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param lda    Leading dimension of the C array storing A.
 *               If row_major=0 (column-major), lda >= max(1, n) (number of rows).
 *               If row_major=1 (row-major), lda >= max(1, n) (number of columns).
 * @param b      [in,out] Input/state matrix B. Dimensions (n x m).
 *               On exit, transformed if inita = 'G'.
 *               Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param ldb    Leading dimension of the C array storing B.
 *               If row_major=0 (column-major), ldb >= max(1, n) (number of rows).
 *               If row_major=1 (row-major), ldb >= max(1, m) (number of columns).
 * @param c      [in,out] State/output matrix C. Dimensions (p x n).
 *               On exit, transformed if inita = 'G'.
 *               Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param ldc    Leading dimension of the C array storing C.
 *               If row_major=0 (column-major), ldc >= max(1, p) (number of rows).
 *               If row_major=1 (row-major), ldc >= max(1, n) (number of columns).
 * @param rcond  [out] Estimate of the reciprocal condition number (if baleig='C' or 'A').
 * @param g      [out] Frequency response matrix G(freq). Dimensions (p x m).
 *               Complex values are stored interleaved (re0,im0,re1,im1,...).
 *               Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param ldg    Leading dimension of the C array storing G.
 *               If row_major=0 (column-major), ldg >= max(1, p) (number of rows).
 *               If row_major=1 (row-major), ldg >= max(1, m) (number of columns).
 * @param evre   [out] Real parts of eigenvalues of A (if applicable).
 *               Must have size >= n.
 * @param evim   [out] Imaginary parts of eigenvalues of A (if applicable).
 *               Must have size >= n.
 * @param hinvb  [out] Product H^-1 * B. Dimensions (n x m).
 *               Complex values are stored interleaved (re0,im0,re1,im1,...).
 *               Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param ldhinv Leading dimension of the C array storing HINVB.
 *               If row_major=0 (column-major), ldhinv >= max(1, n) (number of rows).
 *               If row_major=1 (row-major), ldhinv >= max(1, m) (number of columns).
 * @param row_major Specifies storage layout for input/output matrices:
 *               0 for column-major (Fortran default),
 *               1 for row-major (C default).
 *  
 * @return info Error indicator:
 *         = 0: successful exit
 *         < 0: if info = -i, the i-th argument had an illegal value
 *         = 1: if more than 30*N iterations required to isolate eigenvalues
 *         = 2: if either FREQ is too near to an eigenvalue of A, or RCOND < EPS
 *         = SLICOT_MEMORY_ERROR (-1010): internal memory allocation failed.
 */  
SLICOT_EXPORT  
int slicot_tb05ad(char baleig, char inita, int n, int m, int p,
                 double freq_re, double freq_im,
                 double* a, int lda, double* b, int ldb, double* c, int ldc,
                 double* rcond, double* g, int ldg,
                 double* evre, double* evim,
                 double* hinvb, int ldhinv, int row_major);

#ifdef __cplusplus  
}  
#endif

#endif /* TB05AD_H */
