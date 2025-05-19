/**
 * @file tg01ad.h
 * @brief Header for C wrapper of SLICOT routine TG01AD.
 * @details This routine balances the matrices of a system pencil
 * corresponding to a descriptor triple (A-lambda E,B,C).
 */

#ifndef SLICOT_WRAPPER_TG01AD_H
#define SLICOT_WRAPPER_TG01AD_H

#include "slicot_utils.h" // Provides SLICOT_EXPORT macro

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Balances the matrices of the system pencil S = [A B; C 0] - lambda*[E 0; 0 0].
 * @details This is a C wrapper for the SLICOT Fortran routine TG01AD.
 * Workspace (DWORK) is allocated internally by this wrapper.
 * Matrices A, E, B, C are modified in-place.
 *
 * @param job (Input) CHARACTER*1. Specifies which matrices are involved in balancing:
 * 'A': All matrices A, E, B, C.
 * 'B': Matrices B, A, E.
 * 'C': Matrices C, A, E.
 * 'N': Only matrices A, E (B and C are not involved).
 * @param l (Input) INTEGER. The number of rows of matrices A, B, and E. l >= 0.
 * @param n (Input) INTEGER. The number of columns of matrices A, E, and C. n >= 0.
 * @param m (Input) INTEGER. The number of columns of matrix B. m >= 0.
 * @param p (Input) INTEGER. The number of rows of matrix C. p >= 0.
 * @param thresh (Input) DOUBLE PRECISION. Threshold for balancing. Elements with magnitude <= thresh are ignored. thresh >= 0.
 * @param a (Input/Output) DOUBLE PRECISION array, state dynamics matrix A.
 * Dimensions (l x n). Stored column-wise if row_major=0, row-wise if row_major=1.
 * On exit, contains the balanced matrix Dl*A*Dr.
 * @param lda (Input) INTEGER. Leading dimension of the C array storing A.
 * If row_major=0, lda >= max(1, l).
 * If row_major=1, lda >= max(1, n) (number of columns).
 * @param e (Input/Output) DOUBLE PRECISION array, descriptor matrix E.
 * Dimensions (l x n). Stored column-wise if row_major=0, row-wise if row_major=1.
 * On exit, contains the balanced matrix Dl*E*Dr.
 * @param lde (Input) INTEGER. Leading dimension of the C array storing E.
 * If row_major=0, lde >= max(1, l).
 * If row_major=1, lde >= max(1, n) (number of columns).
 * @param b (Input/Output) DOUBLE PRECISION array, input matrix B.
 * Dimensions (l x m). Stored column-wise if row_major=0, row-wise if row_major=1.
 * Not referenced if m = 0.
 * On exit, if m > 0, contains the balanced matrix Dl*B.
 * @param ldb (Input) INTEGER. Leading dimension of the C array storing B.
 * If m > 0: If row_major=0, ldb >= max(1, l). If row_major=1, ldb >= max(1, m).
 * If m = 0: ldb >= 1.
 * @param c (Input/Output) DOUBLE PRECISION array, output matrix C.
 * Dimensions (p x n). Stored column-wise if row_major=0, row-wise if row_major=1.
 * Not referenced if p = 0.
 * On exit, if p > 0, contains the balanced matrix C*Dr.
 * @param ldc (Input) INTEGER. Leading dimension of the C array storing C.
 * If p > 0: If row_major=0, ldc >= max(1, p). If row_major=1, ldc >= max(1, n).
 * If p = 0: ldc >= 1.
 * @param lscale (Output) DOUBLE PRECISION array, dimension (l).
 * The left scaling factors Dl. lscale[j-1] = Dl(j).
 * @param rscale (Output) DOUBLE PRECISION array, dimension (n).
 * The right scaling factors Dr. rscale[j-1] = Dr(j).
 * @param row_major (Input) INTEGER. Specifies matrix storage for A, E, B, C:
 * 0 for column-major (Fortran default).
 * 1 for row-major (C default).
 *
 * @return info Error indicator:
 * = 0: successful exit.
 * < 0: if info = -i, the i-th argument had an illegal value.
 * Fortran argument indices: 1:JOB, 2:L, 3:N, 4:M, 5:P, 6:THRESH,
 * 7:A, 8:LDA, 9:E, 10:LDE, 11:B, 12:LDB, 13:C, 14:LDC,
 * 15:LSCALE, 16:RSCALE, 17:DWORK.
 * = SLICOT_MEMORY_ERROR (-1010): internal memory allocation failed.
 */
SLICOT_EXPORT
int slicot_tg01ad(char job, int l, int n, int m, int p, double thresh,
                  double* a, int lda,
                  double* e, int lde,
                  double* b, int ldb,
                  double* c, int ldc,
                  double* lscale, double* rscale,
                  int row_major);

#ifdef __cplusplus
}
#endif

#endif // SLICOT_WRAPPER_TG01AD_H
