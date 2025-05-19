/**
 * @file tf01rd.h
 * @brief Header for C wrapper of SLICOT routine TF01RD.
 * @details This routine computes N Markov parameters M(1), M(2),..., M(N)
 * from the parameters (A,B,C) of a linear time-invariant system.
 */

#ifndef TF01RD_H
#define TF01RD_H

#include "slicot_utils.h" // Provides SLICOT_EXPORT macro

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Computes Markov parameters of a multivariable system.
 * @details This is a C wrapper for the SLICOT Fortran routine TF01RD.
 * Workspace (DWORK) is allocated internally by this wrapper.
 *
 * @param na (Input) The order of the state matrix A. na >= 0. (Corresponds to 'n' in Python wrapper)
 * @param nb (Input) The number of system inputs. nb >= 0. (Corresponds to 'm' in Python wrapper)
 * @param nc (Input) The number of system outputs. nc >= 0. (Corresponds to 'p' in Python wrapper)
 * @param N_in (Input) The number of Markov parameters M(k) to be computed. N_in >= 0. (Corresponds to 'N' in Python wrapper)
 * @param a (Input) DOUBLE PRECISION array, state matrix A.
 * Dimensions (na x na). Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param lda (Input) Leading dimension of the C array storing A.
 * If row_major=0, lda >= max(1, na).
 * If row_major=1, lda >= max(1, na) (number of columns).
 * @param b (Input) DOUBLE PRECISION array, input matrix B.
 * Dimensions (na x nb). Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param ldb (Input) Leading dimension of the C array storing B.
 * If row_major=0, ldb >= max(1, na).
 * If row_major=1, ldb >= max(1, nb) (number of columns).
 * @param c (Input) DOUBLE PRECISION array, output matrix C.
 * Dimensions (nc x na). Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param ldc (Input) Leading dimension of the C array storing C.
 * If row_major=0, ldc >= max(1, nc).
 * If row_major=1, ldc >= max(1, na) (number of columns).
 * @param h (Output) DOUBLE PRECISION array, storing the Markov parameters.
 * Dimensions (nc x (N_in*nb)). The k-th Markov parameter M(k) (an nc x nb matrix)
 * is stored in H(:,(k-1)*nb : k*nb -1) using 0-based indexing for columns.
 * Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param ldh (Input) Leading dimension of the C array storing H.
 * If row_major=0, ldh >= max(1, nc).
 * If row_major=1, ldh >= max(1, N_in*nb) (number of columns).
 * @param row_major (Input) Specifies matrix storage for A, B, C, and H:
 * 0 for column-major (Fortran default).
 * 1 for row-major (C default).
 *
 * @return info Error indicator:
 * = 0: successful exit.
 * < 0: if info = -i, the i-th argument had an illegal value.
 * Fortran argument indices:
 * 1:NA, 2:NB, 3:NC, 4:N, 5:A, 6:LDA, 7:B, 8:LDB, 9:C, 10:LDC,
 * 11:H, 12:LDH, 13:DWORK, 14:LDWORK.
 * = SLICOT_MEMORY_ERROR (-1010): internal memory allocation failed.
 */
SLICOT_EXPORT
int slicot_tf01rd(int na, int nb, int nc, int N_in,
                  const double* a, int lda,
                  const double* b, int ldb,
                  const double* c, int ldc,
                  double* h, int ldh,
                  int row_major);

#ifdef __cplusplus
}
#endif

#endif // TF01RD_H
