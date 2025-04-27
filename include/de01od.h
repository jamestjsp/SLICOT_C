/**
 * @file de01od.h
 * @brief Header for C wrapper of SLICOT routine DE01OD.
 */

#ifndef SLICOT_WRAPPER_DE01OD_H
#define SLICOT_WRAPPER_DE01OD_H

#include "slicot_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Computes the convolution or deconvolution of two real signals.
 * @details This is a C wrapper for the SLICOT Fortran routine DE01OD.
 * The algorithm uses an FFT approach. No workspace is required.
 *
 * @param conv Specifies whether convolution or deconvolution is performed:
 *             'C' or 'c' for convolution
 *             'D' or 'd' for deconvolution
 * @param n [in] Number of samples. Must be a power of 2, n >= 2.
 * @param a [in,out] Array of dimension n containing the first signal.
 *        On exit, contains the convolution or deconvolution result.
 *        Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param b [in] Array of dimension n containing the second signal.
 *        NOTE: This array is overwritten during computation.
 *        Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param row_major Specifies matrix storage for input/output arrays:
 *        0 for column-major (Fortran default),
 *        1 for row-major (C default).
 *
 * @return info Error indicator:
 *         = 0: successful exit
 *         < 0: if info = -i, the i-th argument had an illegal value
 *         = SLICOT_MEMORY_ERROR (-1010): internal memory allocation failed.
 */
SLICOT_EXPORT
int slicot_de01od(char conv, int n, double* a, double* b, int row_major);

#ifdef __cplusplus
}
#endif

#endif /* SLICOT_WRAPPER_DE01OD_H */
