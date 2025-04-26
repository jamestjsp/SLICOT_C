/**  
 * @file dg01nd.h  
 * @brief Header for C wrapper of SLICOT routine DG01ND.  
 */

#ifndef SLICOT_WRAPPER_DG01ND_H
#define SLICOT_WRAPPER_DG01ND_H

#include "slicot_utils.h" // Provides SLICOT_EXPORT macro

#ifdef __cplusplus  
extern "C" {  
#endif

/**  
 * @brief Computes the discrete Fourier transform, or inverse Fourier transform, of a real signal.
 * @details This is a C wrapper for the SLICOT Fortran routine DG01ND.
 * **Workspace is allocated internally.**
 *
 * @param indi Indicates whether a Fourier transform or inverse Fourier transform is to be performed:
 *             'D': (Direct) Fourier transform;
 *             'I': Inverse Fourier transform.
 * @param n Half the number of real samples. N must be a power of 2. N >= 2.
 * @param xr [in/out] For INDI = 'D', the first N elements contain the odd part of the input signal;
 *           for INDI = 'I', the first N+1 elements contain the real part of the input discrete
 *           Fourier transform.
 *           On exit with INDI = 'D', the first N+1 elements contain the real part of the
 *           output signal (computed discrete Fourier transform).
 *           On exit with INDI = 'I', the first N elements contain the odd part of the
 *           output signal (computed inverse discrete Fourier transform).
 * @param xi [in/out] For INDI = 'D', the first N elements contain the even part of the input signal;
 *           for INDI = 'I', the first N+1 elements contain the imaginary part of the input
 *           discrete Fourier transform.
 *           On exit with INDI = 'D', the first N+1 elements contain the imaginary part of the
 *           output signal (computed discrete Fourier transform).
 *           On exit with INDI = 'I', the first N elements contain the even part of the
 *           output signal (computed inverse discrete Fourier transform).
 * @param row_major Specifies data storage format for the input/output arrays:
 *                 0 for column-major (Fortran default),
 *                 1 for row-major (C default).
 *                 Note: For this function, the parameter affects memory handling but not
 *                 the array layout since we're working with vectors, not matrices.
 *  
 * @return info Error indicator:  
 * = 0: successful exit  
 * < 0: if info = -i, the i-th argument had an illegal value (wrapper or Fortran validation)
 * = SLICOT_MEMORY_ERROR (-1010): internal memory allocation failed.  
 */  
SLICOT_EXPORT  
int slicot_dg01nd(char indi, int n, double* xr, double* xi, int row_major);

#ifdef __cplusplus  
}  
#endif

#endif /* SLICOT_WRAPPER_DG01ND_H */
