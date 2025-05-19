/**  
 * @file dk01md.h  
 * @brief Header for C wrapper of SLICOT routine DK01MD.  
 */

#ifndef DK01MD_H
#define DK01MD_H

#include "slicot_utils.h" // Provides SLICOT_EXPORT macro

#ifdef __cplusplus  
extern "C" {  
#endif

/**  
 * @brief Apply an anti-aliasing window to a real signal.
 * @details This is a C wrapper for the SLICOT Fortran routine DK01MD.
 * **Workspace is allocated internally.**
 *
 * The routine applies one of three window types to a signal:
 * - Hamming window: A(i) = (0.54 + 0.46*cos(pi*(i-1)/(N-1)))*A(i)
 * - Hann window: A(i) = 0.5*(1 + cos(pi*(i-1)/(N-1)))*A(i)
 * - Quadratic window: (see documentation for formula)
 *
 * @param type Indicates the type of window to be applied to the signal:
 *             'M': Hamming window
 *             'N': Hann window
 *             'Q': Quadratic window
 * @param n The number of samples in vector A. N >= 1.
 * @param a [in,out] On entry, the vector containing the signal to be processed.
 *          On exit, the windowing function applied to the signal. Dimension N.
 * @param row_major Specifies data storage format:
 *                 0 for column-major (Fortran default),
 *                 1 for row-major (C default).
 *                 For vectors, this has minimal impact but affects memory handling.
 *  
 * @return info Error indicator:  
 * = 0: successful exit  
 * < 0: if info = -i, the i-th argument had an illegal value (wrapper or Fortran validation)
 * = SLICOT_MEMORY_ERROR (-1010): internal memory allocation failed.  
 */  
SLICOT_EXPORT  
int slicot_dk01md(char type, int n, double* a, int row_major);

#ifdef __cplusplus  
}  
#endif

#endif /* DK01MD_H */
