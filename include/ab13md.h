/**
 * @file ab13md.h
 * @brief C wrapper for SLICOT routine AB13MD
 *
 * This file provides a C interface to the SLICOT routine AB13MD,
 * which computes an upper bound on the structured singular value for a
 * square complex matrix Z with a given block uncertainty structure.
 */

#ifndef AB13MD_H
#define AB13MD_H

#include "slicot_utils.h" // Include definition for slicot_complex_double

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Computes an upper bound on the structured singular value (mu).
 *
 * Computes an upper bound for the structured singular value of a complex
 * matrix Z, given a block structure for the uncertainty.
 *
 * @param[in] fact      Specifies if X contains information from previous call:
 * = 'F': Use information in X.
 * = 'N': Do not use information in X.
 * @param[in] n         The order of the complex matrix Z, n >= 0.
 * @param[in] z         Complex array for matrix Z. Dimension (ldz, n) or (n, ldz).
 * @param[in] ldz       The leading dimension of array Z. >= max(1,n).
 * @param[in] m         The number of diagonal uncertainty blocks, m >= 1.
 * @param[in] nblock    Integer array, dimension (m). Block sizes (sum must equal n).
 * @param[in] itype     Integer array, dimension (m). Block types (1=real, 2=complex).
 * Real blocks must have size 1 (nblock[i]=1 if itype[i]=1).
 * @param[in,out] x     Double array, dimension (m + number of real blocks - 1).
 * If fact='F' and nblock[0]<n, contains info from previous call on entry.
 * On exit, if nblock[0]<n, contains info for next call. Not used if nblock[0]=n.
 * @param[out] bound    The computed upper bound for the structured singular value.
 * @param[out] d        Double array, dimension (n). Diagonal elements of scaling matrix D.
 * @param[out] g        Double array, dimension (n). Diagonal elements of scaling matrix G.
 * @param[in] row_major Integer flag:
 * = 0: Array z is column-major (Fortran style).
 * = 1: Array z is row-major (C style).
 * (Arrays nblock, itype, x, d, g are 1D).
 *
 * @return info         Error indicator:
 * = 0: successful exit
 * < 0: if info = -i, the i-th argument had an illegal value.
 * = 1: Block sizes must be positive.
 * = 2: Sum of block sizes must equal N.
 * = 3: Real block size must be 1.
 * = 4: Block type must be 1 or 2.
 * = 5: Error solving linear equations or inverting matrix.
 * = 6: Error computing eigenvalues or singular values.
 * Memory allocation errors may also be returned.
 */
int slicot_ab13md(char fact, int n, const slicot_complex_double* z, int ldz,
                  int m, const int* nblock, const int* itype,
                  double* x, double* bound, double* d, double* g,
                  int row_major);

#ifdef __cplusplus
}
#endif

#endif /* AB13MD_H */
