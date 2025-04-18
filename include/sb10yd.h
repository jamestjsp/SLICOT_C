/**
 * @file sb10yd.h
 * @brief C wrapper for SLICOT routine SB10YD
 *
 * This file provides a C interface to the SLICOT routine SB10YD,
 * which fits frequency response data with a stable, minimum phase
 * SISO system.
 */

 #ifndef SB10YD_H
 #define SB10YD_H
 
 #include <stddef.h> // For size_t
 #include "slicot_utils.h" // For slicot_complex_double definition
 
 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Fits frequency response data with a stable, minimum phase SISO system.
  *
  * Fits supplied frequency response data (real/imaginary parts vs. frequency)
  * with a stable, minimum phase SISO system (A, B, C, D). Handles both
  * continuous-time and discrete-time cases. Optionally enforces
  * stability/minimum phase properties.
  *
  * @param[in] discfl    Indicates system type: 0 = continuous, 1 = discrete.
  * @param[in] flag      Constraint flag: 0 = no constraints, 1 = enforce stable/min phase.
  * @param[in] lendat    Length of frequency response data vectors, lendat >= 2.
  * @param[in] rfrdat    Double array, dimension (lendat). Real part of frequency data.
  * @param[in] ifrdat    Double array, dimension (lendat). Imaginary part of frequency data.
  * @param[in] omega     Double array, dimension (lendat). Frequencies (non-negative, increasing).
  * For discrete-time, 0 <= omega[i] <= pi.
  * @param[in,out] n     On entry, desired system order, n <= lendat-1.
  * On exit, order of the obtained system (may change if flag=1).
  * @param[out] a        Double array, dimension (lda, n) or (n, lda). State matrix A.
  * If flag=1, A is upper Hessenberg.
  * @param[in] lda       Leading dimension of A. >= max(1,n).
  * @param[out] b        Double array, dimension (n). Input vector B.
  * @param[out] c        Double array, dimension (n). Output vector C.
  * If flag=1, first n-1 elements are zero.
  * @param[out] d        Double array, dimension (1). Feedthrough scalar D.
  * @param[in] tol       Tolerance for rank determination. If <= 0, default used.
  * @param[in] row_major Integer flag:
  * = 0: Array a is column-major (Fortran style).
  * = 1: Array a is row-major (C style).
  * (Arrays rfrdat, ifrdat, omega, b, c, d are 1D).
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * = 1: Discrete -> continuous transformation failed.
  * = 2: Could not find system poles.
  * = 3: Could not find inverse system (D is near zero).
  * = 4: Could not find system zeros.
  * = 5: Could not find state-space representation of T(s).
  * = 6: Continuous -> discrete transformation failed.
  * Memory allocation errors may also be returned.
  */
 int slicot_sb10yd(int discfl, int flag, int lendat,
                   const double* rfrdat, const double* ifrdat, const double* omega,
                   int* n, double* a, int lda, double* b, double* c, double* d,
                   double tol, int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* SB10YD_H */
 