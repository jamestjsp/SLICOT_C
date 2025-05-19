/**
 * @file sb10yd.h
 * @brief Header for C wrapper of SLICOT routine SB10YD.
 */

#ifndef SB10YD_H
#define SB10YD_H

#include "slicot_utils.h" // Provides SLICOT_EXPORT macro and slicot_complex_double

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Fits supplied frequency response data with a stable, minimum phase SISO system.
 * @details This is a C wrapper for the SLICOT Fortran routine SB10YD.
 * Workspace is allocated internally.
 *
 * @param discfl (input) System type: 0 for continuous-time, 1 for discrete-time.
 * @param flag   (input) Constraint flag: 0 for no constraints, 1 for stable/min-phase.
 * @param lendat (input) Length of RFRDAT, IFRDAT, OMEGA. LENDAT >= 2.
 * @param rfrdat (input) Real part of frequency data. Dimension (LENDAT).
 * @param ifrdat (input) Imaginary part of frequency data. Dimension (LENDAT).
 * @param omega  (input) Frequencies. Dimension (LENDAT). Nonnegative, increasing. For discrete, 0 to PI.
 * @param n_io   (input/output) On entry, desired system order N. N <= LENDAT-1.
 * On exit, actual order of the obtained system.
 * @param a      (output) Matrix A of the identified system. Dimension (LDA, N_out).
 * If row_major=0, LDA >= max(1,N_out). If row_major=1, LDA >= max(1,N_out) (cols).
 * @param lda    (input) Leading dimension of A.
 * @param b      (output) Vector B of the identified system. Dimension (N_out).
 * @param c      (output) Vector C of the identified system. Dimension (N_out).
 * @param d      (output) Scalar D of the identified system. Dimension (1).
 * @param tol_param (input) Tolerance for rank determination. If <= 0, a default is used.
 * @param row_major (input) Specifies matrix storage for A: 0 for column-major, 1 for row-major.
 *
 * @return info Error indicator:
 * = 0:  successful exit
 * < 0:  if info = -i, the i-th argument had an illegal value
 * = 1:  discrete to continuous transformation failed
 * = 2:  system poles cannot be found
 * = 3:  inverse system cannot be found (D is near zero)
 * = 4:  system zeros cannot be found
 * = 5:  state-space representation of new T(s) cannot be found
 * = 6:  continuous to discrete transformation failed
 * = SLICOT_MEMORY_ERROR (-1010): internal memory allocation failed.
 */
SLICOT_EXPORT
int slicot_sb10yd(int discfl, int flag, int lendat,
                  const double* rfrdat, const double* ifrdat, const double* omega,
                  int* n_io,
                  double* a, int lda,
                  double* b, double* c, double* d,
                  double tol_param,
                  int row_major);

#ifdef __cplusplus
}
#endif

#endif /* SB10YD_H */
