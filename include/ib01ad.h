/**
 * @file ib01ad.h
 * @brief Header for C wrapper of SLICOT routine IB01AD
 *
 * This header provides the interface for the C wrapper of the
 * SLICOT routine IB01AD, which estimates the system order and
 * computes the triangular factor of the concatenated block-Hankel
 * matrices for subsequent system identification.
 *
 * @note This version corresponds to the wrapper that handles workspace
 * allocation internally.
 */

 #ifndef IB01AD_H
 #define IB01AD_H
 
 #include "slicot_utils.h" // For SLICOT_EXPORT macro
 
 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Estimates system order and computes the triangular factor of block-Hankel matrices
  *
  * This routine implements a combined MOESP and N4SID subspace identification algorithm.
  * It estimates the system order and computes the triangular factor of the concatenated
  * block-Hankel matrices for later use in system identification.
  * Workspace arrays (IWORK, DWORK) are allocated internally by the wrapper.
  *
  * @param meth      Method to use: 'M' for MOESP, 'N' for N4SID algorithm
  * @param alg       Algorithm: 'C' for Cholesky, 'F' for Fast QR, 'Q' for QR algorithm
  * @param jobd      Specifies whether or not the matrices B and D should be computed using MOESP
  * approach: 'M' to later use MOESP, 'N' don't use MOESP
  * @param batch     Data processing approach: 'F' First batch, 'I' Intermediate batch,
  * 'L' Last batch, 'O' One batch only
  * @param conct     Connection flag: 'C' current data block is a continuation,
  * 'N' no connection between data blocks
  * @param ctrl      System order confirmation: 'C' user confirmation, 'N' no confirmation
  * @param nobr      Number of block rows (s)
  * @param m         Number of system inputs
  * @param l         Number of system outputs
  * @param nsmp      Number of input/output samples (t)
  * @param u         Input data array of size ldu-by-m (col-major) or nsmp-by-m (row-major)
  * @param ldu       Leading dimension of u (>= nsmp if col-major, >= m if row-major)
  * @param y         Input data array containing output measurements, size ldy-by-l (col-major) or nsmp-by-l (row-major)
  * @param ldy       Leading dimension of y (>= nsmp if col-major, >= l if row-major)
  * @param n         [in/out] On entry: 0 for automatic order detection or >0 for fixed order
  * On exit: The determined or user-specified system order
  * @param r         Triangular factor array (column-major format), size ldr-by-(2*(m+l)*nobr)
  * @param ldr       Leading dimension of r (see Fortran docs for minimum requirement)
  * @param sv        Output array of singular values used for order determination, size (l*nobr)
  * @param rcond     Tolerance used for estimating the rank of matrices in N4SID (METH='N').
  * @param tol       Tolerance used for estimating the system order based on singular values.
  * @param iwarn     Output warning flag, set on exit if a warning occurred.
  * @param row_major Non-zero if input/output matrices U, Y are in row-major format (C standard),
  * zero if in column-major format (Fortran standard). Note: R is always output in column-major.
  *
  * @return Zero on success, negative value indicates an error in arguments,
  * positive value indicates a computational error (see Fortran docs),
  * SLICOT_MEMORY_ERROR (-1010) if internal memory allocation failed.
  */
 SLICOT_EXPORT
 int slicot_ib01ad(char meth, char alg, char jobd, char batch, char conct, char ctrl,
                   int nobr, int m, int l, int nsmp,
                   double *u, int ldu, double *y, int ldy,
                   int *n, double *r, int ldr, double *sv, double rcond,
                   double tol, int *iwarn, int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* IB01AD_H */
 