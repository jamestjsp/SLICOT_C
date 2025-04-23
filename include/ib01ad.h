/**
 * @file ib01ad.h
 * @brief Header for C wrapper of SLICOT routine IB01AD
 *
 * This header provides the interface for the C wrapper of the
 * SLICOT routine IB01AD, which estimates the system order and
 * computes the triangular factor of the concatenated block-Hankel
 * matrices for subsequent system identification.
 */

#ifndef SLICOT_WRAPPER_IB01AD_H
#define SLICOT_WRAPPER_IB01AD_H

#include "slicot_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Estimates system order and computes the triangular factor of block-Hankel matrices
 * 
 * This routine implements a combined MOESP and N4SID subspace identification algorithm.
 * It estimates the system order and computes the triangular factor of the concatenated 
 * block-Hankel matrices for later use in system identification.
 * 
 * @param meth    Method to use: 'M' for MOESP, 'N' for N4SID algorithm
 * @param alg     Algorithm: 'C' for Cholesky, 'F' for Fast QR, 'Q' for QR algorithm
 * @param jobd    Specifies whether or not the matrices B and D should be computed using MOESP
 *                approach: 'M' to later use MOESP, 'N' don't use MOESP
 * @param batch   Data processing approach: 'F' First batch, 'I' Intermediate batch, 
 *                'L' Last batch, 'O' One batch only
 * @param conct   Connection flag: 'C' current data block is a continuation, 
 *                'N' no connection between data blocks
 * @param ctrl    System order confirmation: 'C' user confirmation, 'N' no confirmation
 * @param nobr    Number of block rows
 * @param m       Number of system inputs
 * @param l       Number of system outputs
 * @param nsmp    Number of input/output samples
 * @param u       Input data array of size ldu-by-m (col-major) or nsmp-by-m (row-major)
 * @param ldu     Leading dimension of u (>= nsmp if col-major, >= m if row-major)
 * @param y       Input data array containing output measurements, size ldy-by-l (col-major) or nsmp-by-l (row-major)  
 * @param ldy     Leading dimension of y (>= nsmp if col-major, >= l if row-major)
 * @param n       [in/out] On entry: 0 for automatic order detection or >0 for fixed order
 *                On exit: The determined or user-specified system order
 * @param r       Triangular factor array (column-major format) 
 * @param ldr     Leading dimension of r
 * @param sv      Array of singular values used for order determination
 * @param rcond   Used for estimating the rank of matrices
 * @param tol     Tolerance used for estimating the rank of matrices
 * @param iwork   Integer workspace array
 * @param dwork   Double workspace array
 * @param ldwork  Size of dwork array
 * @param iwarn   Output warning flag, set on exit if a warning occurred
 * @param row_major Non-zero if input/output matrices are in row-major format (C standard), 
 *                  zero if in column-major format (Fortran standard)
 * 
 * @return Zero on success, otherwise an error code
 */
SLICOT_C_WRAPPER_API
int slicot_ib01ad(char meth, char alg, char jobd, char batch, char conct, char ctrl, 
                 int nobr, int m, int l, int nsmp, 
                 double *u, int ldu, double *y, int ldy,
                 int *n, double *r, int ldr, double *sv, double rcond, 
                 double tol, int *iwork, double *dwork, int ldwork, 
                 int *iwarn, int row_major);

#ifdef __cplusplus
}
#endif

#endif /* SLICOT_WRAPPER_IB01AD_H */
