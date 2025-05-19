/**
 * @file tc04ad.h
 * @brief Header for C wrapper of SLICOT routine TC04AD.
 */

#ifndef TC04AD_H
#define TC04AD_H

#include "slicot_utils.h" // Provides SLICOT_EXPORT macro

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Converts a polynomial matrix representation to a state-space representation.
 * @details Computes (A,B,C,D) such that C*inv(sI-A)*B + D = P(s)^-1Q(s) or Q(s)P(s)^-1.
 * This is a C wrapper for the SLICOT Fortran routine TC04AD.
 * Workspace is allocated internally.
 *
 * @param leri   (input) Character, specifies if a left or right matrix fraction is input:
 * = 'L': Left matrix fraction P(s)^-1Q(s) is input.
 * = 'R': Right matrix fraction Q(s)P(s)^-1 is input.
 * @param m_c    (input) The number of system inputs. m_c >= 0.
 * @param p_c    (input) The number of system outputs. p_c >= 0.
 * @param index_c (input) Integer array.
 * If leri='L', dimension p_c, degrees of rows of P(s).
 * If leri='R', dimension m_c, degrees of columns of P(s).
 * @param pcoeff_c (input) Double precision 3D array for P(s) coefficients.
 * Stored slice-by-slice. Dimensions depend on leri and index_c.
 * If leri='R', it's modified by Fortran but restored on exit.
 * @param ldpcoeff_c_rows (input) C leading dimension (rows) of each 2D slice of pcoeff_c.
 * @param ldpcoeff_c_cols (input) C second dimension (cols) of each 2D slice of pcoeff_c.
 * @param qcoeff_c (input) Double precision 3D array for Q(s) coefficients.
 * Stored slice-by-slice. Dimensions depend on leri and index_c.
 * If leri='R', it's modified by Fortran but restored on exit.
 * @param ldqcoeff_c_rows (input) C leading dimension (rows) of each 2D slice of qcoeff_c.
 * @param ldqcoeff_c_cols (input) C second dimension (cols) of each 2D slice of qcoeff_c.
 * @param[out] n_out Pointer to store the order of the state-space representation (sum of index_c).
 * @param[out] rcond_out Pointer to store the estimated reciprocal condition number of P(s)'s leading coefficient matrix.
 * @param[out] a_out Double precision array to store the state matrix A (n_out x n_out).
 * @param lda_c   Leading dimension of a_out in C.
 * @param[out] b_out Double precision array to store the input matrix B (n_out x m_c).
 * Fortran uses a larger buffer (n_out x MAX(m_c,p_c)).
 * @param ldb_c   Leading dimension of b_out in C.
 * @param[out] c_out Double precision array to store the output matrix C (p_c x n_out).
 * Fortran uses a larger buffer (MAX(m_c,p_c) x n_out).
 * @param ldc_c   Leading dimension of c_out in C.
 * @param[out] d_out Double precision array to store the direct transmission matrix D (p_c x m_c).
 * Fortran uses a larger buffer (MAX(m_c,p_c) x MAX(m_c,p_c)).
 * @param ldd_c   Leading dimension of d_out in C.
 * @param row_major_flag (input) Specifies matrix storage for all 2D/3D arrays:
 * 0 for column-major C layout (Fortran-like slices).
 * 1 for row-major C layout (standard C slices).
 *
 * @return info Error indicator:
 * = 0:  successful exit
 * < 0:  if info = -i, the i-th argument had an illegal value.
 * = 1:  P(s) is not row/column proper.
 * = SLICOT_MEMORY_ERROR (-1010): internal memory allocation failed.
 */
SLICOT_EXPORT
int slicot_tc04ad(char leri, int m_c, int p_c, const int* index_c,
                  const double* pcoeff_c, int ldpcoeff_c_rows, int ldpcoeff_c_cols,
                  const double* qcoeff_c, int ldqcoeff_c_rows, int ldqcoeff_c_cols,
                  int* n_out, double* rcond_out,
                  double* a_out, int lda_c,
                  double* b_out, int ldb_c,
                  double* c_out, int ldc_c,
                  double* d_out, int ldd_c,
                  int row_major_flag);

#ifdef __cplusplus
}
#endif

#endif // TC04AD_H
