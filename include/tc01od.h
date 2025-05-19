/**
 * @file tc01od.h
 * @brief Header for C wrapper of SLICOT routine TC01OD.
 */

#ifndef TC01OD_H
#define TC01OD_H

#include "slicot_utils.h" // Provides SLICOT_EXPORT macro

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Computes the dual of a left/right polynomial matrix representation.
 * @details Given Q(s)P(s)^-1 (right) or P(s)^-1Q(s) (left), computes the dual.
 * The coefficient arrays PCOEFF and QCOEFF are modified in-place.
 * This is a C wrapper for the SLICOT Fortran routine TC01OD.
 *
 * @param leri   (input) Character, specifies if a left or right matrix fraction is input:
 * = 'L': Left matrix fraction P(s)^-1Q(s) is input.
 * PCOEFF is PxP, QCOEFF is PxM. Output QCOEFF is MxP.
 * = 'R': Right matrix fraction Q(s)P(s)^-1 is input.
 * PCOEFF is MxM, QCOEFF is PxM. Output QCOEFF is MxP.
 * @param m      (input) The number of system inputs. m >= 0.
 * @param p      (input) The number of system outputs. p >= 0.
 * @param indlim (input) The highest K for which PCOEFF(.,.,K) and QCOEFF(.,.,K) are transposed.
 * Effectively, number of coefficient matrices (max_degree + 1). indlim >= 1.
 * @param pcoeff (input/output) Double precision 3D array for P(s) coefficients.
 * Dimensions depend on leri:
 * If leri='L', pcoeff is P x P x INDLIM.
 * If leri='R', pcoeff is M x M x INDLIM.
 * On exit, contains coefficients of P'(s) (transposed P(s)).
 * Stored slice-by-slice (k index varies slowest in C flat array).
 * @param ldpcoeff_c_rows (input) Leading dimension of pcoeff in C (number of rows of each 2D slice).
 * If row_major=0 (column-major C), this is the Fortran LD1 (>= porm).
 * If row_major=1 (row-major C), this is the number of rows in the C slice (= porm).
 * @param ldpcoeff_c_cols (input) Second dimension of pcoeff in C (number of columns of each 2D slice).
 * If row_major=0 (column-major C), this is the Fortran LD2 (>= porm).
 * If row_major=1 (row-major C), this is the number of columns in the C slice (= porm).
 * @param qcoeff (input/output) Double precision 3D array for Q(s) coefficients.
 * Input:  P x M x INDLIM.
 * Output: M x P x INDLIM (coefficients of Q'(s)).
 * Stored slice-by-slice (k index varies slowest in C flat array).
 * Caller must ensure qcoeff buffer is large enough for the output dimensions if M != P.
 * @param ldqcoeff_c_rows (input) Leading dimension of qcoeff in C (number of rows of each 2D slice).
 * If row_major=0, must be >= MAX(M,P) (Fortran LD1).
 * If row_major=1, for C buffer, this is the number of rows allocated for a slice (>= MAX(P,M)).
 * @param ldqcoeff_c_cols (input) Second dimension of qcoeff in C (number of columns of each 2D slice).
 * If row_major=0, must be >= MAX(M,P) (Fortran LD2).
 * If row_major=1, for C buffer, this is the number of columns allocated for a slice (>= MAX(P,M)).
 * @param row_major (input) Specifies matrix storage for 2D slices within PCOEFF and QCOEFF:
 * 0 for column-major (Fortran default for slices).
 * 1 for row-major (C default for slices).
 *
 * @return info Error indicator:
 * = 0:  successful exit
 * < 0:  if info = -i, the i-th argument had an illegal value.
 * = SLICOT_MEMORY_ERROR (-1010): internal memory allocation failed.
 */
SLICOT_EXPORT
int slicot_tc01od(char leri, int m, int p, int indlim,
                  double* pcoeff, int ldpcoeff_c_rows, int ldpcoeff_c_cols,
                  double* qcoeff, int ldqcoeff_c_rows, int ldqcoeff_c_cols,
                  int row_major);

#ifdef __cplusplus
}
#endif

#endif // TC01OD_H
