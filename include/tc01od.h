/**
 * @file tc01od.h
 * @brief C wrapper for SLICOT routine TC01OD
 *
 * This file provides a C interface to the SLICOT routine TC01OD,
 * which finds the dual right (left) polynomial matrix representation
 * of a given left (right) polynomial matrix representation.
 */

 #ifndef TC01OD_H
 #define TC01OD_H
 
 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Finds the dual of a left/right polynomial matrix representation.
  *
  * Given inv(P(s))*Q(s) (left), finds Q'(s)*inv(P'(s)) (right dual).
  * Given Q(s)*inv(P(s)) (right), finds inv(P'(s))*Q'(s) (left dual).
  * The operation essentially transposes the coefficient matrices.
  *
  * @param[in] leri      Specifies input type: 'L' (left fraction) or 'R' (right fraction).
  * @param[in] m         Number of system inputs (columns of Q, dimension of P if leri='R'), m >= 0.
  * @param[in] p         Number of system outputs (rows of Q, dimension of P if leri='L'), p >= 0.
  * @param[in] indlim    Highest power index + 1 for coefficients (kpcoef + 1), indlim >= 1.
  * @param[in,out] pcoeff Double array, dimension (ldpco1, ldpco2, indlim).
  * On entry, coefficients of P(s). Dim (p,p,indlim) if 'L', (m,m,indlim) if 'R'.
  * On exit, coefficients of P'(s). Dim (m,m,indlim) if 'L', (p,p,indlim) if 'R'.
  * @param[in] ldpco1    First leading dimension of pcoeff. >= max(1,p) if 'L', >= max(1,m) if 'R'.
  * @param[in] ldpco2    Second leading dimension of pcoeff. >= max(1,p) if 'L', >= max(1,m) if 'R'.
  * @param[in,out] qcoeff Double array, dimension (ldqco1, ldqco2, indlim).
  * On entry, coefficients of Q(s). Dim (p,m,indlim).
  * On exit, coefficients of Q'(s). Dim (m,p,indlim).
  * @param[in] ldqco1    First leading dimension of qcoeff. >= max(1,m,p).
  * @param[in] ldqco2    Second leading dimension of qcoeff. >= max(1,m,p).
  * @param[in] row_major Integer flag: 0 for column-major, 1 for row-major.
  * Affects interpretation of ldpco*, ldqco*.
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value
  * Memory allocation errors may also be returned.
  */
 SLICOT_C_WRAPPER_API
 int slicot_tc01od(char leri, int m, int p, int indlim,
                   double* pcoeff, int ldpco1, int ldpco2,
                   double* qcoeff, int ldqco1, int ldqco2,
                   int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* TC01OD_H */
 