/**
 * @file mc01td.h
 * @brief C wrapper for SLICOT routine MC01TD
 *
 * This file provides a C interface to the SLICOT routine MC01TD,
 * which checks the stability of a given real polynomial.
 */

#ifndef MC01TD_H
#define MC01TD_H

#include <stddef.h> // For size_t
#include "slicot_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Checks stability of a given real polynomial.
 *
 * Determines whether a polynomial P(x) with real coefficients is stable,
 * either in the continuous-time sense (zeros in left half-plane) or
 * discrete-time sense (zeros inside the unit circle).
 *
 * @param[in] dico      Specifies the type of stability check:
 *                      = 'C': Continuous-time stability.
 *                      = 'D': Discrete-time stability.
 * @param[in] dp        The degree of the polynomial P(x), dp >= 0.
 * @param[in] p         Double array, dimension (dp+1). Coefficients of P(x)
 *                      in increasing powers of x (p[0] = constant term, p[dp] = leading coeff).
 * @param[out] stable   Integer flag (0 or 1). Set to 1 (.TRUE.) if P(x) is stable,
 *                      0 (.FALSE.) otherwise.
 * @param[out] nz       The number of unstable zeros found (zeros in RHP or outside unit circle).
 *                      Not determined if INFO = 2.
 * @param[out] iwarn    Warning indicator:
 *                      = 0: no warning.
 *                      = k: degree was reduced by k due to zero leading coefficients.
 *
 * @return info         Error indicator:
 *                      = 0: successful exit
 *                      < 0: if info = -i, the i-th argument had an illegal value.
 *                      = 1: P(x) is the zero polynomial.
 *                      = 2: Stability inconclusive (zeros potentially very close to boundary). NZ not determined.
 */
SLICOT_EXPORT
int slicot_mc01td(char dico, int dp, const double* p,
                  int* stable, int* nz, int* iwarn);

#ifdef __cplusplus
}
#endif

#endif /* MC01TD_H */
