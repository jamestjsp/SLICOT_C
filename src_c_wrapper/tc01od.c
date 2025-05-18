/**
 * @file tc01od.c
 * @brief C wrapper for SLICOT routine TC01OD.
 * @details Computes the dual of a left/right polynomial matrix representation.
 * PCOEFF and QCOEFF are modified in-place.
 */

#include <stdlib.h>
#include <ctype.h>
#include <string.h> // For memcpy

#include "tc01od.h"
#include "slicot_utils.h" // For CHECK_ALLOC, MAX, MIN, transpose functions
#include "slicot_f77.h"   // For F77_FUNC

// Fortran routine declaration
extern void F77_FUNC(tc01od, TC01OD)(
    const char* leri, const int* m, const int* p, const int* indlim,
    double* pcoeff, const int* ldpco1, const int* ldpco2,
    double* qcoeff, const int* ldqco1, const int* ldqco2,
    int* info,
    int leri_len);

SLICOT_EXPORT
int slicot_tc01od(char leri_c, int m_c, int p_c, int indlim_c,
                  double* pcoeff_c, int ldpcoeff_c_rows, int ldpcoeff_c_cols,
                  double* qcoeff_c, int ldqcoeff_c_rows, int ldqcoeff_c_cols,
                  int row_major_flag) {
    int info = 0;
    char leri_f = toupper(leri_c);

    double* pcoeff_f_ptr = pcoeff_c;
    double* qcoeff_f_ptr = qcoeff_c;

    double* pcoeff_cm = NULL;
    double* qcoeff_cm = NULL;

    // Validate inputs
    if (leri_f != 'L' && leri_f != 'R') { info = -1; goto cleanup; }
    if (m_c < 0) { info = -2; goto cleanup; }
    if (p_c < 0) { info = -3; goto cleanup; }
    if (indlim_c < 1) { info = -4; goto cleanup; }

    int porm_f = (leri_f == 'L') ? p_c : m_c; // Dimension for PCOEFF (PxP or MxM)
    
    // Fortran leading dimensions for PCOEFF (LD1_F = rows, LD2_F = cols for Fortran view of PCOEFF)
    int ldpco1_f_expected = MAX(1, porm_f); 
    int ldpco2_f_expected = MAX(1, porm_f);

    // Fortran leading dimensions for QCOEFF (LD1_F = rows, LD2_F = cols for Fortran view of QCOEFF)
    // Fortran's QCOEFF is dimensioned (LDQCO1, LDQCO2, INDLIM) where LDQCO1 >= MAX(1,M,P), LDQCO2 >= MAX(1,M,P)
    int ldqco1_f_expected = MAX(1, MAX(m_c, p_c));
    int ldqco2_f_expected = MAX(1, MAX(m_c, p_c));

    // Validate C leading dimensions
    if (row_major_flag) {
        // For PCOEFF, C slice is porm_f x porm_f. C ldpcoeff_c_cols is the number of columns in C's RM slice.
        if (pcoeff_c != NULL && porm_f > 0 && (ldpcoeff_c_rows < porm_f || ldpcoeff_c_cols < porm_f)) { info = -6; goto cleanup; }
        // For QCOEFF, C buffer slice must be at least MAX(P,M) x MAX(P,M) to hold input/output.
        // C ldqcoeff_c_cols is the number of columns in C's RM slice.
        if (qcoeff_c != NULL && (p_c > 0 || m_c > 0) && (ldqcoeff_c_rows < MAX(p_c, m_c) || ldqcoeff_c_cols < MAX(p_c, m_c))) { info = -9; goto cleanup; }
    } else { // Column-major C
        // C ldpcoeff_c_rows is Fortran's LD1, ldpcoeff_c_cols is Fortran's LD2
        if (pcoeff_c != NULL && porm_f > 0 && (ldpcoeff_c_rows < ldpco1_f_expected || ldpcoeff_c_cols < ldpco2_f_expected)) { info = -6; goto cleanup; }
        // C ldqcoeff_c_rows is Fortran's LD1, ldqcoeff_c_cols is Fortran's LD2
        if (qcoeff_c != NULL && (p_c > 0 || m_c > 0) && (ldqcoeff_c_rows < ldqco1_f_expected || ldqcoeff_c_cols < ldqco2_f_expected)) { info = -9; goto cleanup; }
    }
    
    if (pcoeff_c == NULL && porm_f > 0) { info = -5; goto cleanup; }
    if (qcoeff_c == NULL && (p_c > 0 || m_c > 0) ) { info = -8; goto cleanup; }

    if (info != 0) goto cleanup;

    if (row_major_flag) {
        // --- PCOEFF handling ---
        if (porm_f > 0 && pcoeff_c != NULL) {
            size_t pcoeff_slice_elems_f = (size_t)porm_f * porm_f; // Elements in one Fortran slice (porm x porm)
            pcoeff_cm = (double*)malloc(pcoeff_slice_elems_f * indlim_c * sizeof(double));
            CHECK_ALLOC(pcoeff_cm);
            pcoeff_f_ptr = pcoeff_cm;

            for (int k = 0; k < indlim_c; ++k) {
                double* c_slice_k_start = pcoeff_c + k * (size_t)ldpcoeff_c_rows * ldpcoeff_c_cols;
                double* f_slice_k_start = pcoeff_cm + k * pcoeff_slice_elems_f;
                // C slice (porm_f x porm_f), C ld_src (cols) = ldpcoeff_c_cols
                // Fortran slice (porm_f x porm_f), Fortran ld_dest (rows) = ldpco1_f_expected
                slicot_transpose_to_fortran_with_ld(c_slice_k_start, f_slice_k_start, 
                                                    porm_f, porm_f, 
                                                    ldpcoeff_c_cols,    // ld_src for C (cols)
                                                    ldpco1_f_expected,  // ld_dest for Fortran (rows)
                                                    sizeof(double));
            }
        } else if (porm_f == 0) { // If pcoeff is logically empty
            pcoeff_f_ptr = NULL; 
        }

        // --- QCOEFF handling ---
        // Input is PxM, output is MxP. Fortran buffer is LDQCO1_F x LDQCO2_F.
        size_t qcoeff_slice_elems_f = (size_t)ldqco1_f_expected * ldqco2_f_expected;

        if ((p_c > 0 || m_c > 0) && qcoeff_c != NULL) { 
            qcoeff_cm = (double*)malloc(qcoeff_slice_elems_f * indlim_c * sizeof(double));
            CHECK_ALLOC(qcoeff_cm);
            qcoeff_f_ptr = qcoeff_cm;
            memset(qcoeff_cm, 0, qcoeff_slice_elems_f * indlim_c * sizeof(double)); 

            // Transpose PxM input from user's C row-major qcoeff_c to top-left of qcoeff_cm slices
            if (p_c > 0 && m_c > 0) { // Only if there's actual input data
                 for (int k = 0; k < indlim_c; ++k) {
                    double* c_slice_k_start = qcoeff_c + k * (size_t)ldqcoeff_c_rows * ldqcoeff_c_cols;
                    double* f_slice_k_start = qcoeff_cm + k * qcoeff_slice_elems_f;
                    // C slice is p_c x m_c, C ld_src (cols) = ldqcoeff_c_cols
                    // Fortran slice is ldqco1_f_expected x ldqco2_f_expected, Fortran ld_dest (rows) = ldqco1_f_expected
                    slicot_transpose_to_fortran_with_ld(c_slice_k_start, f_slice_k_start,
                                                        p_c, m_c,
                                                        ldqcoeff_c_cols,    // ld_src for C (cols)
                                                        ldqco1_f_expected,  // ld_dest for Fortran (rows)
                                                        sizeof(double));
                }
            }
        } else { // If qcoeff is logically empty
             qcoeff_f_ptr = NULL;
        }
    } else { // Column-major C
        // Fortran LDs are taken directly from C LDs if C is column-major
        ldpco1_f_expected = ldpcoeff_c_rows;
        ldpco2_f_expected = ldpcoeff_c_cols;
        ldqco1_f_expected = ldqcoeff_c_rows;
        ldqco2_f_expected = ldqcoeff_c_cols;
        
        if (porm_f == 0) pcoeff_f_ptr = NULL;
        if (p_c == 0 && m_c == 0) qcoeff_f_ptr = NULL;
    }

    // Call Fortran routine
    F77_FUNC(tc01od, TC01OD)(&leri_f, &m_c, &p_c, &indlim_c,
                             pcoeff_f_ptr, &ldpco1_f_expected, &ldpco2_f_expected,
                             qcoeff_f_ptr, &ldqco1_f_expected, &ldqco2_f_expected,
                             &info, 1);

    if (row_major_flag && info == 0) {
        // --- PCOEFF transpose back ---
        if (porm_f > 0 && pcoeff_cm != NULL && pcoeff_c != NULL) {
            size_t pcoeff_slice_elems_f = (size_t)porm_f * porm_f;
            for (int k = 0; k < indlim_c; ++k) {
                double* f_slice_k_start = pcoeff_cm + k * pcoeff_slice_elems_f;
                double* c_slice_k_start = pcoeff_c + k * (size_t)ldpcoeff_c_rows * ldpcoeff_c_cols;
                // Fortran slice is porm_f x porm_f (col-major), Fortran ld_src (rows) = ldpco1_f_expected
                // C slice is porm_f x porm_f (row-major), C ld_dest (cols) = ldpcoeff_c_cols
                slicot_transpose_to_c_with_ld(f_slice_k_start, c_slice_k_start,
                                              porm_f, porm_f,
                                              ldpco1_f_expected, // ld_src for Fortran (rows)
                                              ldpcoeff_c_cols,   // ld_dest for C (cols)
                                              sizeof(double));
            }
        }

        // --- QCOEFF transpose back ---
        // Output is MxP. Fortran buffer was LDQCO1_F x LDQCO2_F.
        size_t qcoeff_slice_elems_f = (size_t)ldqco1_f_expected * ldqco2_f_expected;

        if ((m_c > 0 || p_c > 0) && qcoeff_cm != NULL && qcoeff_c != NULL) { 
             if (m_c > 0 && p_c > 0) { // Only if output has non-zero dimensions
                for (int k = 0; k < indlim_c; ++k) {
                    double* f_slice_k_start = qcoeff_cm + k * qcoeff_slice_elems_f;
                    double* c_slice_k_start = qcoeff_c + k * (size_t)ldqcoeff_c_rows * ldqcoeff_c_cols;
                    // Fortran slice (output part) is m_c x p_c (col-major), Fortran ld_src (rows) = ldqco1_f_expected
                    // C slice is m_c x p_c (row-major), C ld_dest (cols) = ldqcoeff_c_cols
                    slicot_transpose_to_c_with_ld(f_slice_k_start, c_slice_k_start,
                                                  m_c, p_c,
                                                  ldqco1_f_expected, // ld_src for Fortran (rows)
                                                  ldqcoeff_c_cols,   // ld_dest for C (cols)
                                                  sizeof(double));
                }
            }
        }
    }

cleanup:
    free(pcoeff_cm);
    free(qcoeff_cm);
    return info;
}
