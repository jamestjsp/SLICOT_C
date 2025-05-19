/**
 * @file tg01ad.c
 * @brief C wrapper for SLICOT routine TG01AD.
 * @details Balances the matrices of a system pencil.
 * Workspace is allocated internally. Matrices A, E, B, C are modified in-place.
 */

#include <stdlib.h>
#include <string.h> // For memcpy
#include <ctype.h>  // For toupper
#include <stdio.h>  // For debugging if needed

#include "tg01ad.h"
#include "slicot_utils.h" // Provides MAX, CHECK_ALLOC, transpose functions
#include "slicot_f77.h"   // Provides F77_FUNC

/* External Fortran routine declaration */
extern void F77_FUNC(tg01ad, TG01AD)(
    const char* job, const int* l, const int* n, const int* m, const int* p,
    const double* thresh,
    double* a, const int* lda,
    double* e, const int* lde,
    double* b, const int* ldb,
    double* c, const int* ldc,
    double* lscale, double* rscale,
    double* dwork, int* info,
    size_t job_len);

SLICOT_EXPORT
int slicot_tg01ad(char job_in, int l_in, int n_in, int m_in, int p_in, double thresh_in,
                  double* a_c, int lda_c,
                  double* e_c, int lde_c,
                  double* b_c, int ldb_c,
                  double* c_c, int ldc_c,
                  double* lscale_out, double* rscale_out,
                  int row_major)
{
    int info = 0;
    char job_upper;

    // Fortran-equivalent parameters (inputs are not modified by Fortran, except matrices)
    int l_f = l_in;
    int n_f = n_in;
    int m_f = m_in;
    int p_f = p_in;

    // Pointers for Fortran call (these will point to original or _cm buffers)
    double* a_ptr = a_c;
    double* e_ptr = e_c;
    double* b_ptr = b_c;
    double* c_ptr = c_c;
    // lscale_out and rscale_out are directly passed as they are 1D

    // Leading dimensions for Fortran call
    int lda_f = lda_c;
    int lde_f = lde_c;
    int ldb_f = ldb_c;
    int ldc_f = ldc_c;

    // Temporary column-major copies if row_major is used
    double *a_cm = NULL;
    double *e_cm = NULL;
    double *b_cm = NULL;
    double *c_cm = NULL;

    // Workspace
    double *dwork = NULL;
    int ldwork_calc;
    
    double dummy_double_array[1] = {0.0}; // For safe passing if dimension is 0

    // --- 1. Validate input parameters ---
    job_upper = toupper(job_in);
    if (job_upper != 'A' && job_upper != 'B' && job_upper != 'C' && job_upper != 'N') {
        info = -1; goto cleanup;
    }

    if (l_f < 0) { info = -2; goto cleanup; }
    if (n_f < 0) { info = -3; goto cleanup; }
    if (m_f < 0) { info = -4; goto cleanup; }
    if (p_f < 0) { info = -5; goto cleanup; }
    if (thresh_in < 0.0) { info = -6; goto cleanup; }

    // Validate A, LDA, E, LDE (only if l_f > 0 and n_f > 0)
    if (l_f > 0 && n_f > 0) {
        if (a_c == NULL) { info = -7; goto cleanup; }
        if (row_major) { if (lda_c < MAX(1, n_f)) { info = -8; goto cleanup; } }
        else { if (lda_c < MAX(1, l_f)) { info = -8; goto cleanup; } }

        if (e_c == NULL) { info = -9; goto cleanup; }
        if (row_major) { if (lde_c < MAX(1, n_f)) { info = -10; goto cleanup; } }
        else { if (lde_c < MAX(1, l_f)) { info = -10; goto cleanup; } }
    }

    // Validate B, LDB (only if l_f > 0 and m_f > 0)
    if (m_f > 0) { // B is referenced if M > 0
        if (l_f > 0 ) { // B has L rows
             if (b_c == NULL) { info = -11; goto cleanup; }
            if (row_major) { if (ldb_c < MAX(1, m_f)) { info = -12; goto cleanup; } }
            else { if (ldb_c < MAX(1, l_f)) { info = -12; goto cleanup; } }
        }
    } else { // M = 0, B not referenced, LDB must be >= 1
        if (ldb_c < 1 && b_c != NULL) {info = -12; goto cleanup;} // Fortran requires LDB >=1 even if M=0
    }


    // Validate C, LDC (only if p_f > 0 and n_f > 0)
    if (p_f > 0 ) { // C is referenced if P > 0
        if (n_f > 0) { // C has N columns
            if (c_c == NULL) { info = -13; goto cleanup; }
            if (row_major) { if (ldc_c < MAX(1, n_f)) { info = -14; goto cleanup; } }
            else { if (ldc_c < MAX(1, p_f)) { info = -14; goto cleanup; } }
        }
    } else { // P = 0, C not referenced, LDC must be >=1
         if (ldc_c < 1 && c_c != NULL) { info = -14; goto cleanup;}
    }


    if (l_f > 0 && lscale_out == NULL) { info = -15; goto cleanup; }
    if (n_f > 0 && rscale_out == NULL) { info = -16; goto cleanup; }
    
    if (info != 0) { goto cleanup; }

    // --- 2. Calculate workspace size and allocate DWORK ---
    // DWORK dimension (3*(L+N))
    ldwork_calc = MAX(1, 3 * (l_f + n_f));
    dwork = (double*)malloc((size_t)ldwork_calc * sizeof(double));
    CHECK_ALLOC(dwork);

    // --- 3. Handle row_major conversions ---
    size_t a_e_size_elems = (size_t)l_f * n_f;
    size_t b_size_elems = (size_t)l_f * m_f;
    size_t c_size_elems = (size_t)p_f * n_f;

    if (row_major) {
        // Fortran LDs are number of rows
        lda_f = MAX(1, l_f);
        lde_f = MAX(1, l_f);
        ldb_f = MAX(1, l_f); // if m_f > 0
        ldc_f = MAX(1, p_f); // if p_f > 0

        if (l_f > 0 && n_f > 0) {
            if (a_c != NULL) {
                a_cm = (double*)malloc(a_e_size_elems * sizeof(double)); CHECK_ALLOC(a_cm);
                slicot_transpose_to_fortran_with_ld(a_c, a_cm, l_f, n_f, lda_c, lda_f, sizeof(double));
                a_ptr = a_cm;
            } else { a_ptr = NULL; } // Should have been caught by validation

            if (e_c != NULL) {
                e_cm = (double*)malloc(a_e_size_elems * sizeof(double)); CHECK_ALLOC(e_cm);
                slicot_transpose_to_fortran_with_ld(e_c, e_cm, l_f, n_f, lde_c, lde_f, sizeof(double));
                e_ptr = e_cm;
            } else { e_ptr = NULL; }
        } else {
            a_ptr = NULL; e_ptr = NULL;
        }

        if (l_f > 0 && m_f > 0 && b_c != NULL) {
            b_cm = (double*)malloc(b_size_elems * sizeof(double)); CHECK_ALLOC(b_cm);
            slicot_transpose_to_fortran_with_ld(b_c, b_cm, l_f, m_f, ldb_c, ldb_f, sizeof(double));
            b_ptr = b_cm;
        } else { b_ptr = NULL; if (m_f == 0) ldb_f = 1; }


        if (p_f > 0 && n_f > 0 && c_c != NULL) {
            c_cm = (double*)malloc(c_size_elems * sizeof(double)); CHECK_ALLOC(c_cm);
            slicot_transpose_to_fortran_with_ld(c_c, c_cm, p_f, n_f, ldc_c, ldc_f, sizeof(double));
            c_ptr = c_cm;
        } else { c_ptr = NULL; if (p_f == 0) ldc_f = 1; }

    } else { // Column-major C input
        // Assign pointers directly, ensure NULL if dimension is zero
        a_ptr = (l_f > 0 && n_f > 0 && a_c != NULL) ? a_c : NULL;
        e_ptr = (l_f > 0 && n_f > 0 && e_c != NULL) ? e_c : NULL;
        b_ptr = (l_f > 0 && m_f > 0 && b_c != NULL) ? b_c : NULL;
        c_ptr = (p_f > 0 && n_f > 0 && c_c != NULL) ? c_c : NULL;
        
        lda_f = lda_c; lde_f = lde_c; 
        ldb_f = (m_f > 0) ? ldb_c : 1;
        ldc_f = (p_f > 0) ? ldc_c : 1;
    }
    
    // Ensure Fortran LDs are at least 1
    lda_f = MAX(1, lda_f); lde_f = MAX(1, lde_f);
    ldb_f = MAX(1, ldb_f); ldc_f = MAX(1, ldc_f);

    // Use dummy arrays for Fortran if logically zero-sized but Fortran expects non-NULL
    double* final_a_ptr = (l_f > 0 && n_f > 0 && a_ptr != NULL) ? a_ptr : dummy_double_array;
    double* final_e_ptr = (l_f > 0 && n_f > 0 && e_ptr != NULL) ? e_ptr : dummy_double_array;
    double* final_b_ptr = (l_f > 0 && m_f > 0 && b_ptr != NULL) ? b_ptr : dummy_double_array;
    double* final_c_ptr = (p_f > 0 && n_f > 0 && c_ptr != NULL) ? c_ptr : dummy_double_array;
    double* final_lscale_ptr = (l_f > 0 && lscale_out != NULL) ? lscale_out : dummy_double_array;
    double* final_rscale_ptr = (n_f > 0 && rscale_out != NULL) ? rscale_out : dummy_double_array;


    // --- 4. Call Fortran routine ---
    F77_FUNC(tg01ad, TG01AD)(&job_upper, &l_f, &n_f, &m_f, &p_f, &thresh_in,
                             final_a_ptr, &lda_f, final_e_ptr, &lde_f,
                             final_b_ptr, &ldb_f, final_c_ptr, &ldc_f,
                             final_lscale_ptr, final_rscale_ptr,
                             dwork, &info, (size_t)1);

    // --- 5. Transpose results back if row_major and successful ---
    if (row_major && info == 0) {
        if (l_f > 0 && n_f > 0) {
            if (a_cm != NULL && a_c != NULL) slicot_transpose_to_c_with_ld(a_cm, a_c, l_f, n_f, lda_f, lda_c, sizeof(double));
            if (e_cm != NULL && e_c != NULL) slicot_transpose_to_c_with_ld(e_cm, e_c, l_f, n_f, lde_f, lde_c, sizeof(double));
        }
        if (l_f > 0 && m_f > 0 && b_cm != NULL && b_c != NULL) {
            slicot_transpose_to_c_with_ld(b_cm, b_c, l_f, m_f, ldb_f, ldb_c, sizeof(double));
        }
        if (p_f > 0 && n_f > 0 && c_cm != NULL && c_c != NULL) {
            slicot_transpose_to_c_with_ld(c_cm, c_c, p_f, n_f, ldc_f, ldc_c, sizeof(double));
        }
        // lscale and rscale are 1D, no transpose needed.
    }

cleanup:
    // --- 6. Free allocated memory ---
    free(dwork);
    if (row_major) {
        free(a_cm);
        free(e_cm);
        free(b_cm);
        free(c_cm);
    }
    return info;
}
