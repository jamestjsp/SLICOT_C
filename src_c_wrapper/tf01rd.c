/**
 * @file tf01rd.c
 * @brief C wrapper for SLICOT routine TF01RD.
 * @details Computes N Markov parameters from (A,B,C).
 * Workspace is allocated internally.
 */

#include <stdlib.h>
#include <string.h> // For memcpy
#include <stdio.h>  // For debugging printf if needed

#include "tf01rd.h"
#include "slicot_utils.h" // Provides MAX, CHECK_ALLOC, transpose functions
#include "slicot_f77.h"   // Provides F77_FUNC

/* External Fortran routine declaration */
extern void F77_FUNC(tf01rd, TF01RD)(
    const int* na, const int* nb, const int* nc, const int* N_in,
    const double* a, const int* lda,
    const double* b, const int* ldb,
    const double* c, const int* ldc,
    double* h, const int* ldh,
    double* dwork, const int* ldwork,
    int* info);

SLICOT_EXPORT
int slicot_tf01rd(int na_in, int nb_in, int nc_in, int N_in,
                  const double* a_c, int lda_c,
                  const double* b_c, int ldb_c,
                  const double* c_c, int ldc_c,
                  double* h_c, int ldh_c,
                  int row_major)
{
    int info = 0;
    
    // Fortran-equivalent parameters
    int na_f = na_in;
    int nb_f = nb_in;
    int nc_f = nc_in;
    int N_f = N_in;

    // Pointers for Fortran call
    const double* a_ptr = a_c;
    const double* b_ptr = b_c;
    const double* c_ptr = c_c;
    double* h_ptr = h_c;

    // Leading dimensions for Fortran call
    int lda_f = lda_c;
    int ldb_f = ldb_c;
    int ldc_f = ldc_c;
    int ldh_f = ldh_c;

    // Temporary column-major copies if row_major is used
    double *a_cm = NULL;
    double *b_cm = NULL;
    double *c_cm = NULL;
    double *h_cm = NULL; // For output H

    // Workspace
    double *dwork = NULL;
    int ldwork_calc;

    // Dummy array for zero-sized inputs if Fortran expects non-NULL
    double dummy_double_array[1] = {0.0};


    // --- 1. Validate input parameters ---
    if (na_f < 0) { info = -1; goto cleanup; }
    if (nb_f < 0) { info = -2; goto cleanup; }
    if (nc_f < 0) { info = -3; goto cleanup; }
    if (N_f < 0)  { info = -4; goto cleanup; }

    // Validate A and LDA
    if (na_f > 0) {
        if (a_c == NULL) { info = -5; goto cleanup; }
        if (row_major) { // lda_c is num_cols for C row-major
            if (lda_c < MAX(1, na_f)) { info = -6; goto cleanup; }
        } else { // lda_c is num_rows for C col-major
            if (lda_c < MAX(1, na_f)) { info = -6; goto cleanup; }
        }
    }

    // Validate B and LDB
    if (na_f > 0 && nb_f > 0) {
        if (b_c == NULL) { info = -7; goto cleanup; }
        if (row_major) { // ldb_c is num_cols
            if (ldb_c < MAX(1, nb_f)) { info = -8; goto cleanup; }
        } else { // ldb_c is num_rows
            if (ldb_c < MAX(1, na_f)) { info = -8; goto cleanup; }
        }
    }
    
    // Validate C and LDC
    if (nc_f > 0 && na_f > 0) {
        if (c_c == NULL) { info = -9; goto cleanup; }
        if (row_major) { // ldc_c is num_cols
            if (ldc_c < MAX(1, na_f)) { info = -10; goto cleanup; }
        } else { // ldc_c is num_rows
            if (ldc_c < MAX(1, nc_f)) { info = -10; goto cleanup; }
        }
    }

    // Validate H and LDH
    // H is nc_f x (N_f * nb_f)
    int h_cols_f = N_f * nb_f;
    if (nc_f > 0 && h_cols_f > 0) { // Only check if H is non-empty
        if (h_c == NULL) { info = -11; goto cleanup; }
        if (row_major) { // ldh_c is num_cols
            if (ldh_c < MAX(1, h_cols_f)) { info = -12; goto cleanup; }
        } else { // ldh_c is num_rows
            if (ldh_c < MAX(1, nc_f)) { info = -12; goto cleanup; }
        }
    }
    
    if (info != 0) { goto cleanup; }

    // --- 2. Calculate workspace size and allocate DWORK ---
    // LDWORK >= MAX(1, 2*NA*NC)
    ldwork_calc = MAX(1, 2 * na_f * nc_f);
    dwork = (double*)malloc((size_t)ldwork_calc * sizeof(double));
    CHECK_ALLOC(dwork);

    // --- 3. Handle row_major conversions ---
    size_t a_size_elems = (size_t)na_f * na_f;
    size_t b_size_elems = (size_t)na_f * nb_f;
    size_t c_size_elems = (size_t)nc_f * na_f;
    size_t h_size_elems = (size_t)nc_f * h_cols_f;

    if (row_major) {
        // Fortran LDs are number of rows
        lda_f = MAX(1, na_f);
        ldb_f = MAX(1, na_f);
        ldc_f = MAX(1, nc_f);
        ldh_f = MAX(1, nc_f);

        if (na_f > 0 && a_c != NULL) {
            a_cm = (double*)malloc(a_size_elems * sizeof(double)); CHECK_ALLOC(a_cm);
            slicot_transpose_to_fortran_with_ld(a_c, a_cm, na_f, na_f, lda_c, lda_f, sizeof(double));
            a_ptr = a_cm;
        } else { a_ptr = NULL; } // Or pass dummy_double_array if Fortran expects non-NULL

        if (na_f > 0 && nb_f > 0 && b_c != NULL) {
            b_cm = (double*)malloc(b_size_elems * sizeof(double)); CHECK_ALLOC(b_cm);
            slicot_transpose_to_fortran_with_ld(b_c, b_cm, na_f, nb_f, ldb_c, ldb_f, sizeof(double));
            b_ptr = b_cm;
        } else { b_ptr = NULL; }

        if (nc_f > 0 && na_f > 0 && c_c != NULL) {
            c_cm = (double*)malloc(c_size_elems * sizeof(double)); CHECK_ALLOC(c_cm);
            slicot_transpose_to_fortran_with_ld(c_c, c_cm, nc_f, na_f, ldc_c, ldc_f, sizeof(double));
            c_ptr = c_cm;
        } else { c_ptr = NULL; }

        if (nc_f > 0 && h_cols_f > 0 && h_c != NULL) {
            h_cm = (double*)malloc(h_size_elems * sizeof(double)); CHECK_ALLOC(h_cm);
            h_ptr = h_cm; // Fortran will write to this column-major buffer
        } else { h_ptr = NULL; }

    } else { // Column-major C input
        // Ensure pointers are NULL if dimensions are zero, Fortran LDs are C LDs
        if (na_f == 0) { a_ptr = NULL; b_ptr = NULL; c_ptr = NULL; } // A, B, C depend on NA
        else {
             a_ptr = a_c; lda_f = lda_c;
             if (nb_f == 0) b_ptr = NULL; else {b_ptr = b_c; ldb_f = ldb_c;}
             if (nc_f == 0) c_ptr = NULL; else {c_ptr = c_c; ldc_f = ldc_c;} // C also depends on NC
        }
        if (nc_f == 0 || h_cols_f == 0) h_ptr = NULL; else {h_ptr = h_c; ldh_f = ldh_c;}
    }
    
    // Safety for Fortran: pass dummy if ptr is NULL but Fortran might expect non-NULL
    const double* final_a_ptr = (na_f > 0 && a_ptr != NULL) ? a_ptr : dummy_double_array;
    const double* final_b_ptr = (na_f > 0 && nb_f > 0 && b_ptr != NULL) ? b_ptr : dummy_double_array;
    const double* final_c_ptr = (nc_f > 0 && na_f > 0 && c_ptr != NULL) ? c_ptr : dummy_double_array;
    double* final_h_ptr       = (nc_f > 0 && h_cols_f > 0 && h_ptr != NULL) ? h_ptr : dummy_double_array;
    
    // Ensure LDs are at least 1 for Fortran call
    lda_f = MAX(1, lda_f);
    ldb_f = MAX(1, ldb_f);
    ldc_f = MAX(1, ldc_f);
    ldh_f = MAX(1, ldh_f);


    // --- 4. Call Fortran routine ---
    F77_FUNC(tf01rd, TF01RD)(&na_f, &nb_f, &nc_f, &N_f,
                             final_a_ptr, &lda_f,
                             final_b_ptr, &ldb_f,
                             final_c_ptr, &ldc_f,
                             final_h_ptr, &ldh_f,
                             dwork, &ldwork_calc,
                             &info);

    // --- 5. Transpose H back if row_major and successful ---
    if (row_major && info == 0 && nc_f > 0 && h_cols_f > 0 && h_cm != NULL && h_c != NULL) {
        slicot_transpose_to_c_with_ld(h_cm, h_c, nc_f, h_cols_f, ldh_f, ldh_c, sizeof(double));
    }

cleanup:
    // --- 6. Free allocated memory ---
    free(dwork);
    if (row_major) {
        free(a_cm);
        free(b_cm);
        free(c_cm);
        free(h_cm);
    }
    return info;
}
