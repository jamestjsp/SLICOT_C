/**
 * @file tb04ad.c
 * @brief C wrapper for SLICOT routine TB04AD.
 * @details Converts a state-space system (A,B,C,D) to a transfer matrix T(s).
 * Workspace (IWORK, DWORK) is allocated internally by this wrapper.
 * Input/output matrix format is handled via the row_major parameter.
 */

#include <stdlib.h>
#include <ctype.h>
#include <stddef.h>
#include <math.h>
#include <stdio.h> // For error logging (optional)
#include <string.h> // For memcpy

#include "tb04ad.h"       // Public header for this wrapper
#include "slicot_utils.h" // Provides CHECK_ALLOC, SLICOT_MEMORY_ERROR, MAX/MIN, transpose functions etc.
#include "slicot_f77.h"   // Provides F77_FUNC macro for Fortran name mangling

/* External Fortran routine declaration */
extern void F77_FUNC(tb04ad, TB04AD)(
    const char* rowcol, const int* n, const int* m, const int* p,
    double* a, const int* lda,
    double* b, const int* ldb,
    double* c, const int* ldc,
    double* d, const int* ldd, /* D is modified internally by Fortran if ROWCOL='C' */
    int* nr, int* index, double* dcoeff, const int* lddcoe,
    double* ucoeff, const int* lduco1, const int* lduco2,
    const double* tol1, const double* tol2,
    int* iwork, double* dwork, const int* ldwork,
    int* info,
    int rowcol_len);

/* C wrapper function definition */
SLICOT_EXPORT
int slicot_tb04ad(const char* rowcol_in, int n, int m, int p,
                  double* a, int lda_c,
                  double* b, int ldb_c,
                  double* c, int ldc_c,
                  const double* d_in, int ldd_c,
                  int* nr, int* index, double* dcoeff_out, int lddcoe_c,
                  double* ucoeff_out, int lduco1_c, int lduco2_c,
                  double tol1, double tol2,
                  int row_major)
{
    // 1. Variable declarations
    int info = 0;
    int *iwork = NULL;
    double *dwork = NULL;
    int liwork = 0;
    int ldwork = 0;

    double* a_cm = NULL;
    double* b_cm = NULL;
    double* c_cm = NULL;
    double* d_cm = NULL;

    char rowcol_upper_char;
    int porm, porp;

    // Initial check for rowcol_in
    if (rowcol_in == NULL || (rowcol_in[0] != 'R' && rowcol_in[0] != 'r' && rowcol_in[0] != 'C' && rowcol_in[0] != 'c')) {
        info = -1; goto cleanup;
    }
    rowcol_upper_char = toupper(rowcol_in[0]);

    if (rowcol_upper_char == 'R') {
        porm = p;
        porp = m;
    } else {
        porm = m;
        porp = p;
    }
    
    // 2. Input parameter validation
    if (n < 0) { info = -2; goto cleanup; }
    if (m < 0) { info = -3; goto cleanup; }
    if (p < 0) { info = -4; goto cleanup; }

    if (a == NULL && n > 0) { info = -5; goto cleanup; }
    if (b == NULL && n > 0 && m > 0) { info = -7; goto cleanup; }
    if (c == NULL && p > 0 && n > 0) { info = -9; goto cleanup; }
    if (d_in == NULL && p > 0 && m > 0) { info = -11; goto cleanup; }
    if (nr == NULL) { info = -12; goto cleanup; }
    if (index == NULL && porm > 0) { info = -13; goto cleanup; }
    if (dcoeff_out == NULL && porm > 0) { info = -14; goto cleanup; }
    if (ucoeff_out == NULL && porm > 0 && porp > 0) { info = -16; goto cleanup; }

    int min_lda_f = MAX(1, n);
    int min_ldb_f = MAX(1, n);
    int min_ldc_f_val = (rowcol_upper_char == 'R') ? MAX(1,p) : MAX(1, MAX(m,p));
    int min_ldd_f_val = (rowcol_upper_char == 'R') ? MAX(1,p) : MAX(1, MAX(m,p));

    if (row_major) {
        if (n > 0 && lda_c < n) { info = -6; goto cleanup; }
        if (n > 0 && m > 0 && ldb_c < m) { info = -8; goto cleanup; }
        if (p > 0 && n > 0 && ldc_c < n) { info = -10; goto cleanup; }
        if (p > 0 && m > 0 && ldd_c < m) { info = -11; goto cleanup; }
    } else { 
        if (n > 0 && lda_c < min_lda_f) { info = -6; goto cleanup; }
        if (n > 0 && m > 0 && ldb_c < min_ldb_f) { info = -8; goto cleanup; }
        if ( (p > 0 || (rowcol_upper_char == 'C' && m > 0)) && n > 0 && ldc_c < min_ldc_f_val ) { info = -10; goto cleanup; }
        if ( (p > 0 && m > 0) || (rowcol_upper_char == 'C' && (m>0 || p>0)) ) {
             if (ldd_c < min_ldd_f_val) {info = -11; goto cleanup;}
        }
    }

    if (porm > 0 && lddcoe_c < porm) { info = -15; goto cleanup; }
    if (porm > 0 && lduco1_c < porm) { info = -17; goto cleanup; }
    if (porp > 0 && lduco2_c < porp) { info = -18; goto cleanup; }

    if (info != 0) { goto cleanup; }

    // 3. Internal Workspace Allocation (Method B: Formula based)
    int mp_ws = (rowcol_upper_char == 'R') ? m : p;
    int pm_ws = (rowcol_upper_char == 'R') ? p : m;
    ldwork = MAX(1, n*(n+1) + MAX(n*mp_ws + 2*n + MAX(n,mp_ws), MAX(3*mp_ws, pm_ws)));
    if (n==0 && m==0 && p==0 && ldwork < 1) ldwork=1; // Ensure ldwork >= 1 for all-zero case

    liwork = n + MAX(m, p);
    liwork = MAX(1, liwork);

    dwork = (double*)malloc((size_t)ldwork * sizeof(double));
    CHECK_ALLOC(dwork);
    iwork = (int*)malloc((size_t)liwork * sizeof(int));
    CHECK_ALLOC(iwork);

    // 4. Memory allocation for column-major copies and Fortran pointers
    size_t a_size_bytes = (n > 0) ? (size_t)n * n * sizeof(double) : 0;
    size_t b_size_bytes = (n > 0 && m > 0) ? (size_t)n * m * sizeof(double) : 0;
    
    int lda_f = min_lda_f;
    int ldb_f = min_ldb_f;
    int ldc_f = min_ldc_f_val;
    int ldd_f = min_ldd_f_val;

    double* a_ptr = (n > 0 && a != NULL) ? a : NULL;
    double* b_ptr = (n > 0 && m > 0 && b != NULL) ? b : NULL;
    double* c_ptr = ( (p > 0 || (rowcol_upper_char == 'C' && m > 0)) && n > 0 && c != NULL) ? c : NULL;
    double* d_ptr_fortran = NULL; // To be determined

    size_t c_fortran_buf_size = 0;
     if ( (p > 0 || (rowcol_upper_char == 'C' && m > 0)) && n > 0 ) {
        c_fortran_buf_size = (size_t)ldc_f * n * sizeof(double);
     }

    size_t d_fortran_buf_size = 0;
    if (rowcol_upper_char == 'R' && p > 0 && m > 0) {
        d_fortran_buf_size = (size_t)ldd_f * m * sizeof(double);
    } else if (rowcol_upper_char == 'C' && (m > 0 || p > 0)) {
        d_fortran_buf_size = (size_t)ldd_f * MAX(1,MAX(m,p)) * sizeof(double);
    }

    if (row_major) {
        lda_f = MAX(1, n); 
        ldb_f = MAX(1, n); 

        if (a_size_bytes > 0 && a != NULL) {
            a_cm = (double*)malloc(a_size_bytes); CHECK_ALLOC(a_cm);
            slicot_transpose_to_fortran_with_ld(a, a_cm, n, n, lda_c, lda_f, sizeof(double));
            a_ptr = a_cm;
        }

        if (b_size_bytes > 0 && b != NULL) {
            b_cm = (double*)malloc(b_size_bytes); CHECK_ALLOC(b_cm);
            slicot_transpose_to_fortran_with_ld(b, b_cm, n, m, ldb_c, ldb_f, sizeof(double));
            b_ptr = b_cm;
        }
        
        if (c_fortran_buf_size > 0 && c != NULL) {
            c_cm = (double*)malloc(c_fortran_buf_size); CHECK_ALLOC(c_cm);
            memset(c_cm, 0, c_fortran_buf_size);
            if (p > 0 && n > 0) { 
                 slicot_transpose_to_fortran_with_ld(c, c_cm, p, n, ldc_c, ldc_f, sizeof(double));
            }
            c_ptr = c_cm;
        }

        if (d_fortran_buf_size > 0) {
            d_cm = (double*)malloc(d_fortran_buf_size); CHECK_ALLOC(d_cm);
            memset(d_cm, 0, d_fortran_buf_size);
            if (d_in != NULL && p > 0 && m > 0) {
                 if (rowcol_upper_char == 'R') {
                    slicot_transpose_to_fortran_with_ld(d_in, d_cm, p, m, ldd_c, ldd_f, sizeof(double));
                } else { 
                    double* temp_d_pxm_cm = (double*)malloc((size_t)p * m * sizeof(double));
                    CHECK_ALLOC(temp_d_pxm_cm);
                    slicot_transpose_to_fortran_with_ld(d_in, temp_d_pxm_cm, p, m, ldd_c, MAX(1,p), sizeof(double));
                    for (int col_j = 0; col_j < m; ++col_j) {
                        for (int row_i = 0; row_i < p; ++row_i) {
                            d_cm[row_i + col_j * ldd_f] = temp_d_pxm_cm[row_i + col_j * MAX(1,p)];
                        }
                    }
                    free(temp_d_pxm_cm);
                }
            }
            d_ptr_fortran = d_cm;
        } else {
            d_ptr_fortran = NULL; // If d_in is NULL or D has zero logical size for 'R'
        }

    } else { // Column Major C
        if (rowcol_upper_char == 'C' && d_fortran_buf_size > 0) {
            d_cm = (double*)malloc(d_fortran_buf_size); CHECK_ALLOC(d_cm);
            memset(d_cm, 0, d_fortran_buf_size);
            if (d_in != NULL && p > 0 && m > 0) { 
                for (int col_j = 0; col_j < m; ++col_j) {
                    for (int row_i = 0; row_i < p; ++row_i) {
                        d_cm[row_i + col_j * ldd_f] = d_in[row_i + col_j * ldd_c];
                    }
                }
            }
            d_ptr_fortran = d_cm;
        } else if (d_fortran_buf_size > 0 && d_in != NULL){ // ROWCOL='R', make a copy
             d_cm = (double*)malloc(d_fortran_buf_size); CHECK_ALLOC(d_cm);
             if (p > 0 && m > 0) {
                for(int j=0; j<m; ++j) { 
                    for(int i=0; i<p; ++i) {
                        d_cm[i + j*ldd_f] = d_in[i + j*ldd_c];
                    }
                }
             }
             d_ptr_fortran = d_cm;
        } else { // d_in is NULL or D has zero logical size
            d_ptr_fortran = NULL;
        }
    }

    int f_lddcoe = MAX(1, porm);
    int f_lduco1 = MAX(1, porm);
    int f_lduco2 = MAX(1, porp);
    
    int* index_pass = (porm > 0 ? index : NULL);
    double* dcoeff_pass = (porm > 0 ? dcoeff_out : NULL);
    double* ucoeff_pass = (porm > 0 && porp > 0 ? ucoeff_out : NULL);
    if (porm == 0 && dcoeff_out != NULL) dcoeff_pass = dcoeff_out; // Allow dummy if porm=0 but user provided buffer
    if ((porm == 0 || porp == 0) && ucoeff_out != NULL) ucoeff_pass = ucoeff_out;


    // 7. Call Fortran function
    F77_FUNC(tb04ad, TB04AD)(&rowcol_upper_char, &n, &m, &p,
                             a_ptr, &lda_f, b_ptr, &ldb_f, c_ptr, &ldc_f, d_ptr_fortran, &ldd_f,
                             nr, index_pass, dcoeff_pass, &f_lddcoe,
                             ucoeff_pass, &f_lduco1, &f_lduco2,
                             &tol1, &tol2,
                             iwork, dwork, &ldwork, &info, 1);

    // 8. Convert results back to row-major
    if (row_major && info == 0) {
        if (n > 0 && a_cm != NULL && a != NULL) {
            slicot_transpose_to_c_with_ld(a_cm, a, n, n, lda_f, lda_c, sizeof(double));
        }
        if (*nr > 0 && m > 0 && b_cm != NULL && b != NULL) {
            slicot_transpose_to_c_with_ld(b_cm, b, *nr, m, ldb_f, ldb_c, sizeof(double));
        }
        if (p > 0 && *nr > 0 && c_cm != NULL && c != NULL) {
            slicot_transpose_to_c_with_ld(c_cm, c, p, *nr, ldc_f, ldc_c, sizeof(double));
        }
    }

cleanup:
    free(iwork);
    free(dwork);
    free(a_cm);
    free(b_cm);
    free(c_cm);
    free(d_cm);

    return info;
}