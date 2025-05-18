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
int slicot_tb04ad(const char* rowcol_in, int n_c, int m_c, int p_c,
                  double* a_c, int lda_c,
                  double* b_c, int ldb_c,
                  double* c_c, int ldc_c,
                  const double* d_c_in, int ldd_c,
                  int* nr_out, int* index_out, double* dcoeff_out, int lddcoe_c,
                  double* ucoeff_out, int lduco1_c, int lduco2_c,
                  double tol1_in, double tol2_in,
                  int row_major_flag)
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

    char rowcol_f;
    int porm_f, porp_f; // Dimensions based on ROWCOL for Fortran outputs INDEX, DCOEFF, UCOEFF

    // Initial check for rowcol_in
    if (rowcol_in == NULL || (rowcol_in[0] != 'R' && rowcol_in[0] != 'r' && rowcol_in[0] != 'C' && rowcol_in[0] != 'c')) {
        info = -1; goto cleanup;
    }
    rowcol_f = toupper(rowcol_in[0]);

    if (rowcol_f == 'R') {
        porm_f = p_c; 
        porp_f = m_c; 
    } else { // ROWCOL = 'C'
        porm_f = m_c; 
        porp_f = p_c; 
    }
    
    // 2. Input parameter validation
    if (n_c < 0) { info = -2; goto cleanup; }
    if (m_c < 0) { info = -3; goto cleanup; }
    if (p_c < 0) { info = -4; goto cleanup; }

    if (a_c == NULL && n_c > 0) { info = -5; goto cleanup; }
    if (b_c == NULL && n_c > 0 && m_c > 0) { info = -7; goto cleanup; } 
    if (c_c == NULL && p_c > 0 && n_c > 0) { info = -9; goto cleanup; } 
    if (d_c_in == NULL && p_c > 0 && m_c > 0) { info = -11; goto cleanup; } 
    if (nr_out == NULL) { info = -12; goto cleanup; }
    
    if (index_out == NULL && porm_f > 0) { info = -13; goto cleanup; }
    if (dcoeff_out == NULL && porm_f > 0) { info = -14; goto cleanup; }
    if (ucoeff_out == NULL && porm_f > 0 && porp_f > 0) { info = -16; goto cleanup; }

    // Fortran-side leading dimensions (number of rows for Fortran)
    int lda_f_expected = MAX(1, n_c);
    int ldb_f_expected = MAX(1, n_c); 
    int ldc_f_expected = (rowcol_f == 'R') ? MAX(1,p_c) : MAX(1, MAX(m_c,p_c));
    int ldd_f_expected = (rowcol_f == 'R') ? MAX(1,p_c) : MAX(1, MAX(m_c,p_c));

    if (row_major_flag) {
        // C LDs are number of columns
        if (n_c > 0 && lda_c < n_c) { info = -6; goto cleanup; } 
        if (n_c > 0 && m_c > 0 && ldb_c < m_c) { info = -8; goto cleanup; } 
        if (p_c > 0 && n_c > 0 && ldc_c < n_c) { info = -10; goto cleanup; } 
        if (p_c > 0 && m_c > 0 && ldd_c < m_c) { info = -11; goto cleanup; } 
    } else { // Column Major C
        if (n_c > 0 && lda_c < lda_f_expected) { info = -6; goto cleanup; }
        if (n_c > 0 && m_c > 0 && ldb_c < ldb_f_expected) { info = -8; goto cleanup; } 
        if (n_c > 0 && (p_c > 0 || (rowcol_f == 'C' && (m_c > 0 || p_c > 0))) && ldc_c < ldc_f_expected ) { info = -10; goto cleanup; }
        if ( (p_c > 0 && m_c > 0) || (rowcol_f == 'C' && (m_c > 0 || p_c > 0)) ) {
             if (ldd_c < ldd_f_expected) {info = -11; goto cleanup;}
        }
    }

    if (porm_f > 0 && lddcoe_c < porm_f) { info = -15; goto cleanup; }
    if (porm_f > 0 && lduco1_c < porm_f) { info = -17; goto cleanup; }
    if (porp_f > 0 && lduco2_c < porp_f) { info = -18; goto cleanup; }

    if (info != 0) { goto cleanup; }

    // 3. Internal Workspace Allocation (Formula based)
    int mp_ws = (rowcol_f == 'R') ? m_c : p_c;
    int pm_ws = (rowcol_f == 'R') ? p_c : m_c;
    ldwork = MAX(1, n_c*(n_c+1) + MAX(n_c*mp_ws + 2*n_c + MAX(n_c,mp_ws), MAX(3*mp_ws, pm_ws)));
    if (n_c==0 && m_c==0 && p_c==0 && ldwork < 1) ldwork=1; 

    liwork = n_c + MAX(m_c, p_c);
    liwork = MAX(1, liwork);

    dwork = (double*)malloc((size_t)ldwork * sizeof(double));
    CHECK_ALLOC(dwork);
    iwork = (int*)malloc((size_t)liwork * sizeof(int));
    CHECK_ALLOC(iwork);

    // 4. Prepare Fortran pointers and temporary column-major buffers if needed
    double* a_f_ptr = a_c; // Default to C pointer
    double* b_f_ptr = b_c;
    double* c_f_ptr = c_c;
    double* d_f_ptr = (double*)d_c_in; // Fortran D is not const if ROWCOL='C'

    // Fortran buffer element counts (rows_f * cols_f)
    // These are the total elements Fortran might access in its view of the array
    size_t a_f_buf_elems = (n_c > 0) ? (size_t)lda_f_expected * n_c : 0;
    size_t b_f_buf_elems = (n_c > 0) ? (size_t)ldb_f_expected * MAX(1, (rowcol_f == 'R') ? m_c : MAX(m_c,p_c)) : 0;
    size_t c_f_buf_elems = (n_c > 0) ? (size_t)ldc_f_expected * n_c : 0;
    size_t d_f_buf_elems = (p_c > 0 || m_c > 0 || rowcol_f == 'C') ? /*D is PxM or MAX(M,P)xMAX(M,P) workspace*/
                           (size_t)ldd_f_expected * MAX(1, (rowcol_f == 'R') ? m_c : MAX(m_c,p_c)) : 0;


    if (row_major_flag) {
        // For row_major, Fortran LDs are the number of rows in the conceptual Fortran matrix
        // lda_f_expected, ldb_f_expected, ldc_f_expected, ldd_f_expected are already set.

        if (n_c > 0 && a_c != NULL) {
            a_cm = (double*)malloc(a_f_buf_elems * sizeof(double)); CHECK_ALLOC(a_cm);
            slicot_transpose_to_fortran_with_ld(a_c, a_cm, n_c, n_c, lda_c, lda_f_expected, sizeof(double));
            a_f_ptr = a_cm;
        } else { a_f_ptr = NULL; }

        if (n_c > 0 && b_c != NULL) { 
            b_cm = (double*)malloc(b_f_buf_elems * sizeof(double)); CHECK_ALLOC(b_cm);
            memset(b_cm, 0, b_f_buf_elems * sizeof(double)); 
            if (m_c > 0) { 
                slicot_transpose_to_fortran_with_ld(b_c, b_cm, n_c, m_c, ldb_c, ldb_f_expected, sizeof(double));
            }
            b_f_ptr = b_cm;
        } else { b_f_ptr = NULL; }
        
        if (n_c > 0 && c_c != NULL) { 
            c_cm = (double*)malloc(c_f_buf_elems * sizeof(double)); CHECK_ALLOC(c_cm);
            memset(c_cm, 0, c_f_buf_elems * sizeof(double)); 
            if (p_c > 0) { 
                slicot_transpose_to_fortran_with_ld(c_c, c_cm, p_c, n_c, ldc_c, ldc_f_expected, sizeof(double));
            }
            c_f_ptr = c_cm;
        } else { c_f_ptr = NULL; }

        if (d_f_buf_elems > 0 ) { 
            d_cm = (double*)malloc(d_f_buf_elems * sizeof(double)); CHECK_ALLOC(d_cm);
            memset(d_cm, 0, d_f_buf_elems * sizeof(double));
            if (d_c_in != NULL && p_c > 0 && m_c > 0) {
                if (rowcol_f == 'R') {
                    slicot_transpose_to_fortran_with_ld(d_c_in, d_cm, p_c, m_c, ldd_c, ldd_f_expected, sizeof(double));
                } else { 
                    double* temp_d_pxm_cm = (double*)malloc((size_t)p_c * m_c * sizeof(double));
                    CHECK_ALLOC(temp_d_pxm_cm);
                    slicot_transpose_to_fortran_with_ld(d_c_in, temp_d_pxm_cm, p_c, m_c, ldd_c, MAX(1,p_c), sizeof(double));
                    for (int col_j = 0; col_j < m_c; ++col_j) {
                        for (int row_i = 0; row_i < p_c; ++row_i) {
                            d_cm[row_i + col_j * ldd_f_expected] = temp_d_pxm_cm[row_i + col_j * MAX(1,p_c)];
                        }
                    }
                    free(temp_d_pxm_cm);
                }
            }
            d_f_ptr = d_cm;
        } else {
            d_f_ptr = NULL;
        }

    } else { // Column Major C
        if (rowcol_f == 'C' && d_f_buf_elems > 0) {
            d_cm = (double*)malloc(d_f_buf_elems * sizeof(double)); CHECK_ALLOC(d_cm);
            memset(d_cm, 0, d_f_buf_elems * sizeof(double));
            if (d_c_in != NULL && p_c > 0 && m_c > 0) { 
                for (int col_j = 0; col_j < m_c; ++col_j) {
                    for (int row_i = 0; row_i < p_c; ++row_i) {
                        d_cm[row_i + col_j * ldd_f_expected] = d_c_in[row_i + col_j * ldd_c];
                    }
                }
            }
            d_f_ptr = d_cm;
        } else if (d_f_buf_elems > 0 && d_c_in != NULL){ 
             d_cm = (double*)malloc(d_f_buf_elems * sizeof(double)); CHECK_ALLOC(d_cm);
             if (p_c > 0 && m_c > 0) { 
                for(int j=0; j<m_c; ++j) { 
                    for(int i=0; i<p_c; ++i) {
                        d_cm[i + j*ldd_f_expected] = d_c_in[i + j*ldd_c];
                    }
                }
             }
             d_f_ptr = d_cm;
        } else { 
            d_f_ptr = NULL;
        }
        // Ensure NULL is passed if logical size is 0 and input pointer is NULL
        if (n_c == 0 && a_c == NULL) a_f_ptr = NULL;
        if ((n_c == 0 || m_c == 0) && b_c == NULL) b_f_ptr = NULL;
        if ((p_c == 0 || n_c == 0) && c_c == NULL) c_f_ptr = NULL;
    }

    // Fortran leading dimensions for output coefficient arrays
    int f_lddcoe_fortran = MAX(1, porm_f);
    int f_lduco1_fortran = MAX(1, porm_f);
    int f_lduco2_fortran = MAX(1, porp_f);
    
    int* index_f_pass = (porm_f > 0 ? index_out : NULL);
    double* dcoeff_f_pass = (porm_f > 0 ? dcoeff_out : NULL);
    double* ucoeff_f_pass = (porm_f > 0 && porp_f > 0 ? ucoeff_out : NULL);

    // 7. Call Fortran function
    F77_FUNC(tb04ad, TB04AD)(&rowcol_f, &n_c, &m_c, &p_c,
                             a_f_ptr, &lda_f_expected, 
                             b_f_ptr, &ldb_f_expected, 
                             c_f_ptr, &ldc_f_expected, 
                             d_f_ptr, &ldd_f_expected,
                             nr_out, index_f_pass, dcoeff_f_pass, &f_lddcoe_fortran,
                             ucoeff_f_pass, &f_lduco1_fortran, &f_lduco2_fortran,
                             &tol1_in, &tol2_in,
                             iwork, dwork, &ldwork, &info, 1);

    // 8. Convert results back to row-major
    if (row_major_flag && info == 0) {
        if (n_c > 0 && a_cm != NULL && a_c != NULL) { 
            slicot_transpose_to_c_with_ld(a_cm, a_c, n_c, n_c, lda_f_expected, lda_c, sizeof(double));
        }
        if (*nr_out > 0 && m_c > 0 && b_cm != NULL && b_c != NULL) {
            slicot_transpose_to_c_with_ld(b_cm, b_c, *nr_out, m_c, ldb_f_expected, ldb_c, sizeof(double));
        }
        if (p_c > 0 && *nr_out > 0 && c_cm != NULL && c_c != NULL) {
            slicot_transpose_to_c_with_ld(c_cm, c_c, p_c, *nr_out, ldc_f_expected, ldc_c, sizeof(double));
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
