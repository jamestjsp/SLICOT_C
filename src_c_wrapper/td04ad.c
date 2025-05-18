/**
 * @file td04ad.c
 * @brief C wrapper for SLICOT routine TD04AD.
 */

#include <stdlib.h> 
#include <ctype.h>  
#include <stddef.h> 
#include <string.h> 
#include <stdio.h> 

#include "td04ad.h"       
#include "slicot_utils.h" 
#include "slicot_f77.h"   

extern void F77_FUNC(td04ad, TD04AD)(
    const char* rowcol, const int* m, const int* p, const int* index,
    const double* dcoeff, const int* lddcoe,
    double* ucoeff, const int* lduco1, const int* lduco2, 
    int* nr, double* a, const int* lda,
    double* b, const int* ldb, double* c, const int* ldc,
    double* d, const int* ldd, const double* tol,
    int* iwork, double* dwork, const int* ldwork, int* info,
    size_t rowcol_len);

SLICOT_EXPORT
int slicot_td04ad(char rowcol_char, int m_in, int p_in, const int* index_in,
                  const double* dcoeff_c, int lddcoe_c,
                  const double* ucoeff_c, int lduco1_c, int lduco2_c,
                  double tol_in,
                  int* nr_out, double* a_c, int lda_c,
                  double* b_c, int ldb_c, double* c_c, int ldc_c,
                  double* d_c, int ldd_c, int row_major)
{
    int info = 0;
    char rowcol_upper;
    int porm; 
    int n_sum = 0;   
    int kdcoef = 1;  
    int max_idx_val = 0; 

    int *iwork = NULL;
    double *dwork = NULL;
    int liwork_calc = 0;
    int ldwork_calc = 0;
    // ldwork_f will be set by calculation, not query
    int ldwork_f = 0; 

    const double* dcoeff_ptr = NULL; 
    double* ucoeff_ptr = NULL;       
    double* a_ptr = NULL;            
    double* b_ptr = NULL;
    double* c_ptr = NULL;
    double* d_ptr = NULL;

    int lddcoe_f = 1; 
    int lduco1_f = 1;
    int lduco2_f = 1;
    int lda_f = 1;
    int ldb_f = 1;
    int ldc_f = 1;
    int ldd_f = 1;

    double *dcoeff_cm = NULL;
    double *ucoeff_cm = NULL;
    double *a_cm = NULL;
    double *b_cm = NULL;
    double *c_cm = NULL;
    double *d_cm = NULL;

    int m_f = m_in; 
    int p_f = p_in; 

    int dummy_int_array[1] = {0};
    double dummy_double_array[1] = {0.0};


    int max_mp = MAX(m_f, p_f); 
    int max_m_p_1 = MAX(1, max_mp); 

    rowcol_upper = toupper(rowcol_char);
    if (rowcol_upper != 'R' && rowcol_upper != 'C') { info = -1; goto cleanup; }

    if (m_f < 0) { info = -2; goto cleanup; }
    if (p_f < 0) { info = -3; goto cleanup; }

    porm = (rowcol_upper == 'R') ? p_f : m_f;

    if (index_in == NULL && porm > 0) { info = -4; goto cleanup; } 

    if (porm > 0 && index_in != NULL) {
        for (int i = 0; i < porm; ++i) {
            if (index_in[i] < 0) { info = -4; goto cleanup; }
            n_sum += index_in[i];
            if (index_in[i] > max_idx_val) max_idx_val = index_in[i];
        }
        kdcoef = max_idx_val + 1;
    } else { 
        n_sum = 0;
        max_idx_val = 0;
        kdcoef = 1; 
    }
    
    if (p_f == 1 && m_f == 1 && index_in && index_in[0] == -1 && rowcol_upper == 'R') {
         fprintf(stderr, "DEBUG InvalidIndexContent path: porm=%d, kdcoef=%d, dcoeff_c is %sNULL, index_in[0]=%d. Current info=%d (before -5 check)\n",
                 porm, kdcoef, dcoeff_c ? "NOT " : "", index_in[0], info);
    }

    if (dcoeff_c == NULL && porm > 0 && kdcoef > 0) { info = -5; goto cleanup; }
    if (porm > 0 && kdcoef > 0 && dcoeff_c != NULL) { 
        if (row_major) { if (lddcoe_c < MAX(1, kdcoef)) { info = -6; goto cleanup; } } 
        else { if (lddcoe_c < MAX(1, porm)) { info = -6; goto cleanup; } }
    }

    if (ucoeff_c == NULL && p_f > 0 && m_f > 0 && kdcoef > 0) { info = -7; goto cleanup; }
    if (p_f > 0 && m_f > 0 && kdcoef > 0 && ucoeff_c != NULL) {
        if (row_major) {
            if (lduco1_c < MAX(1, m_f)) { info = -8; goto cleanup; }
            if (lduco2_c < MAX(1, kdcoef)) { info = -9; goto cleanup; }
        } else {
            int min_lduco1_f_val = (rowcol_upper == 'R') ? MAX(1,p_f) : MAX(1,max_mp);
            int min_lduco2_f_val = (rowcol_upper == 'R') ? MAX(1,m_f) : MAX(1,max_mp);
            if (lduco1_c < min_lduco1_f_val) { info = -8; goto cleanup; }
            if (lduco2_c < min_lduco2_f_val) { info = -9; goto cleanup; }
        }
    }

    if (nr_out == NULL) { info = -10; goto cleanup; }
    if (a_c == NULL && n_sum > 0) { info = -11; goto cleanup; }
    if (n_sum > 0 && a_c != NULL) { if (lda_c < MAX(1, n_sum)) { info = -12; goto cleanup; } }
    
    if (b_c == NULL && n_sum > 0 && m_f > 0) { info = -13; goto cleanup; }
    if (n_sum > 0 && m_f > 0 && b_c != NULL) { 
        if (row_major) { if (ldb_c < MAX(1, m_f)) { info = -14; goto cleanup; } } 
        else { if (ldb_c < MAX(1, n_sum)) { info = -14; goto cleanup; } }
    }
    
    if (c_c == NULL && p_f > 0 && n_sum > 0) { info = -15; goto cleanup; }
    if (p_f > 0 && n_sum > 0 && c_c != NULL) {
        if (row_major) { if (ldc_c < MAX(1, n_sum)) { info = -16; goto cleanup; } } 
        else { if (ldc_c < MAX(1, p_f)) { info = -16; goto cleanup; } }
    }

    if (d_c == NULL && p_f > 0 && m_f > 0) { info = -17; goto cleanup; }
    if (p_f > 0 && m_f > 0 && d_c != NULL) {
        if (row_major) { if (ldd_c < MAX(1, m_f)) { info = -18; goto cleanup; } } 
        else { if (ldd_c < MAX(1, p_f)) { info = -18; goto cleanup; } }
    }

    if (info != 0) { goto cleanup; }

    // Calculate IWORK size
    liwork_calc = n_sum + max_m_p_1;
    if (liwork_calc == 0) { iwork = NULL; } 
    else { iwork = (int*)malloc((size_t)MAX(1, liwork_calc) * sizeof(int)); CHECK_ALLOC(iwork); }

    // Calculate DWORK size using the formula from Python wrapper / SLICOT docs
    // ldwork = max(1, n_sum + max(n_sum, max(3*m, 3*p)))
    ldwork_calc = MAX(1, n_sum + MAX(n_sum, MAX(3 * m_f, 3 * p_f)));
    dwork = (double*)malloc((size_t)ldwork_calc * sizeof(double));
    CHECK_ALLOC(dwork);
    ldwork_f = ldwork_calc;


    size_t dcoeff_rows_f_dim = (rowcol_upper == 'R') ? MAX(1,p_f) : MAX(1,m_f);
    size_t dcoeff_cols_f_dim = MAX(1,kdcoef);
    size_t dcoeff_size_elems = (porm > 0 && kdcoef > 0) ? (dcoeff_rows_f_dim * dcoeff_cols_f_dim) : 0;

    int ucoeff_f_ld1_dim = (rowcol_upper == 'R') ? MAX(1,p_f) : max_m_p_1;
    int ucoeff_f_ld2_dim = (rowcol_upper == 'R') ? MAX(1,m_f) : max_m_p_1;
    size_t ucoeff_f_depth_dim = MAX(1,kdcoef);
    size_t ucoeff_cm_size_elems = (p_f > 0 && m_f > 0 && kdcoef > 0) ? 
                                  ((size_t)ucoeff_f_ld1_dim * ucoeff_f_ld2_dim * ucoeff_f_depth_dim) : 0;
    if (rowcol_upper == 'C' && (p_f == 0 || m_f == 0)) { 
         if (max_mp > 0 && kdcoef > 0) { 
            ucoeff_cm_size_elems = (size_t)max_m_p_1 * max_m_p_1 * ucoeff_f_depth_dim;
         } else {
            ucoeff_cm_size_elems = 0;
         }
    }

    size_t a_f_dim1_max = MAX(1, n_sum); size_t a_f_dim2_max = MAX(1, n_sum);
    size_t a_size_elems = (n_sum > 0) ? (a_f_dim1_max * a_f_dim2_max) : 0;
    size_t b_f_dim1_max = MAX(1, n_sum); size_t b_f_dim2_max = MAX(1, m_f);
    size_t b_size_elems = (n_sum > 0 && m_f > 0) ? (b_f_dim1_max * b_f_dim2_max) : 0;
    size_t c_f_dim1_max = MAX(1, p_f); size_t c_f_dim2_max = MAX(1, n_sum);
    size_t c_size_elems = (p_f > 0 && n_sum > 0) ? (c_f_dim1_max * c_f_dim2_max) : 0;
    size_t d_f_dim1_max = MAX(1, p_f); size_t d_f_dim2_max = MAX(1, m_f);
    size_t d_size_elems = (p_f > 0 && m_f > 0) ? (d_f_dim1_max * d_f_dim2_max) : 0;

    lddcoe_f = dcoeff_rows_f_dim; lduco1_f = ucoeff_f_ld1_dim; lduco2_f = ucoeff_f_ld2_dim;
    lda_f = a_f_dim1_max; ldb_f = b_f_dim1_max; ldc_f = c_f_dim1_max; ldd_f = d_f_dim1_max;

    if (row_major) {
        if (dcoeff_size_elems > 0 && dcoeff_c != NULL) {
            dcoeff_cm = (double*)malloc(dcoeff_size_elems * sizeof(double)); CHECK_ALLOC(dcoeff_cm);
            slicot_transpose_to_fortran_with_ld(dcoeff_c, dcoeff_cm, porm, kdcoef, lddcoe_c, lddcoe_f, sizeof(double));
            dcoeff_ptr = dcoeff_cm;
        } else { dcoeff_ptr = NULL; }

        if (ucoeff_cm_size_elems > 0 && ucoeff_c != NULL) {
            ucoeff_cm = (double*)malloc(ucoeff_cm_size_elems * sizeof(double)); CHECK_ALLOC(ucoeff_cm);
            memset(ucoeff_cm, 0, ucoeff_cm_size_elems * sizeof(double));
            if (p_f > 0 && m_f > 0) { 
                for (int k_idx = 0; k_idx < kdcoef; ++k_idx) {
                    for (int m_idx_loop = 0; m_idx_loop < m_f; ++m_idx_loop) { 
                        for (int p_idx_loop = 0; p_idx_loop < p_f; ++p_idx_loop) { 
                             ucoeff_cm[p_idx_loop + m_idx_loop*ucoeff_f_ld1_dim + k_idx*ucoeff_f_ld1_dim*ucoeff_f_ld2_dim] =
                                ucoeff_c[p_idx_loop * m_f * kdcoef + m_idx_loop * kdcoef + k_idx];
                        }
                    }
                }
            }
            ucoeff_ptr = ucoeff_cm;
        } else { ucoeff_ptr = NULL; }

        if (a_size_elems > 0 && a_c != NULL) { a_cm = (double*)malloc(a_size_elems * sizeof(double)); CHECK_ALLOC(a_cm); a_ptr = a_cm; } else { a_ptr = NULL; }
        if (b_size_elems > 0 && b_c != NULL) { b_cm = (double*)malloc(b_size_elems * sizeof(double)); CHECK_ALLOC(b_cm); b_ptr = b_cm; } else { b_ptr = NULL; }
        if (c_size_elems > 0 && c_c != NULL) { c_cm = (double*)malloc(c_size_elems * sizeof(double)); CHECK_ALLOC(c_cm); c_ptr = c_cm; } else { c_ptr = NULL; }
        if (d_size_elems > 0 && d_c != NULL) { d_cm = (double*)malloc(d_size_elems * sizeof(double)); CHECK_ALLOC(d_cm); d_ptr = d_cm; } else { d_ptr = NULL; }
    } else { 
        dcoeff_ptr = (dcoeff_size_elems > 0 && dcoeff_c != NULL) ? dcoeff_c : NULL;
        ucoeff_ptr = (ucoeff_cm_size_elems > 0 && ucoeff_c != NULL) ? (double*)ucoeff_c : NULL; 
        a_ptr = (a_size_elems > 0 && a_c != NULL) ? a_c : NULL;
        b_ptr = (b_size_elems > 0 && b_c != NULL) ? b_c : NULL;
        c_ptr = (c_size_elems > 0 && c_c != NULL) ? c_c : NULL;
        d_ptr = (d_size_elems > 0 && d_c != NULL) ? d_c : NULL;
    }
    
    const int* final_index_ptr = (porm > 0 && index_in != NULL) ? index_in : dummy_int_array;
    
    F77_FUNC(td04ad, TD04AD)(&rowcol_upper, &m_f, &p_f, final_index_ptr,
                             dcoeff_ptr ? dcoeff_ptr : dummy_double_array, &lddcoe_f,
                             ucoeff_ptr ? ucoeff_ptr : dummy_double_array, &lduco1_f, &lduco2_f,
                             nr_out, 
                             a_ptr ? a_ptr : dummy_double_array, &lda_f, 
                             b_ptr ? b_ptr : dummy_double_array, &ldb_f,
                             c_ptr ? c_ptr : dummy_double_array, &ldc_f, 
                             d_ptr ? d_ptr : dummy_double_array, &ldd_f,
                             &tol_in, iwork, dwork, &ldwork_f, &info, (size_t)1);

    if (row_major && info == 0) {
        if (a_ptr && nr_out && *nr_out > 0 && a_c != NULL) {
            slicot_transpose_to_c_with_ld(a_ptr, a_c, *nr_out, *nr_out, lda_f, lda_c, sizeof(double));
        }
        if (b_ptr && nr_out && *nr_out > 0 && m_f > 0 && b_c != NULL) {
            slicot_transpose_to_c_with_ld(b_ptr, b_c, *nr_out, m_f, ldb_f, ldb_c, sizeof(double));
        }
        if (c_ptr && p_f > 0 && nr_out && *nr_out > 0 && c_c != NULL) {
            slicot_transpose_to_c_with_ld(c_ptr, c_c, p_f, *nr_out, ldc_f, ldc_c, sizeof(double));
        }
        if (d_ptr && p_f > 0 && m_f > 0 && d_c != NULL) {
            slicot_transpose_to_c_with_ld(d_ptr, d_c, p_f, m_f, ldd_f, ldd_c, sizeof(double));
        }
        if (rowcol_upper == 'C' && ucoeff_ptr && p_f > 0 && m_f > 0 && kdcoef > 0 && ucoeff_c != NULL) {
            for (int k_idx = 0; k_idx < kdcoef; ++k_idx) {
                for (int m_idx_loop = 0; m_idx_loop < m_f; ++m_idx_loop) {
                    for (int p_idx_loop = 0; p_idx_loop < p_f; ++p_idx_loop) {
                        ((double*)ucoeff_c)[p_idx_loop * m_f * kdcoef + m_idx_loop * kdcoef + k_idx] =
                            ucoeff_ptr[p_idx_loop + m_idx_loop*ucoeff_f_ld1_dim + k_idx*ucoeff_f_ld1_dim*ucoeff_f_ld2_dim];
                    }
                }
            }
        }
    }

cleanup:
    free(iwork); free(dwork);
    if (row_major) {
        free(dcoeff_cm); free(ucoeff_cm);
        free(a_cm); free(b_cm); free(c_cm); free(d_cm);
    }
    return info;
}
