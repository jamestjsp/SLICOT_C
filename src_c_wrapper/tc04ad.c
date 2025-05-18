/**
 * @file tc04ad.c
 * @brief C wrapper for SLICOT routine TC04AD.
 * @details Converts a polynomial matrix representation to state-space.
 */

#include <stdlib.h>
#include <ctype.h>
#include <string.h> // For memcpy
#include <stdio.h>  // For debugging

#include "tc04ad.h"
#include "slicot_utils.h"
#include "slicot_f77.h"

// Fortran routine declaration
extern void F77_FUNC(tc04ad, TC04AD)(
    const char* leri, const int* m, const int* p, const int* index,
    double* pcoeff, const int* ldpco1, const int* ldpco2,
    double* qcoeff, const int* ldqco1, const int* ldqco2,
    int* n, double* rcond,
    double* a, const int* lda,
    double* b, const int* ldb,
    double* c, const int* ldc,
    double* d, const int* ldd,
    int* iwork, double* dwork, const int* ldwork,
    int* info,
    int leri_len);

SLICOT_EXPORT
int slicot_tc04ad(char leri_c_in, int m_c, int p_c, const int* index_c,
                  const double* pcoeff_c_in, int ldpcoeff_c_rows, int ldpcoeff_c_cols,
                  const double* qcoeff_c_in, int ldqcoeff_c_rows, int ldqcoeff_c_cols,
                  int* n_out, double* rcond_out,
                  double* a_out, int lda_c,
                  double* b_out, int ldb_c,
                  double* c_out, int ldc_c,
                  double* d_out, int ldd_c,
                  int row_major_flag) {
    int info = 0;
    char leri_f = toupper(leri_c_in);

    // Workspace arrays
    int *iwork = NULL;
    double *dwork = NULL;
    int ldwork = 0;
    int liwork = 0; 

    // Column-major temporary arrays for Fortran
    double *pcoeff_cm = NULL;
    double *qcoeff_cm = NULL;
    double *a_cm = NULL;
    double *b_cm = NULL;
    double *c_cm = NULL;
    double *d_cm = NULL;
    
    // --- Parameter Validation ---
    if (leri_f != 'L' && leri_f != 'R') { info = -1; goto cleanup; }
    if (m_c < 0) { info = -2; goto cleanup; }
    if (p_c < 0) { info = -3; goto cleanup; }
    if (index_c == NULL && ( (leri_f == 'L' && p_c > 0) || (leri_f == 'R' && m_c > 0) ) ) { info = -4; goto cleanup; }
    
    int N_val = 0;
    int max_deg = 0;
    int index_dim = (leri_f == 'L') ? p_c : m_c;
    if (index_dim > 0 && index_c != NULL) {
        for (int i = 0; i < index_dim; ++i) {
            if (index_c[i] < 0) { info = -4; goto cleanup; } 
            N_val += index_c[i];
            max_deg = MAX(max_deg, index_c[i]);
        }
    } else if (index_dim > 0 && index_c == NULL) { 
         info = -4; goto cleanup;
    }
    int kpcoef = max_deg + 1;
    if (index_dim == 0 && (p_c > 0 || m_c > 0) ) { 
         kpcoef = 1;
    } else if (index_dim == 0 && p_c == 0 && m_c == 0) { 
        kpcoef = 1; 
    }

    int porm_f = (leri_f == 'L') ? p_c : m_c;
    int porp_f = (leri_f == 'L') ? m_c : p_c;

    if (pcoeff_c_in == NULL && porm_f > 0 && kpcoef > 0) { info = -5; goto cleanup; }
    if (qcoeff_c_in == NULL && porm_f > 0 && porp_f > 0 && kpcoef > 0) { info = -8; goto cleanup; }
    
    if (n_out == NULL) { info = -11; goto cleanup; } 
    if (rcond_out == NULL) { info = -12; goto cleanup; }
    if (a_out == NULL && N_val > 0) { info = -13; goto cleanup; }
    if (b_out == NULL && N_val > 0 && m_c > 0) { info = -15; goto cleanup; }
    if (c_out == NULL && p_c > 0 && N_val > 0) { info = -17; goto cleanup; }
    if (d_out == NULL && p_c > 0 && m_c > 0) { info = -19; goto cleanup; }

    int ldpco1_f_exp = MAX(1, porm_f);
    int ldpco2_f_exp = MAX(1, porm_f);
    int ldqco1_f_exp = (leri_f == 'L') ? MAX(1, p_c) : MAX(1, MAX(m_c, p_c));
    int ldqco2_f_exp = (leri_f == 'L') ? MAX(1, m_c) : MAX(1, MAX(m_c, p_c));
    int lda_f_exp = MAX(1, N_val);
    int ldb_f_exp = MAX(1, N_val); 
    int ldc_f_exp = MAX(1, MAX(m_c, p_c)); 
    int ldd_f_exp = MAX(1, MAX(m_c, p_c)); 

    if (row_major_flag) {
        if (porm_f > 0 && kpcoef > 0 && pcoeff_c_in != NULL && (ldpcoeff_c_rows < porm_f || ldpcoeff_c_cols < porm_f)) { info = -6; goto cleanup; }
        if (qcoeff_c_in != NULL ) { 
            int q_c_rows_logical = (leri_f == 'L') ? p_c : m_c; 
            int q_c_cols_logical = (leri_f == 'L') ? m_c : p_c;
            if (leri_f == 'R') { 
                q_c_rows_logical = MAX(m_c,p_c); q_c_cols_logical = MAX(m_c,p_c);
            }
            if ( (porm_f > 0 && porp_f > 0 && kpcoef > 0) && 
                 (ldqcoeff_c_rows < q_c_rows_logical || ldqcoeff_c_cols < q_c_cols_logical ) ) { info = -9; goto cleanup; }
        }
        if (N_val > 0 && a_out != NULL && lda_c < N_val) { info = -14; goto cleanup; } 
        if (N_val > 0 && m_c > 0 && b_out != NULL && ldb_c < m_c) { info = -16; goto cleanup; } 
        if (p_c > 0 && N_val > 0 && c_out != NULL && ldc_c < N_val) { info = -18; goto cleanup; } 
        if (p_c > 0 && m_c > 0 && d_out != NULL && ldd_c < m_c) { info = -20; goto cleanup; } 
    } else { // Column-major C
        if (porm_f > 0 && kpcoef > 0 && pcoeff_c_in != NULL && (ldpcoeff_c_rows < ldpco1_f_exp || ldpcoeff_c_cols < ldpco2_f_exp)) { info = -6; goto cleanup; }
        if (porm_f > 0 && porp_f > 0 && kpcoef > 0 && qcoeff_c_in != NULL && (ldqcoeff_c_rows < ldqco1_f_exp || ldqcoeff_c_cols < ldqco2_f_exp)) { info = -9; goto cleanup; }
        if (N_val > 0 && a_out != NULL && lda_c < lda_f_exp) { info = -14; goto cleanup; }
        if (N_val > 0 && b_out != NULL && ldb_c < ldb_f_exp) { info = -16; goto cleanup; } 
        if (N_val > 0 && c_out != NULL && ldc_c < ldc_f_exp) { info = -18; goto cleanup; } 
        if (d_out != NULL && (p_c > 0 || m_c > 0 || MAX(m_c, p_c) > 0) && ldd_c < ldd_f_exp) { info = -20; goto cleanup; } 
    }
    if (info != 0) goto cleanup;

    // --- Workspace Allocation (Formula based) ---
    int max_mp = MAX(1, MAX(m_c, p_c)); 
    ldwork = MAX(1, max_mp * (max_mp + 4));
    liwork = 2 * max_mp;
    liwork = MAX(1, liwork); 

    iwork = (int*)malloc((size_t)liwork * sizeof(int)); CHECK_ALLOC(iwork);
    dwork = (double*)malloc((size_t)ldwork * sizeof(double)); CHECK_ALLOC(dwork);

    // --- Prepare Fortran pointers and temporary CM arrays ---
    double* pcoeff_f_ptr = (double*)pcoeff_c_in; 
    double* qcoeff_f_ptr = (double*)qcoeff_c_in; 
    double* a_f_ptr = a_out;
    double* b_f_ptr = b_out;
    double* c_f_ptr = c_out;
    double* d_f_ptr = d_out;

    size_t pcoeff_f_elems = (porm_f > 0 && kpcoef > 0) ? (size_t)ldpco1_f_exp * ldpco2_f_exp * kpcoef : 0;
    size_t qcoeff_f_elems = (porm_f > 0 && porp_f > 0 && kpcoef > 0) ? (size_t)ldqco1_f_exp * ldqco2_f_exp * kpcoef : 0;
    size_t a_f_elems = (N_val > 0) ? (size_t)lda_f_exp * N_val : 0;
    size_t b_f_elems = (N_val > 0) ? (size_t)ldb_f_exp * max_mp : 0; 
    size_t c_f_elems = (N_val > 0) ? (size_t)ldc_f_exp * N_val : 0; 
    size_t d_f_elems = (max_mp > 0) ? (size_t)ldd_f_exp * max_mp : 0;

    if (row_major_flag) {
        if (pcoeff_f_elems > 0 && pcoeff_c_in != NULL) {
            pcoeff_cm = (double*)malloc(pcoeff_f_elems * sizeof(double)); CHECK_ALLOC(pcoeff_cm);
            pcoeff_f_ptr = pcoeff_cm;
            for (int k = 0; k < kpcoef; ++k) {
                slicot_transpose_to_fortran_with_ld(
                    pcoeff_c_in + k * (size_t)ldpcoeff_c_rows * ldpcoeff_c_cols,
                    pcoeff_cm + k * (size_t)ldpco1_f_exp * ldpco2_f_exp, 
                    porm_f, porm_f, ldpcoeff_c_cols, ldpco1_f_exp, sizeof(double));
            }
        } else { pcoeff_f_ptr = NULL; }

        if (qcoeff_f_elems > 0 && qcoeff_c_in != NULL) {
            qcoeff_cm = (double*)malloc(qcoeff_f_elems * sizeof(double)); CHECK_ALLOC(qcoeff_cm);
            qcoeff_f_ptr = qcoeff_cm;
            memset(qcoeff_cm, 0, qcoeff_f_elems * sizeof(double));
            int q_c_rows_logical = (leri_f == 'L') ? p_c : m_c; 
            int q_c_cols_logical = (leri_f == 'L') ? m_c : p_c; 
            if (leri_f == 'R') { 
                 q_c_rows_logical = porm_f; 
                 q_c_cols_logical = porp_f; 
            }

            if (q_c_rows_logical > 0 && q_c_cols_logical > 0) {
                 for (int k = 0; k < kpcoef; ++k) {
                     slicot_transpose_to_fortran_with_ld(
                        qcoeff_c_in + k * (size_t)ldqcoeff_c_rows * ldqcoeff_c_cols,
                        qcoeff_cm + k * (size_t)ldqco1_f_exp * ldqco2_f_exp, 
                        q_c_rows_logical, q_c_cols_logical, ldqcoeff_c_cols, ldqco1_f_exp, sizeof(double));
                }
            }
        } else { qcoeff_f_ptr = NULL; }

        if (a_f_elems > 0 && a_out != NULL) { a_cm = (double*)malloc(a_f_elems * sizeof(double)); CHECK_ALLOC(a_cm); a_f_ptr = a_cm; } else {a_f_ptr = NULL;}
        if (b_f_elems > 0 && b_out != NULL) { b_cm = (double*)malloc(b_f_elems * sizeof(double)); CHECK_ALLOC(b_cm); b_f_ptr = b_cm; memset(b_cm, 0, b_f_elems * sizeof(double));} else {b_f_ptr = NULL;}
        if (c_f_elems > 0 && c_out != NULL) { c_cm = (double*)malloc(c_f_elems * sizeof(double)); CHECK_ALLOC(c_cm); c_f_ptr = c_cm; memset(c_cm, 0, c_f_elems * sizeof(double));} else {c_f_ptr = NULL;}
        if (d_f_elems > 0 && d_out != NULL) { d_cm = (double*)malloc(d_f_elems * sizeof(double)); CHECK_ALLOC(d_cm); d_f_ptr = d_cm; memset(d_cm, 0, d_f_elems * sizeof(double));} else {d_f_ptr = NULL;}

    } else { // Column-major C
        if (leri_f == 'R') { 
            if (pcoeff_f_elems > 0 && pcoeff_c_in != NULL) {
                pcoeff_cm = (double*)malloc(pcoeff_f_elems * sizeof(double)); CHECK_ALLOC(pcoeff_cm);
                memcpy(pcoeff_cm, pcoeff_c_in, pcoeff_f_elems * sizeof(double));
                pcoeff_f_ptr = pcoeff_cm;
            } else { pcoeff_f_ptr = NULL; }

            if (qcoeff_f_elems > 0 && qcoeff_c_in != NULL) {
                 qcoeff_cm = (double*)malloc(qcoeff_f_elems * sizeof(double)); CHECK_ALLOC(qcoeff_cm);
                 memcpy(qcoeff_cm, qcoeff_c_in, qcoeff_f_elems * sizeof(double));
                 qcoeff_f_ptr = qcoeff_cm;
            } else { qcoeff_f_ptr = NULL; }
        }
        if (N_val == 0) a_f_ptr = NULL;
        if (N_val == 0 || max_mp == 0) b_f_ptr = NULL; 
        if (max_mp == 0 || N_val == 0) c_f_ptr = NULL; 
        if (max_mp == 0) d_f_ptr = NULL; 
    }

    // Call Fortran routine
    F77_FUNC(tc04ad, TC04AD)(&leri_f, &m_c, &p_c, (int*)index_c, 
                             pcoeff_f_ptr, &ldpco1_f_exp, &ldpco2_f_exp,
                             qcoeff_f_ptr, &ldqco1_f_exp, &ldqco2_f_exp,
                             &N_val, rcond_out,
                             a_f_ptr, &lda_f_exp, b_f_ptr, &ldb_f_exp,
                             c_f_ptr, &ldc_f_exp, d_f_ptr, &ldd_f_exp,
                             iwork, dwork, &ldwork, &info, 1);

    if (info == 0) {
        *n_out = N_val; 
        if (row_major_flag) {
            if (N_val > 0 && a_cm != NULL && a_out != NULL) {
                slicot_transpose_to_c_with_ld(a_cm, a_out, N_val, N_val, lda_f_exp, lda_c, sizeof(double));
            }
            if (N_val > 0 && m_c > 0 && b_cm != NULL && b_out != NULL) { 
                slicot_transpose_to_c_with_ld(b_cm, b_out, N_val, m_c, ldb_f_exp, ldb_c, sizeof(double));
            }
            if (p_c > 0 && N_val > 0 && c_cm != NULL && c_out != NULL) { 
                slicot_transpose_to_c_with_ld(c_cm, c_out, p_c, N_val, ldc_f_exp, ldc_c, sizeof(double));
            }
            if (p_c > 0 && m_c > 0 && d_cm != NULL && d_out != NULL) { 
                slicot_transpose_to_c_with_ld(d_cm, d_out, p_c, m_c, ldd_f_exp, ldd_c, sizeof(double));
            }
        }
    }

cleanup:
    // free(iwork_dummy_query); // This was removed as query is skipped
    free(iwork);
    free(dwork);
    free(pcoeff_cm);
    free(qcoeff_cm);
    free(a_cm);
    free(b_cm);
    free(c_cm);
    free(d_cm);

    return info;
}
