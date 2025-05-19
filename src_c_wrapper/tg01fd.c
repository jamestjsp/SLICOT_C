/**
 * @file tg01fd.c
 * @brief C wrapper for SLICOT routine TG01FD.
 */

#include <stdlib.h>
#include <string.h> // For memcpy, memset
#include <ctype.h>  // For toupper
#include <stdio.h>  // For debugging if needed

#include "tg01fd.h"
#include "slicot_utils.h" // Provides MAX, MIN, CHECK_ALLOC, transpose functions
#include "slicot_f77.h"   // Provides F77_FUNC

/* External Fortran routine declaration */
extern void F77_FUNC(tg01fd, TG01FD)(
    const char* compq, const char* compz, const char* joba,
    const int* l, const int* n, const int* m, const int* p,
    double* a, const int* lda,
    double* e, const int* lde,
    double* b, const int* ldb,
    double* c, const int* ldc,
    double* q, const int* ldq,
    double* z, const int* ldz,
    int* ranke, int* rnka22,
    const double* tol, int* iwork, double* dwork, const int* ldwork,
    int* info,
    size_t compq_len, size_t compz_len, size_t joba_len);

SLICOT_EXPORT
int slicot_tg01fd(char compq_in, char compz_in, char joba_in,
                  int l_in, int n_in, int m_in, int p_in,
                  double* a_c, int lda_c,
                  double* e_c, int lde_c,
                  double* b_c, int ldb_c,
                  double* c_c, int ldc_c,
                  double* q_c, int ldq_c,
                  double* z_c, int ldz_c,
                  int* ranke_out, int* rnka22_out,
                  double tol_in, int row_major)
{
    int info = 0;
    char compq_upper, compz_upper, joba_upper;

    int l_f = l_in; int n_f = n_in; int m_f = m_in; int p_f = p_in;

    double* a_ptr = a_c; double* e_ptr = e_c;
    double* b_ptr = b_c; double* c_ptr = c_c;
    double* q_ptr = q_c; double* z_ptr = z_c;

    int lda_f = lda_c; int lde_f = lde_c;
    int ldb_f = ldb_c; int ldc_f = ldc_c;
    int ldq_f = ldq_c; int ldz_f = ldz_c;

    double *a_cm = NULL; double *e_cm = NULL;
    double *b_cm = NULL; double *c_cm = NULL;
    double *q_cm = NULL; double *z_cm = NULL;

    int *iwork = NULL;
    double *dwork = NULL;
    int ldwork_calc;
    
    double dummy_double_array[1] = {0.0}; 
    int dummy_int_array[1] = {0};


    // --- 1. Validate input parameters ---
    compq_upper = toupper(compq_in); compz_upper = toupper(compz_in); joba_upper  = toupper(joba_in);

    if (compq_upper!='N' && compq_upper!='I' && compq_upper!='U') { info = -1; goto cleanup; }
    if (compz_upper!='N' && compz_upper!='I' && compz_upper!='U') { info = -2; goto cleanup; }
    if (joba_upper!='N' && joba_upper!='R' && joba_upper!='T') { info = -3; goto cleanup; }

    if (l_f < 0) { info = -4; goto cleanup; }
    if (n_f < 0) { info = -5; goto cleanup; }
    if (m_f < 0) { info = -6; goto cleanup; }
    if (p_f < 0) { info = -7; goto cleanup; }
    // TOL is checked by Fortran (tol < 1.0 by SLICOT docs, but not for >= 0 in TG01FD specifically)

    // A, E (l x n)
    if (l_f > 0 && n_f > 0) {
        if (a_c == NULL) { info = -8; goto cleanup; }
        if (row_major) { if (lda_c < MAX(1,n_f)) { info = -9; goto cleanup; } }
        else { if (lda_c < MAX(1,l_f)) { info = -9; goto cleanup; } }
        if (e_c == NULL) { info = -10; goto cleanup; }
        if (row_major) { if (lde_c < MAX(1,n_f)) { info = -11; goto cleanup; } }
        else { if (lde_c < MAX(1,l_f)) { info = -11; goto cleanup; } }
    } else { 
        if (a_c != NULL && lda_c < 1) { info = -9; goto cleanup; }
        if (e_c != NULL && lde_c < 1) { info = -11; goto cleanup; }
    }
    
    // B (l x m)
    if (m_f > 0) { 
        if (l_f > 0 && b_c == NULL) { info = -12; goto cleanup; }
        if (l_f > 0 && b_c != NULL) {
            if (row_major) { if (ldb_c < MAX(1,m_f)) { info = -13; goto cleanup; } }
            else { if (ldb_c < MAX(1,l_f)) { info = -13; goto cleanup; } }
        } else if (l_f == 0 && b_c != NULL && ldb_c < 1 ) {info = -13; goto cleanup;}
    } else { 
        if (b_c != NULL && ldb_c < 1) { info = -13; goto cleanup; } 
    }

    // C (p x n)
    if (p_f > 0) { 
        if (n_f > 0 && c_c == NULL) { info = -14; goto cleanup; }
        if (n_f > 0 && c_c != NULL) {
            if (row_major) { if (ldc_c < MAX(1,n_f)) { info = -15; goto cleanup; } }
            else { if (ldc_c < MAX(1,p_f)) { info = -15; goto cleanup; } }
        } else if (n_f == 0 && c_c != NULL && ldc_c < 1) {info = -15; goto cleanup;}
    } else { 
        if (c_c != NULL && ldc_c < 1) { info = -15; goto cleanup; }
    }
    
    // Q (l x l)
    if (compq_upper != 'N') {
        if (l_f > 0 && q_c == NULL) { info = -16; goto cleanup; }
        if (l_f > 0 && q_c != NULL) {
             if (row_major) { if (ldq_c < MAX(1,l_f)) { info = -17; goto cleanup; } }
             else { if (ldq_c < MAX(1,l_f)) { info = -17; goto cleanup; } }
        } else if (l_f == 0 && q_c != NULL && ldq_c < 1) { info = -17; goto cleanup;}
    } else { if (q_c != NULL && ldq_c < 1) { info = -17; goto cleanup;} }

    // Z (n x n)
    if (compz_upper != 'N') {
        if (n_f > 0 && z_c == NULL) { info = -18; goto cleanup; }
        if (n_f > 0 && z_c != NULL) {
            if (row_major) { if (ldz_c < MAX(1,n_f)) { info = -19; goto cleanup; } }
            else { if (ldz_c < MAX(1,n_f)) { info = -19; goto cleanup; } }
        } else if (n_f == 0 && z_c != NULL && ldz_c < 1 ) { info = -19; goto cleanup;}
    } else { if (z_c != NULL && ldz_c < 1) { info = -19; goto cleanup;} }

    if (ranke_out == NULL) { info = -20; goto cleanup; }
    if (joba_upper != 'N' && rnka22_out == NULL) { info = -21; goto cleanup; }
    
    if (info != 0) { goto cleanup; }

    // --- 2. Allocate Workspace ---
    if (n_f > 0) {
        iwork = (int*)malloc((size_t)n_f * sizeof(int)); CHECK_ALLOC(iwork);
    } else { iwork = dummy_int_array; } // Pass dummy if N=0
    
    ldwork_calc = MAX(1, n_f + p_f);
    ldwork_calc = MAX(ldwork_calc, MIN(l_f, n_f) + MAX(MAX(1, 3*n_f-1), MAX(m_f, l_f)));
    ldwork_calc = MAX(1, ldwork_calc); 

    dwork = (double*)malloc((size_t)ldwork_calc * sizeof(double));
    CHECK_ALLOC(dwork);

    // --- 3. Handle row_major conversions & Fortran LDs ---
    size_t a_e_size_elems = (size_t)l_f * n_f;
    size_t b_size_elems = (size_t)l_f * m_f;
    size_t c_size_elems = (size_t)p_f * n_f;
    size_t q_size_elems = (size_t)l_f * l_f;
    size_t z_size_elems = (size_t)n_f * n_f;

    // Set Fortran leading dimensions
    lda_f = MAX(1, l_f); lde_f = MAX(1, l_f);
    ldb_f = (m_f > 0) ? MAX(1, l_f) : 1;
    ldc_f = (p_f > 0) ? MAX(1, p_f) : 1;
    ldq_f = (compq_upper != 'N' && l_f > 0) ? MAX(1, l_f) : 1;
    ldz_f = (compz_upper != 'N' && n_f > 0) ? MAX(1, n_f) : 1;


    if (row_major) {
        if (l_f > 0 && n_f > 0) {
            if (a_c) { a_cm = (double*)malloc(a_e_size_elems * sizeof(double)); CHECK_ALLOC(a_cm);
                       slicot_transpose_to_fortran_with_ld(a_c, a_cm, l_f, n_f, lda_c, lda_f, sizeof(double)); a_ptr = a_cm; }
            if (e_c) { e_cm = (double*)malloc(a_e_size_elems * sizeof(double)); CHECK_ALLOC(e_cm);
                       slicot_transpose_to_fortran_with_ld(e_c, e_cm, l_f, n_f, lde_c, lde_f, sizeof(double)); e_ptr = e_cm; }
        }
        if (l_f > 0 && m_f > 0 && b_c) {
            b_cm = (double*)malloc(b_size_elems * sizeof(double)); CHECK_ALLOC(b_cm);
            slicot_transpose_to_fortran_with_ld(b_c, b_cm, l_f, m_f, ldb_c, ldb_f, sizeof(double)); b_ptr = b_cm;
        }
        if (p_f > 0 && n_f > 0 && c_c) {
            c_cm = (double*)malloc(c_size_elems * sizeof(double)); CHECK_ALLOC(c_cm);
            slicot_transpose_to_fortran_with_ld(c_c, c_cm, p_f, n_f, ldc_c, ldc_f, sizeof(double)); c_ptr = c_cm;
        }
        if (compq_upper != 'N' && l_f > 0 && q_c) {
            q_cm = (double*)malloc(q_size_elems * sizeof(double)); CHECK_ALLOC(q_cm);
            if (compq_upper == 'U') { 
                slicot_transpose_to_fortran_with_ld(q_c, q_cm, l_f, l_f, ldq_c, ldq_f, sizeof(double));
            } 
            q_ptr = q_cm;
        }
        if (compz_upper != 'N' && n_f > 0 && z_c) {
            z_cm = (double*)malloc(z_size_elems * sizeof(double)); CHECK_ALLOC(z_cm);
            if (compz_upper == 'U') { 
                slicot_transpose_to_fortran_with_ld(z_c, z_cm, n_f, n_f, ldz_c, ldz_f, sizeof(double));
            }
            z_ptr = z_cm;
        }
    } else { 
        // For column-major C, pointers and LDs are already set, just ensure safety for zero-dim
        if (!(l_f > 0 && n_f > 0)) { a_ptr = NULL; e_ptr = NULL; }
        if (!(l_f > 0 && m_f > 0)) { b_ptr = NULL; }
        if (!(p_f > 0 && n_f > 0)) { c_ptr = NULL; }
        if (!(compq_upper != 'N' && l_f > 0)) { q_ptr = NULL; }
        if (!(compz_upper != 'N' && n_f > 0)) { z_ptr = NULL; }
    }
    
    // Use dummy arrays for Fortran if logically zero-sized or not referenced and ptr is NULL
    double* final_a_ptr = (l_f > 0 && n_f > 0 && a_ptr != NULL) ? a_ptr : dummy_double_array;
    double* final_e_ptr = (l_f > 0 && n_f > 0 && e_ptr != NULL) ? e_ptr : dummy_double_array;
    double* final_b_ptr = (m_f > 0 && l_f > 0 && b_ptr != NULL) ? b_ptr : dummy_double_array; // B not ref if m=0
    double* final_c_ptr = (p_f > 0 && n_f > 0 && c_ptr != NULL) ? c_ptr : dummy_double_array; // C not ref if p=0
    double* final_q_ptr = (compq_upper != 'N' && l_f > 0 && q_ptr != NULL) ? q_ptr : dummy_double_array; // Q not ref if compq='N'
    double* final_z_ptr = (compz_upper != 'N' && n_f > 0 && z_ptr != NULL) ? z_ptr : dummy_double_array; // Z not ref if compz='N'
    int* final_iwork_ptr = (n_f > 0 && iwork != NULL) ? iwork : dummy_int_array;


    // --- 4. Call Fortran routine ---
    F77_FUNC(tg01fd, TG01FD)(&compq_upper, &compz_upper, &joba_upper,
                             &l_f, &n_f, &m_f, &p_f,
                             final_a_ptr, &lda_f, final_e_ptr, &lde_f,
                             final_b_ptr, &ldb_f, final_c_ptr, &ldc_f,
                             final_q_ptr, &ldq_f, final_z_ptr, &ldz_f,
                             ranke_out, rnka22_out, &tol_in,
                             final_iwork_ptr, dwork, &ldwork_calc, &info,
                             (size_t)1, (size_t)1, (size_t)1);

    // --- 5. Transpose results back if row_major and successful ---
    if (row_major && info == 0) {
        if (l_f > 0 && n_f > 0) {
            if (a_cm && a_c) slicot_transpose_to_c_with_ld(a_cm, a_c, l_f, n_f, lda_f, lda_c, sizeof(double));
            if (e_cm && e_c) slicot_transpose_to_c_with_ld(e_cm, e_c, l_f, n_f, lde_f, lde_c, sizeof(double));
        }
        if (l_f > 0 && m_f > 0 && b_cm && b_c) {
            slicot_transpose_to_c_with_ld(b_cm, b_c, l_f, m_f, ldb_f, ldb_c, sizeof(double));
        }
        if (p_f > 0 && n_f > 0 && c_cm && c_c) {
            slicot_transpose_to_c_with_ld(c_cm, c_c, p_f, n_f, ldc_f, ldc_c, sizeof(double));
        }
        if (compq_upper != 'N' && l_f > 0 && q_cm && q_c) {
            slicot_transpose_to_c_with_ld(q_cm, q_c, l_f, l_f, ldq_f, ldq_c, sizeof(double));
        }
        if (compz_upper != 'N' && n_f > 0 && z_cm && z_c) {
            slicot_transpose_to_c_with_ld(z_cm, z_c, n_f, n_f, ldz_f, ldz_c, sizeof(double));
        }
    }

cleanup:
    if (iwork != dummy_int_array) free(iwork); // Only free if not dummy
    free(dwork);
    if (row_major) {
        free(a_cm); free(e_cm);
        free(b_cm); free(c_cm);
        free(q_cm); free(z_cm);
    }
    return info;
}
