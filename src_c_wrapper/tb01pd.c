/**
 * @file tb01pd.c
 * @brief C wrapper for SLICOT routine TB01PD.
 * @details Finds a reduced (controllable, observable, or minimal) state-space
 * representation (Ar,Br,Cr) for a given (A,B,C). Ar is upper block Hessenberg.
 * Matrices A, B, C are modified. NR is output.
 * Workspace (IWORK, DWORK) is allocated internally.
 * Input/output matrix format is handled via the row_major parameter.
 */

#include <stdlib.h> // For malloc, free
#include <string.h> // For memcpy, memset
#include <ctype.h>  // For toupper
#include <math.h>   // For fmax (from C standard library)

#include "tb01pd.h"       // Public header for this wrapper
#include "slicot_utils.h" // Provides CHECK_ALLOC, SLICOT_MEMORY_ERROR, MAX, MIN, transpose functions etc.
#include "slicot_f77.h"   // Provides F77_FUNC macro

/* External Fortran routine declaration */
extern void F77_FUNC(tb01pd, TB01PD)(
    const char* job, const char* equil,
    const int* n, const int* m, const int* p,
    double* a, const int* lda, /* input/output */
    double* b, const int* ldb, /* input/output */
    double* c, const int* ldc, /* input/output */
    int* nr,                   /* output */
    const double* tol,
    int* iwork,   /* workspace */
    double* dwork, const int* ldwork, /* workspace */
    int* info,    /* output */
    size_t job_len, size_t equil_len
);

SLICOT_EXPORT
int slicot_tb01pd(
    char job_param, char equil_param,
    int n_param, int m_param, int p_param,
    double* a_io, int lda,
    double* b_io, int ldb,
    double* c_io, int ldc,
    int* nr_out,
    double tol_param,
    int row_major)
{
    // 1. Variable declarations
    int info = 0;
    int *iwork = NULL;
    double *dwork = NULL;
    int liwork = 0;
    int ldwork = 0;

    double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL;

    int lda_f, ldb_f, ldc_f;
    
    char job_f = toupper(job_param);
    char equil_f = toupper(equil_param);
    double tol_f = tol_param; 

    // 2. Input parameter validation
    if (strchr("MCO", job_f) == NULL) { info = -1; goto cleanup; }
    if (strchr("SN", equil_f) == NULL) { info = -2; goto cleanup; }
    if (n_param < 0) { info = -3; goto cleanup; }
    if (m_param < 0) { info = -4; goto cleanup; }
    if (p_param < 0) { info = -5; goto cleanup; }

    int b_fortran_cols_workspace = (job_f == 'C') ? m_param : MAX(m_param, p_param);
    // Ensure b_fortran_cols_workspace is at least 0
    b_fortran_cols_workspace = MAX(0, b_fortran_cols_workspace);


    if (a_io == NULL && n_param > 0) { info = -6; goto cleanup; }
    if (b_io == NULL && n_param > 0 && b_fortran_cols_workspace > 0) { info = -8; goto cleanup; }
    if (c_io == NULL && n_param > 0 && p_param > 0) { info = -10; goto cleanup; }
    if (nr_out == NULL) { info = -12; goto cleanup; }


    int min_lda_f = MAX(1, n_param);
    int min_ldb_f = MAX(1, n_param); 
    int min_ldc_f = (n_param == 0) ? 1 : MAX(1, MAX(m_param, p_param)); 

    int lda_c_cols = n_param;
    int ldb_c_cols_min_req = b_fortran_cols_workspace; 
    int ldc_c_cols_min_req = n_param;


    if (row_major) { 
        if (n_param > 0 && lda < lda_c_cols) { info = -7; goto cleanup; }
        if (n_param > 0 && b_fortran_cols_workspace > 0 && ldb < ldb_c_cols_min_req) { info = -9; goto cleanup; }
        if (p_param > 0 && n_param > 0 && ldc < ldc_c_cols_min_req) { info = -11; goto cleanup; }
    } else { // Column-major C 
        if (n_param > 0 && lda < min_lda_f) { info = -7; goto cleanup; }
        if (n_param > 0 && b_fortran_cols_workspace > 0 && ldb < min_ldb_f) { info = -9; goto cleanup; } 
        if (n_param > 0 && ldc < min_ldc_f) { info = -11; goto cleanup; } 
        else if (n_param == 0 && p_param == 0 && ldc < 1 && c_io != NULL) {info = -11; goto cleanup;} 
        else if (n_param == 0 && p_param > 0 && ldc < MAX(1,m_param) && c_io != NULL) {info = -11; goto cleanup;}
    }
    if (info != 0) { goto cleanup; }

    // 3. Workspace Allocation
    liwork = n_param + MAX(1, MAX(m_param, p_param)); 
    if (n_param == 0 && m_param == 0 && p_param == 0) liwork = 1; 
    liwork = MAX(1, liwork);
    iwork = (int*)malloc((size_t)liwork * sizeof(int));
    CHECK_ALLOC(iwork);

    // LDWORK: Use formula directly. LDWORK >= MAX(1, N + MAX(N, 3*M, 3*P)).
    int term_max_n_3m_3p = MAX(n_param, MAX(3 * m_param, 3 * p_param));
    ldwork = n_param + MAX(1, term_max_n_3m_3p); 
    if (n_param == 0 && m_param == 0 && p_param == 0 && term_max_n_3m_3p == 0) { // Special case if N,M,P all 0
        ldwork = 1;
    } else if (n_param == 0 && term_max_n_3m_3p == 0) { // N=0, M=0, P=0
         ldwork = 1;
    }
    ldwork = MAX(1, ldwork);

    dwork = (double*)malloc((size_t)ldwork * sizeof(double));
    CHECK_ALLOC(dwork);

    // 4. Memory allocation for column-major copies
    size_t a_nelems = (size_t)n_param * n_param; if(n_param == 0) a_nelems = 0;
    
    size_t b_cm_nelems = (size_t)min_ldb_f * b_fortran_cols_workspace; 
    if(n_param == 0 || b_fortran_cols_workspace == 0) b_cm_nelems = 0;
    
    size_t c_cm_nelems = (size_t)min_ldc_f * n_param;
    if(n_param == 0) c_cm_nelems = 0; 


    if (row_major) {
        if (a_nelems > 0) { a_cm = (double*)malloc(a_nelems * sizeof(double)); CHECK_ALLOC(a_cm); }
        if (b_cm_nelems > 0) { b_cm = (double*)malloc(b_cm_nelems * sizeof(double)); CHECK_ALLOC(b_cm); }
        if (c_cm_nelems > 0) { c_cm = (double*)malloc(c_cm_nelems * sizeof(double)); CHECK_ALLOC(c_cm); }
    }

    // 5. Prepare Fortran parameters and perform conversions
    double* a_ptr = a_io; double* b_ptr = b_io; double* c_ptr = c_io;
    
    lda_f=lda; ldb_f=ldb; ldc_f=ldc;

    if (row_major) {
        lda_f = min_lda_f; 
        ldb_f = min_ldb_f; 
        ldc_f = min_ldc_f;

        if (a_nelems > 0) { slicot_transpose_to_fortran_with_ld(a_io, a_cm, n_param, n_param, lda, lda_f, sizeof(double)); a_ptr = a_cm; } 
        else {a_ptr = NULL;}
        
        if (b_cm_nelems > 0) { 
            // Initialize b_cm to zeros, then copy the M columns if M > 0
            memset(b_cm, 0, b_cm_nelems * sizeof(double));
            if (m_param > 0 && n_param > 0) { 
                slicot_transpose_to_fortran_with_ld(b_io, b_cm, n_param, m_param, ldb, ldb_f, sizeof(double)); 
            }
            b_ptr = b_cm; 
        } else {b_ptr = NULL;}
        
        if (c_cm_nelems > 0) { 
            memset(c_cm, 0, c_cm_nelems * sizeof(double));
            if (p_param > 0 && n_param > 0) { 
                 slicot_transpose_to_fortran_with_ld(c_io, c_cm, p_param, n_param, ldc, ldc_f, sizeof(double)); 
            }
            // If min_ldc_f (Fortran LDC) > p_param, rest of c_cm is workspace (already zeroed)
            c_ptr = c_cm; 
        } else {c_ptr = NULL;}

    } else { // Column-major C
        if(a_nelems == 0) a_ptr = NULL; 
        if(b_cm_nelems == 0) b_ptr = NULL; // Use b_cm_nelems as it reflects Fortran structure
        if(c_cm_nelems == 0) c_ptr = NULL;
    }
    
    int n_f_call = n_param, m_f_call = m_param, p_f_call = p_param;
    int nr_f_call; 
    int ldwork_f_call = ldwork;

    // 7. Call Fortran function
    F77_FUNC(tb01pd, TB01PD)(&job_f, &equil_f, &n_f_call, &m_f_call, &p_f_call,
                             a_ptr, &lda_f, b_ptr, &ldb_f, c_ptr, &ldc_f,
                             &nr_f_call, &tol_f, iwork, 
                             dwork, &ldwork_f_call, &info,
                             1,1); 

    if (nr_out != NULL) {
        *nr_out = nr_f_call;
    }
    int current_nr = nr_f_call; 

    // 8. Convert results back to row-major
    if (row_major && info == 0) {
        if (a_nelems > 0 && current_nr > 0 && a_cm != NULL && a_io != NULL) {
            slicot_transpose_to_c_with_ld(a_cm, a_io, current_nr, current_nr, lda_f, lda, sizeof(double));
        }
        // Output Br is NR x M. It's in the first M columns of b_cm.
        if (b_cm_nelems > 0 && current_nr > 0 && m_param > 0 && b_cm != NULL && b_io != NULL) {
            slicot_transpose_to_c_with_ld(b_cm, b_io, current_nr, m_param, ldb_f, ldb, sizeof(double));
        }
        // Output Cr is P x NR. It's in the first P rows of c_cm (if LDC_F > P).
        if (c_cm_nelems > 0 && p_param > 0 && current_nr > 0 && c_cm != NULL && c_io != NULL) {
            slicot_transpose_to_c_with_ld(c_cm, c_io, p_param, current_nr, ldc_f, ldc, sizeof(double));
        }
    }

cleanup:
    free(iwork);
    free(dwork);
    if(row_major){
        free(a_cm); free(b_cm); free(c_cm);
    }
    
    if (info == SLICOT_MEMORY_ERROR) {
       // Memory allocation error
    }
    return info;
}
