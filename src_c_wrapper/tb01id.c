/**
 * @file tb01id.c
 * @brief C wrapper for SLICOT routine TB01ID.
 * @details Balances a system matrix corresponding to a triplet (A,B,C).
 * Matrices A, B, C are modified. MAXRED is input/output. SCALE is output.
 * No explicit workspace arrays are required by this routine.
 * Input/output matrix format is handled via the row_major parameter.
 */

#include <stdlib.h> // For malloc, free
#include <string.h> // For memcpy
#include <ctype.h>  // For toupper
#include <math.h>   // For fmax (from C standard library)

#include "tb01id.h"       // Public header for this wrapper
#include "slicot_utils.h" // Provides CHECK_ALLOC, SLICOT_MEMORY_ERROR, MAX, MIN, transpose functions etc.
#include "slicot_f77.h"   // Provides F77_FUNC macro

/* External Fortran routine declaration */
extern void F77_FUNC(tb01id, TB01ID)(
    const char* job,
    const int* n, const int* m, const int* p,
    double* maxred, /* input/output */
    double* a, const int* lda, /* input/output */
    double* b, const int* ldb, /* input/output */
    double* c, const int* ldc, /* input/output */
    double* scale, /* output */
    int* info,     /* output */
    size_t job_len
);

SLICOT_EXPORT
int slicot_tb01id(
    char job_param,
    int n_param, int m_param, int p_param,
    double* maxred_io,      /* input/output */
    double* a_io, int lda,  /* input/output */
    double* b_io, int ldb,  /* input/output */
    double* c_io, int ldc,  /* input/output */
    double* scale_out,      /* output */
    int row_major)
{
    // 1. Variable declarations
    int info = 0;

    // Column-major copies for input/output matrices
    double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL;

    // Fortran-style leading dimensions
    int lda_f, ldb_f, ldc_f;
    
    char job_f = toupper(job_param);

    // 2. Input parameter validation
    if (strchr("ABCN", job_f) == NULL) { info = -1; goto cleanup; }
    if (n_param < 0) { info = -2; goto cleanup; }
    if (m_param < 0) { info = -3; goto cleanup; }
    if (p_param < 0) { info = -4; goto cleanup; }
    if (maxred_io == NULL) { info = -5; goto cleanup; } 

    // Validate pointers for matrices based on dimensions
    // A is N x N
    if (a_io == NULL && n_param > 0) { info = -6; goto cleanup; }
    // B is N x M. Not referenced if M=0.
    if (b_io == NULL && n_param > 0 && m_param > 0 && job_f != 'C' && job_f != 'N') { info = -8; goto cleanup; }
    // C is P x N. Not referenced if P=0.
    if (c_io == NULL && p_param > 0 && n_param > 0 && job_f != 'B' && job_f != 'N') { info = -10; goto cleanup; }
    // SCALE is N-dim output.
    if (scale_out == NULL && n_param > 0) { info = -12; goto cleanup; } 


    // Validate leading dimensions
    int min_lda_f = MAX(1, n_param);
    int min_ldb_f = MAX(1, n_param); 
    int min_ldc_f = MAX(1, p_param); 

    if (row_major) { // C LDA is number of columns
        if (n_param > 0 && lda < n_param) { info = -7; goto cleanup; }
        if (n_param > 0 && m_param > 0 && ldb < m_param) { info = -9; goto cleanup; }
        if (p_param > 0 && n_param > 0 && ldc < n_param) { info = -11; goto cleanup; }
    } else { // Column-major C (Fortran-style LDs)
        if (n_param > 0 && lda < min_lda_f) { info = -7; goto cleanup; }
        if (n_param > 0 && m_param > 0 && ldb < min_ldb_f) { info = -9; goto cleanup; }
        if (p_param > 0 && n_param > 0 && ldc < min_ldc_f) { info = -11; goto cleanup; }
    }
    if (info != 0) { goto cleanup; }

    // 3. Workspace Allocation - None explicit for TB01ID

    // 4. Memory allocation for column-major copies
    size_t a_nelems = (size_t)n_param * n_param; if(n_param == 0) a_nelems = 0;
    size_t b_nelems = (size_t)n_param * m_param; if(n_param == 0 || m_param == 0) b_nelems = 0;
    size_t c_nelems = (size_t)p_param * n_param; if(p_param == 0 || n_param == 0) c_nelems = 0;

    if (row_major) {
        if (a_nelems > 0) { a_cm = (double*)malloc(a_nelems * sizeof(double)); CHECK_ALLOC(a_cm); }
        if (b_nelems > 0 && (job_f == 'A' || job_f == 'B')) { b_cm = (double*)malloc(b_nelems * sizeof(double)); CHECK_ALLOC(b_cm); }
        if (c_nelems > 0 && (job_f == 'A' || job_f == 'C')) { c_cm = (double*)malloc(c_nelems * sizeof(double)); CHECK_ALLOC(c_cm); }
    }

    // 5. Prepare Fortran parameters and perform conversions
    double* a_ptr = a_io; double* b_ptr = b_io; double* c_ptr = c_io;
    
    lda_f=lda; ldb_f=ldb; ldc_f=ldc;

    if (row_major) {
        lda_f = min_lda_f; 
        ldb_f = min_ldb_f; 
        ldc_f = min_ldc_f;

        if (a_nelems > 0) { slicot_transpose_to_fortran_with_ld(a_io, a_cm, n_param, n_param, lda, lda_f, sizeof(double)); a_ptr = a_cm; } 
        else {a_ptr = (n_param == 0 ? NULL : a_io);}
        
        if (b_nelems > 0 && (job_f == 'A' || job_f == 'B')) { slicot_transpose_to_fortran_with_ld(b_io, b_cm, n_param, m_param, ldb, ldb_f, sizeof(double)); b_ptr = b_cm; } 
        else {b_ptr = (n_param == 0 || m_param == 0 ? NULL : b_io);}
        
        if (c_nelems > 0 && (job_f == 'A' || job_f == 'C')) { slicot_transpose_to_fortran_with_ld(c_io, c_cm, p_param, n_param, ldc, ldc_f, sizeof(double)); c_ptr = c_cm; } 
        else {c_ptr = (p_param == 0 || n_param == 0 ? NULL : c_io);}

    } else { // Column-major C
        if(a_nelems == 0 && n_param == 0) a_ptr = NULL; 
        if(b_nelems == 0 && (n_param == 0 || m_param == 0)) b_ptr = NULL;
        if(c_nelems == 0 && (p_param == 0 || n_param == 0)) c_ptr = NULL;
    }
    
    int n_f_call = n_param, m_f_call = m_param, p_f_call = p_param;

    // 7. Call Fortran function
    F77_FUNC(tb01id, TB01ID)(&job_f, &n_f_call, &m_f_call, &p_f_call,
                             maxred_io, 
                             a_ptr, &lda_f, b_ptr, &ldb_f, c_ptr, &ldc_f,
                             scale_out, &info,
                             1 /* job_len */); 

    // 8. Convert results back to row-major
    if (row_major && info == 0) {
        if (a_nelems > 0 && a_cm != NULL && a_io != NULL) {
            slicot_transpose_to_c_with_ld(a_cm, a_io, n_param, n_param, lda_f, lda, sizeof(double));
        }
        if (b_nelems > 0 && b_cm != NULL && b_io != NULL && (job_f == 'A' || job_f == 'B')) {
            slicot_transpose_to_c_with_ld(b_cm, b_io, n_param, m_param, ldb_f, ldb, sizeof(double));
        }
        if (c_nelems > 0 && c_cm != NULL && c_io != NULL && (job_f == 'A' || job_f == 'C')) {
            slicot_transpose_to_c_with_ld(c_cm, c_io, p_param, n_param, ldc_f, ldc, sizeof(double));
        }
    }

cleanup:
    if(row_major){
        free(a_cm); 
        if (job_f == 'A' || job_f == 'B') free(b_cm);
        if (job_f == 'A' || job_f == 'C') free(c_cm);
    }
    
    if (info == SLICOT_MEMORY_ERROR) {
       // Memory allocation error caught by CHECK_ALLOC
    }
    return info;
}
