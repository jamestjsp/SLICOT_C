/**
 * @file ab09md.c
 * @brief C wrapper implementation for SLICOT routine AB09MD
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB09MD,
 * which computes a reduced order model (Ar,Br,Cr) for the ALPHA-stable
 * part of an original state-space representation (A,B,C) using
 * Balance & Truncate methods.
 * This version skips workspace query and uses formula-based LDWORK.
 */

#include <stdlib.h> // For malloc, free
#include <ctype.h>  // For toupper
#include <stddef.h> // For size_t
#include <string.h> // For memset
#include <math.h>   // For MAX/MIN if needed (using slicot_utils.h MAX)

#include "ab09md.h"
#include "slicot_utils.h" 
#include "slicot_f77.h"   

/*
 * Declare the external Fortran routine using the F77_FUNC macro.
 */
extern void F77_FUNC(ab09md, AB09MD)(
    const char* dico,     
    const char* job,      
    const char* equil,    
    const char* ordsel,   
    const int* n,         
    const int* m,         
    const int* p,         
    int* nr,              
    const double* alpha,  
    double* a,            
    const int* lda,       
    double* b,            
    const int* ldb,       
    double* c,            
    const int* ldc,       
    int* ns,              
    double* hsv,          
    const double* tol,    
    int* iwork,           
    double* dwork,        
    const int* ldwork,    
    int* iwarn,           
    int* info,            
    int dico_len,         
    int job_len,          
    int equil_len,        
    int ordsel_len        
);


SLICOT_EXPORT
int slicot_ab09md(char dico_in, char job_in, char equil_in, char ordsel_in,
                  int n_in, int m_in, int p_in, int* nr_io, double alpha_in,
                  double* a_io, int lda_in,
                  double* b_io, int ldb_in,
                  double* c_io, int ldc_in,
                  int* ns_out, double* hsv_out, double tol_in, int* iwarn_out,
                  int row_major)
{
    int info = 0;
    int* iwork_allocated_buffer = NULL;
    double* dwork_allocated_buffer = NULL;
    int liwork_actual_size = 0;
    int ldwork_actual = 0;

    const int dico_len = 1, job_len = 1, equil_len = 1, ordsel_len = 1;

    char dico_upper = toupper(dico_in);
    char job_upper = toupper(job_in);
    char equil_upper = toupper(equil_in);
    char ordsel_upper = toupper(ordsel_in);

    double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL;
    double *a_ptr, *b_ptr, *c_ptr;
    double *hsv_ptr;
    int lda_f, ldb_f, ldc_f;

    /* --- Input Parameter Validation --- */
    if (dico_upper != 'C' && dico_upper != 'D') { info = -1; goto cleanup; }
    if (job_upper != 'B' && job_upper != 'N') { info = -2; goto cleanup; }
    if (equil_upper != 'S' && equil_upper != 'N') { info = -3; goto cleanup; }
    if (ordsel_upper != 'F' && ordsel_upper != 'A') { info = -4; goto cleanup; }
    if (n_in < 0) { info = -5; goto cleanup; }
    if (m_in < 0) { info = -6; goto cleanup; }
    if (p_in < 0) { info = -7; goto cleanup; }
    if (nr_io == NULL) { info = -8; goto cleanup; } // NR must be provided
    if (ordsel_upper == 'F' && (*nr_io < 0 || *nr_io > n_in) ) { info = -8; goto cleanup; } // NR validity
    if (dico_upper == 'C' && alpha_in > 0.0) { info = -9; goto cleanup; }
    if (dico_upper == 'D' && (alpha_in < 0.0 || alpha_in > 1.0)) { info = -9; goto cleanup; }

    if (n_in > 0 && a_io == NULL) { info = -10; goto cleanup; }
    // B can be NULL if M=0, even if N>0
    if (n_in > 0 && m_in > 0 && b_io == NULL) { info = -12; goto cleanup; }
    // C can be NULL if P=0, even if N>0
    if (p_in > 0 && n_in > 0 && c_io == NULL) { info = -14; goto cleanup; }
    
    if (ns_out == NULL) { info = -16; goto cleanup; }
    // HSV can be NULL if N=0
    if (n_in > 0 && hsv_out == NULL) { info = -17; goto cleanup; }
    if (iwarn_out == NULL) { info = -22; goto cleanup; } // Parameter 22 in Fortran is IWARN

    // Leading dimension checks
    int min_lda_f_val = MAX(1, n_in);
    int min_ldb_f_val = MAX(1, n_in); // LDB is N-by-M, so leading dim depends on N
    int min_ldc_f_val = MAX(1, p_in); // LDC is P-by-N, so leading dim depends on P

    if (row_major) {
        // For row-major C, input lda_in is number of columns of A (which is N for A, M for B, N for C in test)
        if (n_in > 0 && lda_in < n_in) { info = -11; goto cleanup; } // A is N x N, C lda is N
        if (n_in > 0 && m_in > 0 && ldb_in < m_in) { info = -13; goto cleanup; } // B is N x M, C ldb is M
        if (p_in > 0 && n_in > 0 && ldc_in < n_in) { info = -15; goto cleanup; } // C is P x N, C ldc is N
        else if (p_in > 0 && n_in == 0 && ldc_in < 1) {info = -15; goto cleanup;} // LDC must be >=1 if N=0, P>0
    } else { // Column-major C
        if (n_in > 0 && lda_in < min_lda_f_val) { info = -11; goto cleanup; }
        if (n_in > 0 && m_in > 0 && ldb_in < min_ldb_f_val) { info = -13; goto cleanup; } 
        if (p_in > 0 && n_in > 0 && ldc_in < min_ldc_f_val) { info = -15; goto cleanup; } 
        else if (p_in > 0 && n_in == 0 && ldc_in < min_ldc_f_val ) {info = -15; goto cleanup;}
    }

    /* --- Workspace Allocation --- */
    // LIWORK
    if (job_upper == 'N') {
        liwork_actual_size = MAX(1, n_in); 
        if (n_in > 0) {
            iwork_allocated_buffer = (int*)malloc((size_t)liwork_actual_size * sizeof(int));
            if (iwork_allocated_buffer == NULL && liwork_actual_size > 0) { info = SLICOT_MEMORY_ERROR; goto cleanup; }
            if (iwork_allocated_buffer) memset(iwork_allocated_buffer, 0, (size_t)liwork_actual_size * sizeof(int));
        } else {
            liwork_actual_size = 0; 
            iwork_allocated_buffer = NULL;
        }
    } else { // JOB == 'B'
        liwork_actual_size = 0; 
        iwork_allocated_buffer = NULL;
    }

    // LDWORK: Directly use the formula, no query.
    // LDWORK >= MAX(1,N*(2*N+MAX(N,M,P)+5) + N*(N+1)/2).
    if (n_in == 0) {
        ldwork_actual = 1;
    } else {
        long long term1 = (long long)2 * n_in;
        long long term_max = MAX(n_in, MAX(m_in, p_in));
        long long term_sum_paren = term1 + term_max + 5;
        long long term_mult = (long long)n_in * term_sum_paren;
        long long term_tri = ((long long)n_in * (n_in + 1)) / 2;
        ldwork_actual = (int)(term_mult + term_tri);
    }
    ldwork_actual = MAX(1, ldwork_actual); 
    
    dwork_allocated_buffer = (double*)malloc((size_t)ldwork_actual * sizeof(double));
    if (dwork_allocated_buffer == NULL && ldwork_actual > 0) { info = SLICOT_MEMORY_ERROR; goto cleanup; }
    if (dwork_allocated_buffer) memset(dwork_allocated_buffer, 0, (size_t)ldwork_actual * sizeof(double));


    /* --- Prepare Arrays for Computational Call --- */
    size_t elem_size = sizeof(double);
    if (row_major) {
        // Fortran leading dimensions for temporary column-major arrays
        lda_f = MAX(1, n_in); 
        ldb_f = MAX(1, n_in); // B is N x M, so Fortran LD is N
        ldc_f = MAX(1, p_in); // C is P x N, so Fortran LD is P

        size_t a_size_bytes = (size_t)n_in * n_in * elem_size; // A is N x N
        size_t b_size_bytes = (size_t)n_in * m_in * elem_size; // B is N x M
        size_t c_size_bytes = (size_t)p_in * n_in * elem_size; // C is P x N

        if (n_in > 0) { 
            a_cm = (double*)malloc(a_size_bytes);
            if (a_cm == NULL && a_size_bytes > 0) { info = SLICOT_MEMORY_ERROR; goto cleanup; }
            if (a_cm) slicot_transpose_to_fortran_with_ld(a_io, a_cm, n_in, n_in, lda_in, lda_f, elem_size);
        }
        if (n_in > 0 && m_in > 0) {
            b_cm = (double*)malloc(b_size_bytes);
            if (b_cm == NULL && b_size_bytes > 0) { info = SLICOT_MEMORY_ERROR; goto cleanup; }
            if (b_cm) slicot_transpose_to_fortran_with_ld(b_io, b_cm, n_in, m_in, ldb_in, ldb_f, elem_size);
        }
        if (p_in > 0 && n_in > 0) {
            c_cm = (double*)malloc(c_size_bytes);
            if (c_cm == NULL && c_size_bytes > 0) { info = SLICOT_MEMORY_ERROR; goto cleanup; }
            if (c_cm) slicot_transpose_to_fortran_with_ld(c_io, c_cm, p_in, n_in, ldc_in, ldc_f, elem_size);
        }
        
        a_ptr = (n_in > 0) ? a_cm : NULL; 
        b_ptr = (n_in > 0 && m_in > 0) ? b_cm : NULL; 
        c_ptr = (p_in > 0 && n_in > 0) ? c_cm : NULL;
    } else { // Column-major C
        lda_f = lda_in; 
        ldb_f = ldb_in; 
        ldc_f = ldc_in;
        a_ptr = (n_in > 0) ? a_io : NULL;
        b_ptr = (n_in > 0 && m_in > 0) ? b_io : NULL;
        c_ptr = (p_in > 0 && n_in > 0) ? c_io : NULL;
    }
    hsv_ptr = (n_in > 0) ? hsv_out : NULL;
    
    // Ensure Fortran LDs are at least 1, even if dimensions are 0
    lda_f = MAX(1, lda_f);
    ldb_f = MAX(1, ldb_f);
    ldc_f = MAX(1, ldc_f);


    F77_FUNC(ab09md, AB09MD)(&dico_upper, &job_upper, &equil_upper, &ordsel_upper,
                              &n_in, &m_in, &p_in, nr_io, &alpha_in,
                              a_ptr, &lda_f, b_ptr, &ldb_f, c_ptr, &ldc_f,
                              ns_out, hsv_ptr, &tol_in, 
                              iwork_allocated_buffer, dwork_allocated_buffer, &ldwork_actual,
                              iwarn_out, &info,
                              dico_len, job_len, equil_len, ordsel_len);

    if (row_major && info == 0) {
        int nr_val = *nr_io; // Reduced order
        if (nr_val >= 0) { 
            // A is overwritten with Ar (nr_val x nr_val)
            // Source a_ptr (a_cm) is column-major, nr_val rows, nr_val columns, with leading dimension lda_f.
            // Destination a_io is row-major, nr_val rows, nr_val columns, with leading dimension lda_in.
            if (a_cm != NULL && nr_val > 0) {
                slicot_transpose_to_c_with_ld(a_cm, a_io, nr_val, nr_val, 
                                             lda_f, lda_in, elem_size);
            }
            // B is overwritten with Br (nr_val x m_in)
            // Source b_ptr (b_cm) is column-major, nr_val rows, m_in columns, with leading dimension ldb_f.
            // Destination b_io is row-major, nr_val rows, m_in columns, with leading dimension ldb_in.
            if (b_cm != NULL && nr_val > 0 && m_in > 0) {
                slicot_transpose_to_c_with_ld(b_cm, b_io, nr_val, m_in, 
                                             ldb_f, ldb_in, elem_size);
            }
            // C is overwritten with Cr (p_in x nr_val)
            // Source c_ptr (c_cm) is column-major, p_in rows, nr_val columns, with leading dimension ldc_f.
            // Destination c_io is row-major, p_in rows, nr_val columns, with leading dimension ldc_in.
            if (c_cm != NULL && nr_val > 0 && p_in > 0) {
                slicot_transpose_to_c_with_ld(c_cm, c_io, p_in, nr_val, 
                                             ldc_f, ldc_in, elem_size);
            }
        }
    }

cleanup:
    free(dwork_allocated_buffer);
    free(iwork_allocated_buffer); 
    if (a_cm) free(a_cm); 
    if (b_cm) free(b_cm);
    if (c_cm) free(c_cm);

    return info;
}
