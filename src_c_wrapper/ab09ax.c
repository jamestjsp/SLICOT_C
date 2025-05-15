/**
 * @file ab09ax.c
 * @brief C wrapper implementation for SLICOT routine AB09AX
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB09AX,
 * which computes a reduced order model (Ar,Br,Cr) for a stable system
 * (A,B,C) using Balance & Truncate methods, where A is already in
 * real Schur form. It also returns the truncation matrices T and TI.
 * Workspace (IWORK, DWORK) is allocated internally by this wrapper.
 */

 #include <stdlib.h> // For malloc, free
 #include <ctype.h>  // For toupper
 #include <stddef.h> // For size_t
 #include <math.h>   // For MAX/MIN if needed
 
 // Include the header file for this wrapper
 #include "ab09ax.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Provides CHECK_ALLOC, SLICOT_MEMORY_ERROR, MAX, MIN, transpose functions
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions
 
 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  */
 extern void F77_FUNC(ab09ax, AB09AX)(
     const char* dico,       // CHARACTER*1 DICO
     const char* job,        // CHARACTER*1 JOB
     const char* ordsel,     // CHARACTER*1 ORDSEL
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     const int* p,           // INTEGER P
     int* nr,                // INTEGER NR (in/out)
     double* a,              // DOUBLE PRECISION A(LDA,*) (in/out)
     const int* lda,         // INTEGER LDA
     double* b,              // DOUBLE PRECISION B(LDB,*) (in/out)
     const int* ldb,         // INTEGER LDB
     double* c,              // DOUBLE PRECISION C(LDC,*) (in/out)
     const int* ldc,         // INTEGER LDC
     double* hsv,            // DOUBLE PRECISION HSV(*) (output)
     double* t,              // DOUBLE PRECISION T(LDT,*) (output)
     const int* ldt,         // INTEGER LDT
     double* ti,             // DOUBLE PRECISION TI(LDTI,*) (output)
     const int* ldti,        // INTEGER LDTI
     const double* tol,      // DOUBLE PRECISION TOL
     int* iwork,             // INTEGER IWORK(*)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* iwarn,             // INTEGER IWARN (output)
     int* info,              // INTEGER INFO (output)
     int dico_len,           // Hidden length
     int job_len,            // Hidden length
     int ordsel_len          // Hidden length
 );
 
 
 /* C wrapper function definition */
 SLICOT_EXPORT
 int slicot_ab09ax(char dico_in, char job_in, char ordsel_in,
                   int n_in, int m_in, int p_in, int* nr_io,
                   double* a_io, int lda_in,
                   double* b_io, int ldb_in,
                   double* c_io, int ldc_in,
                   double* hsv_out,
                   double* t_out, int ldt_in,
                   double* ti_out, int ldti_in,
                   double tol_in, int* iwarn_out,
                   int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork_actual = 0;
     double* dwork_allocated_buffer = NULL;
     int* iwork_allocated_buffer = NULL;
     int liwork_actual_size = 0;
 
     const int dico_len = 1, job_len = 1, ordsel_len = 1;
 
     char dico_upper = toupper(dico_in);
     char job_upper = toupper(job_in);
     char ordsel_upper = toupper(ordsel_in);
 
     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL;
     double *t_cm = NULL, *ti_cm = NULL; 
 
     /* Pointers to pass to Fortran */
     double *a_ptr, *b_ptr, *c_ptr;
     double *hsv_ptr;
     double *t_ptr, *ti_ptr;
     int lda_f, ldb_f, ldc_f;
     int ldt_f, ldti_f;
 
     /* --- Input Parameter Validation --- */
     /* Parameter indices for error reporting match Fortran routine's 1-based indexing */
     if (dico_upper != 'C' && dico_upper != 'D') { info = -1; goto cleanup; }
     if (job_upper != 'B' && job_upper != 'N') { info = -2; goto cleanup; }
     if (ordsel_upper != 'F' && ordsel_upper != 'A') { info = -3; goto cleanup; }
     if (n_in < 0) { info = -4; goto cleanup; }
     if (m_in < 0) { info = -5; goto cleanup; }
     if (p_in < 0) { info = -6; goto cleanup; }
     if (nr_io == NULL) { info = -7; goto cleanup; } // NR pointer itself
     if (ordsel_upper == 'F' && (*nr_io < 0 || *nr_io > n_in) ) { info = -7; goto cleanup; }
 
     // Check for NULL pointers for matrices if dimensions are non-zero
     // A, B, C, T, TI can be NULL if N=0. HSV also.
     if (n_in > 0 && a_io == NULL) { info = -8; goto cleanup; }
     if (n_in > 0 && m_in > 0 && b_io == NULL) { info = -10; goto cleanup; }
     if (p_in > 0 && n_in > 0 && c_io == NULL) { info = -12; goto cleanup; }
     if (n_in > 0 && hsv_out == NULL) { info = -14; goto cleanup; }
     // T and TI are outputs; pointers must be valid if NR (on exit) > 0.
     // This check is tricky before NR is known. Assume non-NULL if N > 0 for now.
     if (n_in > 0 && t_out == NULL) { info = -15; goto cleanup; }
     if (n_in > 0 && ti_out == NULL) { info = -17; goto cleanup; }
     if (iwarn_out == NULL) { info = -23; goto cleanup; } // IWARN pointer
 
     // Validate leading dimensions
     int min_lda_f_val = MAX(1, n_in);
     int min_ldb_f_val = MAX(1, n_in);
     int min_ldc_f_val = MAX(1, p_in);
     int min_ldt_f_val = MAX(1, n_in); 
     // LDTI depends on NR (output). For input validation, if N>0, LDTI must be at least 1.
     // A more precise check for LDTI against NR (output) will be done after the Fortran call.
     int min_ldti_f_val = (ordsel_upper == 'F' && *nr_io > 0) ? MAX(1, *nr_io) : 1;
 
 
     if (row_major) {
         if (n_in > 0 && lda_in < n_in) { info = -9; goto cleanup; }
         if (n_in > 0 && m_in > 0 && ldb_in < m_in) { info = -11; goto cleanup; }
         if (p_in > 0 && n_in > 0 && ldc_in < n_in) { info = -13; goto cleanup; }
         if (n_in > 0 && ldt_in < n_in) { info = -16; goto cleanup; } // T is N x NR, C LDA is N cols
         // TI is NR x N. For row-major, C LDTI is N cols.
         if (n_in > 0 && ldti_in < n_in) { info = -18; goto cleanup; }
     } else { // Column-major
         if (lda_in < min_lda_f_val) { info = -9; goto cleanup; }
         if (ldb_in < min_ldb_f_val) { info = -11; goto cleanup; }
         if (ldc_in < min_ldc_f_val) { info = -13; goto cleanup; }
         if (ldt_in < min_ldt_f_val) { info = -16; goto cleanup; }
         if (ldti_in < min_ldti_f_val && n_in > 0 && !(ordsel_upper == 'F' && *nr_io == 0) ) { info = -18; goto cleanup; }
     }
 
 
     /* --- Workspace Allocation --- */
     // LIWORK
     if (job_upper == 'N') {
         liwork_actual_size = MAX(1, n_in); // LIWORK = N if JOB = 'N'
         if (n_in > 0) {
             iwork_allocated_buffer = (int*)malloc((size_t)liwork_actual_size * sizeof(int));
             CHECK_ALLOC(iwork_allocated_buffer);
         } else {
             liwork_actual_size = 0;
             iwork_allocated_buffer = NULL;
         }
     } else { // JOB == 'B'
         liwork_actual_size = 0; // LIWORK = 0 if JOB = 'B'
         iwork_allocated_buffer = NULL;
     }

     // LDWORK: Use the formula provided in the documentation
     if (n_in == 0) {
         ldwork_actual = 1;
     } else {
         ldwork_actual = MAX(1, n_in * (MAX(n_in, MAX(m_in, p_in)) + 5) + (n_in * (n_in + 1)) / 2);
     }

     dwork_allocated_buffer = (double*)malloc((size_t)ldwork_actual * sizeof(double));
     CHECK_ALLOC(dwork_allocated_buffer);
 
 
     /* --- Prepare Arrays for Computational Call --- */
     size_t elem_size = sizeof(double);
     if (row_major) {
         /* --- Row-Major Case: Allocate CM buffers and transpose inputs --- */
         size_t a_size = (size_t)n_in * n_in;
         size_t b_size = (size_t)n_in * m_in;
         size_t c_size = (size_t)p_in * n_in;
         // T is N x NR_out, TI is NR_out x N. For allocation, use N for NR_out initially.
         size_t t_alloc_cols = n_in; // Max possible NR is N
         size_t ti_alloc_rows = n_in; 
         size_t t_size = (size_t)n_in * t_alloc_cols; 
         size_t ti_size = (size_t)ti_alloc_rows * n_in;
 
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
         if (t_size > 0) { t_cm = (double*)malloc(t_size * elem_size); CHECK_ALLOC(t_cm); } // T is output
         if (ti_size > 0) { ti_cm = (double*)malloc(ti_size * elem_size); CHECK_ALLOC(ti_cm); } // TI is output
 
         /* Transpose C (row-major) inputs to Fortran (column-major) copies */
         if (a_size > 0) slicot_transpose_to_fortran_with_ld(a_io, a_cm, n_in, n_in, lda_in, MAX(1,n_in), elem_size);
         if (b_size > 0) slicot_transpose_to_fortran_with_ld(b_io, b_cm, n_in, m_in, ldb_in, MAX(1,n_in), elem_size);
         if (c_size > 0) slicot_transpose_to_fortran_with_ld(c_io, c_cm, p_in, n_in, ldc_in, MAX(1,p_in), elem_size);
 
         lda_f = MAX(1, n_in); ldb_f = MAX(1, n_in); ldc_f = MAX(1, p_in);
         ldt_f = MAX(1, n_in); ldti_f = MAX(1, n_in); // Initial LDTI_F based on N, will be checked against NR_out
 
         a_ptr = (n_in > 0) ? a_cm : NULL; 
         b_ptr = (n_in > 0 && m_in > 0) ? b_cm : NULL; 
         c_ptr = (p_in > 0 && n_in > 0) ? c_cm : NULL;
         t_ptr = (n_in > 0) ? t_cm : NULL;   // Fortran writes to t_cm
         ti_ptr = (n_in > 0) ? ti_cm : NULL; // Fortran writes to ti_cm
     } else {
         /* --- Column-Major Case --- */
         lda_f = lda_in; ldb_f = ldb_in; ldc_f = ldc_in;
         ldt_f = ldt_in; ldti_f = ldti_in;
         a_ptr = (n_in > 0) ? a_io : NULL;
         b_ptr = (n_in > 0 && m_in > 0) ? b_io : NULL;
         c_ptr = (p_in > 0 && n_in > 0) ? c_io : NULL;
         t_ptr = (n_in > 0) ? t_out : NULL;
         ti_ptr = (n_in > 0) ? ti_out : NULL;
     }
     hsv_ptr = (n_in > 0) ? hsv_out : NULL;
 
     /* --- Call the computational routine --- */
     F77_FUNC(ab09ax, AB09AX)(&dico_upper, &job_upper, &ordsel_upper,
                               &n_in, &m_in, &p_in, nr_io, 
                               a_ptr, &lda_f, b_ptr, &ldb_f, c_ptr, &ldc_f,
                               hsv_ptr, t_ptr, &ldt_f, ti_ptr, &ldti_f,
                               &tol_in, iwork_allocated_buffer, dwork_allocated_buffer, &ldwork_actual,
                               iwarn_out, &info,
                               dico_len, job_len, ordsel_len);
 
     /* --- Copy back results from column-major temps to original row-major arrays if row_major --- */
     if (row_major && info == 0) {
         int nr_val = *nr_io; // Get the final reduced order
 
         // Validate LDTI for row-major C array: ldti_in (cols) must be >= N
         // And Fortran LDTI_F (rows) must be >= NR_val
         if (nr_val > 0) {
             if (ldti_f < nr_val) { info = -18; goto cleanup; } // Fortran LDTI too small
             if (ldti_in < n_in)  { info = -18; goto cleanup; } // C LDTI too small
         }
 
 
         if (nr_val >= 0) { 
             if (n_in > 0) { // Only transpose if original dimension was > 0
                 // A is overwritten with Ar (nr_val x nr_val)
                 if (a_cm != NULL) slicot_transpose_to_c_with_ld(a_cm, a_io, nr_val, nr_val, lda_f, lda_in, elem_size);
                 // B is overwritten with Br (nr_val x m_in)
                 if (m_in > 0 && b_cm != NULL) slicot_transpose_to_c_with_ld(b_cm, b_io, nr_val, m_in, ldb_f, ldb_in, elem_size);
                 // C is overwritten with Cr (p_in x nr_val)
                 if (p_in > 0 && c_cm != NULL) slicot_transpose_to_c_with_ld(c_cm, c_io, p_in, nr_val, ldc_f, ldc_in, elem_size);
             
                 // T is N x NR_val
                 if (t_cm != NULL && t_out != NULL) slicot_transpose_to_c_with_ld(t_cm, t_out, n_in, nr_val, ldt_f, ldt_in, elem_size);
                 // TI is NR_val x N
                 if (ti_cm != NULL && ti_out != NULL) slicot_transpose_to_c_with_ld(ti_cm, ti_out, nr_val, n_in, ldti_f, ldti_in, elem_size);
             }
         }
     }
     // HSV is 1D, filled directly.
 
 cleanup:
     /* --- Cleanup --- */
     free(dwork_allocated_buffer);
     free(iwork_allocated_buffer); 
     free(a_cm);
     free(b_cm);
     free(c_cm);
     free(t_cm);
     free(ti_cm);
 
     return info;
 }
