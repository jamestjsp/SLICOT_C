/**
 * @file ab09ad.c
 * @brief C wrapper implementation for SLICOT routine AB09AD
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB09AD,
 * which computes a reduced order model (Ar,Br,Cr) for a stable original
 * state-space representation (A,B,C) using Balance & Truncate methods.
 * Workspace (IWORK, DWORK) is allocated internally by this wrapper.
 * Input/output matrix format is handled via the row_major parameter.
 */

 #include <stdlib.h> // For malloc, free
 #include <ctype.h>  // For toupper
 #include <stddef.h> // For size_t
 #include <math.h>   // For fabs if needed for TOL comparison, though MAX/MIN are from slicot_utils
 
 // Include the header file for this wrapper
 #include "ab09ad.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Provides CHECK_ALLOC, SLICOT_MEMORY_ERROR, MAX, MIN, transpose functions
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions
 
 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * Note A, B, C are input/output. NR is input/output.
  */
 extern void F77_FUNC(ab09ad, AB09AD)(
     const char* dico,       // CHARACTER*1 DICO
     const char* job,        // CHARACTER*1 JOB
     const char* equil,      // CHARACTER*1 EQUIL
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
     const double* tol,      // DOUBLE PRECISION TOL
     int* iwork,             // INTEGER IWORK(*)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* iwarn,             // INTEGER IWARN (output)
     int* info,              // INTEGER INFO (output)
     int dico_len,           // Hidden length
     int job_len,            // Hidden length
     int equil_len,          // Hidden length
     int ordsel_len          // Hidden length
 );
 
 
 /* C wrapper function definition */
 SLICOT_EXPORT
 int slicot_ab09ad(char dico_in, char job_in, char equil_in, char ordsel_in,
                   int n_in, int m_in, int p_in, int* nr_io,
                   double* a_io, int lda_in,
                   double* b_io, int ldb_in,
                   double* c_io, int ldc_in,
                   double* hsv_out, double tol_in, int* iwarn_out,
                   int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork_actual = 0;
     double* dwork_allocated_buffer = NULL;
     int* iwork_allocated_buffer = NULL;
     int liwork_actual_size = 0;
 
     const int dico_len = 1, job_len = 1, equil_len = 1, ordsel_len = 1;
 
     char dico_upper = toupper(dico_in);
     char job_upper = toupper(job_in);
     char equil_upper = toupper(equil_in);
     char ordsel_upper = toupper(ordsel_in);
 
     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL;
     double *a_ptr, *b_ptr, *c_ptr;
     int lda_f, ldb_f, ldc_f; // Fortran leading dimensions
 
     /* --- Input Parameter Validation --- */
     /* Parameter indices for error reporting match Fortran routine's 1-based indexing */
     if (dico_upper != 'C' && dico_upper != 'D') { info = -1; goto cleanup; }
     if (job_upper != 'B' && job_upper != 'N') { info = -2; goto cleanup; }
     if (equil_upper != 'S' && equil_upper != 'N') { info = -3; goto cleanup; }
     if (ordsel_upper != 'F' && ordsel_upper != 'A') { info = -4; goto cleanup; }
     if (n_in < 0) { info = -5; goto cleanup; }
     if (m_in < 0) { info = -6; goto cleanup; }
     if (p_in < 0) { info = -7; goto cleanup; }
     if (nr_io == NULL) { info = -8; goto cleanup; } // NR pointer itself
     if (ordsel_upper == 'F' && (*nr_io < 0 || *nr_io > n_in) ) { info = -8; goto cleanup; }
 
     // Check for NULL pointers for matrices if dimensions are non-zero
     if (n_in > 0 && a_io == NULL) { info = -9; goto cleanup; }
     if (n_in > 0 && m_in > 0 && b_io == NULL) { info = -11; goto cleanup; }
     if (p_in > 0 && n_in > 0 && c_io == NULL) { info = -13; goto cleanup; }
     if (n_in > 0 && hsv_out == NULL) { info = -15; goto cleanup; } // HSV is always output if N > 0
     if (iwarn_out == NULL) { info = -20; goto cleanup; } // IWARN pointer
 
     // Validate leading dimensions
     // Fortran LDA is number of rows
     int min_lda_f_val = MAX(1, n_in);
     int min_ldb_f_val = MAX(1, n_in);
     int min_ldc_f_val = MAX(1, p_in);
 
     if (row_major) {
         // C LDA is number of columns for row-major
         if (n_in > 0 && lda_in < n_in) { info = -10; goto cleanup; } // A is N x N, so LDA (cols) >= N
         if (m_in > 0 && ldb_in < m_in) { info = -12; goto cleanup; } // B is N x M, so LDB (cols) >= M
         if (n_in > 0 && ldc_in < n_in) { info = -14; goto cleanup; } // C is P x N, so LDC (cols) >= N
     } else {
         // C LDA is number of rows for column-major (Fortran style)
         if (lda_in < min_lda_f_val) { info = -10; goto cleanup; }
         if (ldb_in < min_ldb_f_val) { info = -12; goto cleanup; }
         if (ldc_in < min_ldc_f_val) { info = -14; goto cleanup; }
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
         ldwork_actual = MAX(1, n_in * (2 * n_in + MAX(n_in, MAX(m_in, p_in)) + 5) + (n_in * (n_in + 1)) / 2);
     }

     dwork_allocated_buffer = (double*)malloc((size_t)ldwork_actual * sizeof(double));
     CHECK_ALLOC(dwork_allocated_buffer);
 
     /* --- Prepare arrays for column-major format if using row-major --- */
     size_t elem_size = sizeof(double);
     if (row_major) {
         size_t a_rows_orig = n_in; size_t a_cols_orig = n_in; size_t a_size_orig = a_rows_orig * a_cols_orig;
         size_t b_rows_orig = n_in; size_t b_cols_orig = m_in; size_t b_size_orig = b_rows_orig * b_cols_orig;
         size_t c_rows_orig = p_in; size_t c_cols_orig = n_in; size_t c_size_orig = c_rows_orig * c_cols_orig;
 
         lda_f = MAX(1, n_in); // Fortran LDA is number of rows of A
         ldb_f = MAX(1, n_in); // Fortran LDB is number of rows of B
         ldc_f = MAX(1, p_in); // Fortran LDC is number of rows of C
 
         if (a_size_orig > 0) {
             a_cm = (double*)malloc(a_size_orig * elem_size); CHECK_ALLOC(a_cm);
             slicot_transpose_to_fortran_with_ld(a_io, a_cm, a_rows_orig, a_cols_orig, lda_in, lda_f, elem_size);
             a_ptr = a_cm;
         } else {
             a_ptr = NULL; // Pass NULL if N=0
         }
 
         if (b_size_orig > 0) {
             b_cm = (double*)malloc(b_size_orig * elem_size); CHECK_ALLOC(b_cm);
             slicot_transpose_to_fortran_with_ld(b_io, b_cm, b_rows_orig, b_cols_orig, ldb_in, ldb_f, elem_size);
             b_ptr = b_cm;
         } else {
             b_ptr = NULL; // Pass NULL if N=0 or M=0
         }
 
         if (c_size_orig > 0) {
             c_cm = (double*)malloc(c_size_orig * elem_size); CHECK_ALLOC(c_cm);
             slicot_transpose_to_fortran_with_ld(c_io, c_cm, c_rows_orig, c_cols_orig, ldc_in, ldc_f, elem_size);
             c_ptr = c_cm;
         } else {
             c_ptr = NULL; // Pass NULL if P=0 or N=0
         }
     } else {
         /* Column-major case - use original arrays and LDs */
         lda_f = lda_in;
         ldb_f = ldb_in;
         ldc_f = ldc_in;
         a_ptr = (n_in > 0) ? a_io : NULL;
         b_ptr = (n_in > 0 && m_in > 0) ? b_io : NULL;
         c_ptr = (p_in > 0 && n_in > 0) ? c_io : NULL;
     }
     
     double* hsv_ptr = (n_in > 0) ? hsv_out : NULL;
 
 
     /* --- Call the computational routine --- */
     F77_FUNC(ab09ad, AB09AD)(&dico_upper, &job_upper, &equil_upper, &ordsel_upper,
                               &n_in, &m_in, &p_in, nr_io,
                               a_ptr, &lda_f,
                               b_ptr, &ldb_f,
                               c_ptr, &ldc_f,
                               hsv_ptr, &tol_in,
                               iwork_allocated_buffer, dwork_allocated_buffer, &ldwork_actual, /* Pass actual DWORK size */
                               iwarn_out, &info,
                               dico_len, job_len, equil_len, ordsel_len);
 
     /* --- Copy results back to row-major format if needed --- */
     if (row_major && info == 0) {
         int nr_val = *nr_io; // Get the final reduced order
 
         // A is overwritten with Ar (NR x NR)
         if (nr_val > 0 && a_cm != NULL) {
             slicot_transpose_to_c_with_ld(a_cm, a_io, nr_val, nr_val, lda_f, lda_in, elem_size);
         }
         // B is overwritten with Br (NR x M)
         if (nr_val > 0 && m_in > 0 && b_cm != NULL) {
             slicot_transpose_to_c_with_ld(b_cm, b_io, nr_val, m_in, ldb_f, ldb_in, elem_size);
         }
         // C is overwritten with Cr (P x NR)
         if (p_in > 0 && nr_val > 0 && c_cm != NULL) {
             slicot_transpose_to_c_with_ld(c_cm, c_io, p_in, nr_val, ldc_f, ldc_in, elem_size);
         }
         // HSV is 1D, filled directly by Fortran, no transpose needed.
     }
     // In column-major case, A, B, C, HSV are modified in place by the Fortran call.
 
 cleanup:
     /* --- Cleanup --- */
     free(dwork_allocated_buffer);
     free(iwork_allocated_buffer); 
     free(a_cm);
     free(b_cm);
     free(c_cm);
 
     // If info was SLICOT_MEMORY_ERROR from CHECK_ALLOC, it's already set.
     // Otherwise, info contains the Fortran routine's status.
     return info;
 }
