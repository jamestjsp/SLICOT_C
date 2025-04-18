/**
 * @file tb01pd.c
 * @brief C wrapper implementation for SLICOT routine TB01PD
 *
 * This file provides a C wrapper implementation for the SLICOT routine TB01PD,
 * which finds a reduced (minimal, controllable, or observable)
 * state-space representation (Ar,Br,Cr) in upper block Hessenberg form
 * for a given system (A,B,C).
 */

 #include <stdlib.h>
 #include <ctype.h>
 #include <stddef.h> // For size_t
 
 // Include the header file for this wrapper
 #include "tb01pd.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions
 
 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * A, B, C are input/output. NR is output. IWORK contains output info.
  */
 extern void F77_FUNC(tb01pd, TB01PD)(
     const char* job,        // CHARACTER*1 JOB
     const char* equil,      // CHARACTER*1 EQUIL
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     const int* p,           // INTEGER P
     double* a,              // DOUBLE PRECISION A(LDA,*) (in/out)
     const int* lda,         // INTEGER LDA
     double* b,              // DOUBLE PRECISION B(LDB,*) (in/out)
     const int* ldb,         // INTEGER LDB
     double* c,              // DOUBLE PRECISION C(LDC,*) (in/out)
     const int* ldc,         // INTEGER LDC
     int* nr,                // INTEGER NR (output)
     const double* tol,      // DOUBLE PRECISION TOL
     int* iwork,             // INTEGER IWORK(*) (output info)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* info,              // INTEGER INFO (output)
     int job_len,            // Hidden length
     int equil_len           // Hidden length
 );
 
 
 /* C wrapper function definition */
 SLICOT_C_WRAPPER_API
 int slicot_tb01pd(char job, char equil, int n, int m, int p,
                   double* a, int lda, double* b, int ldb,
                   double* c, int ldc, int* nr, double tol,
                   int* iwork, int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = -1; /* Use -1 for workspace query */
     double dwork_query;
     double* dwork = NULL;
     const int job_len = 1, equil_len = 1;
 
     char job_upper = toupper(job);
     char equil_upper = toupper(equil);
 
     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL;
 
     /* --- Input Parameter Validation --- */
 
     if (n < 0) { info = -3; goto cleanup; }
     if (m < 0) { info = -4; goto cleanup; }
     if (p < 0) { info = -5; goto cleanup; }
     if (job_upper != 'M' && job_upper != 'C' && job_upper != 'O') { info = -1; goto cleanup; }
     if (equil_upper != 'S' && equil_upper != 'N') { info = -2; goto cleanup; }
     // TOL check done by Fortran
 
     // Check leading dimensions based on storage order and JOB
     int max_mp = MAX(m, p);
     int min_lda_f = MAX(1, n);
     int min_ldb_f = MAX(1, n);
     int min_ldc_f = (n > 0) ? MAX(1, max_mp) : 1; // LDC >= MAX(1,M,P) if N > 0
 
     int b_fort_cols = (job_upper == 'C') ? m : max_mp; // Effective columns for B workspace/input
     int c_fort_rows = max_mp; // Effective rows for C workspace/input
 
     if (row_major) {
         // For row-major C, LD is the number of columns
         int min_lda_rm_cols = n;
         int min_ldb_rm_cols = b_fort_cols;
         int min_ldc_rm_cols = n;
         if (lda < min_lda_rm_cols) { info = -7; goto cleanup; }
         if (ldb < min_ldb_rm_cols) { info = -9; goto cleanup; }
         if (ldc < min_ldc_rm_cols) { info = -11; goto cleanup; }
     } else {
         // For column-major C, LD is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -7; goto cleanup; }
         if (ldb < min_ldb_f) { info = -9; goto cleanup; }
         if (ldc < min_ldc_f) { info = -11; goto cleanup; }
     }
 
     /* --- Workspace Allocation --- */
 
     // Allocate DWORK based on query
     ldwork = -1; // Query mode
     F77_FUNC(tb01pd, TB01PD)(&job_upper, &equil_upper, &n, &m, &p,
                              a, &lda, b, &ldb, c, &ldc, nr, &tol,
                              iwork, &dwork_query, &ldwork, &info,
                              job_len, equil_len);
 
     if (info < 0 && info != -16) { info = info; goto cleanup; } // Query failed due to invalid argument
     info = 0; // Reset info after query
 
     // Get the required dwork size from query result
     ldwork = (int)dwork_query;
     // Check against minimum documented size: MAX(1, N + MAX(N, 3*M, 3*P))
     int min_ldwork = 1;
     min_ldwork = MAX(min_ldwork, n + MAX(n, MAX(3 * m, 3 * p)));
     ldwork = MAX(ldwork, min_ldwork);
 
     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure
 
     /* --- Prepare Arrays and Call Fortran Routine --- */
     size_t elem_size = sizeof(double);
 
     // Determine sizes for potential copies
     size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
     size_t b_rows = n; size_t b_cols = m; size_t b_size = b_rows * b_cols; // Actual B data size
     size_t c_rows = p; size_t c_cols = n; size_t c_size = c_rows * c_cols; // Actual C data size
 
     // Sizes for temporary arrays (might need larger workspace area in B/C)
     size_t b_work_cols = b_fort_cols; size_t b_work_size = b_rows * b_work_cols;
     size_t c_work_rows = c_fort_rows; size_t c_work_size = c_work_rows * c_cols;
 
 
     if (row_major) {
         /* --- Row-Major Case --- */
 
         /* Allocate memory for column-major copies */
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_work_size > 0) { b_cm = (double*)malloc(b_work_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_work_size > 0) { c_cm = (double*)malloc(c_work_size * elem_size); CHECK_ALLOC(c_cm); }
 
         /* Transpose C inputs to Fortran copies */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         if (b_size > 0) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, elem_size); // Copy only actual B part
         if (c_size > 0) slicot_transpose_to_fortran(c, c_cm, c_rows, c_cols, elem_size); // Copy only actual C part
 
         /* Fortran leading dimensions */
         int lda_f = (a_rows > 0) ? a_rows : 1;
         int ldb_f = (b_rows > 0) ? b_rows : 1;
         int ldc_f = (c_work_rows > 0) ? c_work_rows : 1; // Use workspace row dim for LDC
 
         /* Call the Fortran routine */
         F77_FUNC(tb01pd, TB01PD)(&job_upper, &equil_upper, &n, &m, &p,
                                  a_cm, &lda_f,           // Pass CM A
                                  b_cm, &ldb_f,           // Pass CM B (workspace)
                                  c_cm, &ldc_f,           // Pass CM C (workspace)
                                  nr, &tol, iwork, dwork, &ldwork, &info,
                                  job_len, equil_len);
 
         /* Copy back results from column-major temps to original row-major arrays */
         if (info == 0) {
             int nr_val = *nr;
             if (nr_val >= 0) { // Allow nr = 0
                 if (a_size > 0) slicot_transpose_to_c(a_cm, a, nr_val, nr_val, elem_size); // Ar
                 if (b_size > 0) slicot_transpose_to_c(b_cm, b, nr_val, m, elem_size);       // Br
                 if (c_size > 0) slicot_transpose_to_c(c_cm, c, p, nr_val, elem_size);       // Cr
             }
             // NR, IWORK modified directly
         }
         /* Temps freed in cleanup */
 
     } else {
         /* --- Column-Major Case --- */
 
         /* Call the Fortran routine directly with user-provided arrays */
         F77_FUNC(tb01pd, TB01PD)(&job_upper, &equil_upper, &n, &m, &p,
                                  a, &lda, b, &ldb, c, &ldc, nr, &tol,
                                  iwork, dwork, &ldwork, &info,
                                  job_len, equil_len);
         // A, B, C, NR, IWORK modified in place.
     }
 
 cleanup:
     /* --- Cleanup --- */
     free(dwork);
     // IWORK is output, not allocated here unless needed conditionally? Docs say output. Assume user provides it.
     free(a_cm);
     free(b_cm);
     free(c_cm);
 
     return info;
 }
 