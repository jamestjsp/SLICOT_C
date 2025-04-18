/**
 * @file tb01id.c
 * @brief C wrapper implementation for SLICOT routine TB01ID
 *
 * This file provides a C wrapper implementation for the SLICOT routine TB01ID,
 * which balances a system matrix S = [A B; C 0] corresponding to a
 * state-space triplet (A,B,C) using diagonal similarity transformations.
 */

 #include <stdlib.h>
 #include <ctype.h>
 #include <stddef.h> // For size_t
 
 // Include the header file for this wrapper
 #include "tb01id.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions
 
 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * MAXRED, A, B, C are input/output. SCALE is output.
  */
 extern void F77_FUNC(tb01id, TB01ID)(
     const char* job,        // CHARACTER*1 JOB
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     const int* p,           // INTEGER P
     double* maxred,         // DOUBLE PRECISION MAXRED (in/out)
     double* a,              // DOUBLE PRECISION A(LDA,*) (in/out)
     const int* lda,         // INTEGER LDA
     double* b,              // DOUBLE PRECISION B(LDB,*) (in/out)
     const int* ldb,         // INTEGER LDB
     double* c,              // DOUBLE PRECISION C(LDC,*) (in/out)
     const int* ldc,         // INTEGER LDC
     double* scale,          // DOUBLE PRECISION SCALE(*) (output)
     int* info,              // INTEGER INFO (output)
     int job_len             // Hidden length
 );
 
 
 /* C wrapper function definition */
 SLICOT_C_WRAPPER_API
 int slicot_tb01id(char job, int n, int m, int p, double* maxred,
                   double* a, int lda, double* b, int ldb,
                   double* c, int ldc, double* scale,
                   int row_major)
 {
     /* Local variables */
     int info = 0;
     const int job_len = 1;
     char job_upper = toupper(job);
 
     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL;
 
     /* --- Input Parameter Validation --- */
 
     if (n < 0) { info = -2; goto cleanup; }
     if (m < 0) { info = -3; goto cleanup; }
     if (p < 0) { info = -4; goto cleanup; }
     if (job_upper != 'A' && job_upper != 'B' && job_upper != 'C' && job_upper != 'N') {
         info = -1; goto cleanup;
     }
     // MAXRED check done by Fortran
 
     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, n);
     int min_ldb_f = (m > 0) ? MAX(1, n) : 1;
     int min_ldc_f = MAX(1, p);
 
     if (row_major) {
         // For row-major C, LD is the number of columns
         int min_lda_rm_cols = n;
         int min_ldb_rm_cols = (m > 0) ? m : 1;
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
     // No external workspace needed for this routine
 
     /* --- Prepare Arrays and Call Fortran Routine --- */
     size_t elem_size = sizeof(double);
 
     // Determine sizes for potential copies
     size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
     size_t b_rows = n; size_t b_cols = m; size_t b_size = b_rows * b_cols;
     size_t c_rows = p; size_t c_cols = n; size_t c_size = c_rows * c_cols;
 
     if (row_major) {
         /* --- Row-Major Case --- */
 
         /* Allocate memory for column-major copies */
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (m > 0 && b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (p > 0 && c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
 
         /* Transpose C inputs to Fortran copies */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         if (m > 0 && b_size > 0) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, elem_size);
         if (p > 0 && c_size > 0) slicot_transpose_to_fortran(c, c_cm, c_rows, c_cols, elem_size);
 
         /* Fortran leading dimensions */
         int lda_f = (a_rows > 0) ? a_rows : 1;
         int ldb_f = (b_rows > 0) ? b_rows : 1;
         int ldc_f = (c_rows > 0) ? c_rows : 1;
 
         /* Call the Fortran routine */
         F77_FUNC(tb01id, TB01ID)(&job_upper, &n, &m, &p, maxred,
                                  a_cm, &lda_f,
                                  (m > 0 ? b_cm : NULL), &ldb_f,
                                  (p > 0 ? c_cm : NULL), &ldc_f,
                                  scale, &info, job_len);
 
         /* Copy back results from column-major temps to original row-major arrays */
         if (info == 0) {
             if (a_size > 0) slicot_transpose_to_c(a_cm, a, a_rows, a_cols, elem_size);
             if (m > 0 && b_size > 0) slicot_transpose_to_c(b_cm, b, b_rows, b_cols, elem_size);
             if (p > 0 && c_size > 0) slicot_transpose_to_c(c_cm, c, c_rows, c_cols, elem_size);
             // MAXRED, SCALE modified directly
         }
         /* Temps freed in cleanup */
 
     } else {
         /* --- Column-Major Case --- */
 
         /* Call the Fortran routine directly with user-provided arrays */
         F77_FUNC(tb01id, TB01ID)(&job_upper, &n, &m, &p, maxred,
                                  a, &lda, b, &ldb, c, &ldc,
                                  scale, &info, job_len);
         // MAXRED, A, B, C, SCALE modified in place.
     }
 
 cleanup:
     /* --- Cleanup --- */
     free(a_cm);
     free(b_cm);
     free(c_cm);
 
     return info;
 }
 