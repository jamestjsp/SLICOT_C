/**
 * @file tg01ad.c
 * @brief C wrapper implementation for SLICOT routine TG01AD
 *
 * This file provides a C wrapper implementation for the SLICOT routine TG01AD,
 * which balances the matrices of the system pencil corresponding to a
 * descriptor triple (A-lambda E,B,C) using diagonal similarity
 * transformations.
 */

 #include <stdlib.h>
 #include <ctype.h>
 #include <string.h> // For memcpy
 #include <stddef.h> // For size_t

 // Include the header file for this wrapper
 // #include "tg01ad.h" // Assuming a header file exists
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, MIN, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * Hidden length for CHARACTER argument is added at the end.
  */
 extern void F77_FUNC(tg01ad, TG01AD)(
     const char* job,        // CHARACTER*1 JOB
     const int* l,           // INTEGER L
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     const int* p,           // INTEGER P
     const double* thresh,   // DOUBLE PRECISION THRESH
     double* a,              // DOUBLE PRECISION A(LDA,*) (in/out)
     const int* lda,         // INTEGER LDA
     double* e,              // DOUBLE PRECISION E(LDE,*) (in/out)
     const int* lde,         // INTEGER LDE
     double* b,              // DOUBLE PRECISION B(LDB,*) (in/out)
     const int* ldb,         // INTEGER LDB
     double* c,              // DOUBLE PRECISION C(LDC,*) (in/out)
     const int* ldc,         // INTEGER LDC
     double* lscale,         // DOUBLE PRECISION LSCALE(*) (output)
     double* rscale,         // DOUBLE PRECISION RSCALE(*) (output)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     int* info,              // INTEGER INFO (output)
     int job_len             // Hidden length
 );


 /* C wrapper function definition */
 SLICOT_C_WRAPPER_API
 int slicot_tg01ad(char job, int l, int n, int m, int p, double thresh,
                   double* a, int lda, double* e, int lde,
                   double* b, int ldb, double* c, int ldc,
                   double* lscale, double* rscale, int row_major)
 {
     /* Local variables */
     int info = 0;
     double* dwork = NULL;
     int dwork_size = 0;
     int ldwork = 0; // Use calculated size
     // No iwork needed

     const int job_len = 1;
     char job_upper = toupper(job);

     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *e_cm = NULL, *b_cm = NULL, *c_cm = NULL;
     double *a_ptr, *e_ptr, *b_ptr, *c_ptr;
     int lda_f, lde_f, ldb_f, ldc_f;

     /* --- Input Parameter Validation --- */
     if (job_upper != 'A' && job_upper != 'B' && job_upper != 'C' && job_upper != 'N') {
         info = -1; goto cleanup;
     }
     if (l < 0) { info = -2; goto cleanup; }
     if (n < 0) { info = -3; goto cleanup; }
     if (m < 0) { info = -4; goto cleanup; }
     if (p < 0) { info = -5; goto cleanup; }
     if (thresh < 0.0) { info = -6; goto cleanup; }

     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, l);
     int min_lde_f = MAX(1, l);
     int min_ldb_f = (m > 0) ? MAX(1, l) : 1; // B is L x M
     int min_ldc_f = MAX(1, p);              // C is P x N

     if (row_major) {
         // For row-major C, LD is number of columns
         int min_lda_rm_cols = n;
         int min_lde_rm_cols = n;
         int min_ldb_rm_cols = m;
         int min_ldc_rm_cols = n;
         if (lda < min_lda_rm_cols) { info = -8; goto cleanup; }
         if (lde < min_lde_rm_cols) { info = -10; goto cleanup; }
         if (m > 0 && ldb < min_ldb_rm_cols) { info = -12; goto cleanup; }
         if (ldc < min_ldc_rm_cols) { info = -14; goto cleanup; }
     } else {
         // For column-major C, LD is number of rows (Fortran style)
         if (lda < min_lda_f) { info = -8; goto cleanup; }
         if (lde < min_lde_f) { info = -10; goto cleanup; }
         if (ldb < min_ldb_f) { info = -12; goto cleanup; } // Check even if m=0
         if (ldc < min_ldc_f) { info = -14; goto cleanup; }
     }

     /* --- Prepare arrays for column-major format if using row-major --- */
     size_t elem_size = sizeof(double);
     if (row_major) {
         // Determine sizes for potential copies
         size_t a_rows = l; size_t a_cols = n; size_t a_size = a_rows * a_cols;
         size_t e_rows = l; size_t e_cols = n; size_t e_size = e_rows * e_cols;
         size_t b_rows = l; size_t b_cols = m; size_t b_size = b_rows * b_cols;
         size_t c_rows = p; size_t c_cols = n; size_t c_size = c_rows * c_cols;

         /* Allocate memory for column-major copies */
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (e_size > 0) { e_cm = (double*)malloc(e_size * elem_size); CHECK_ALLOC(e_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }

         /* Transpose C inputs to Fortran copies */
         if (a_cm) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         if (e_cm) slicot_transpose_to_fortran(e, e_cm, e_rows, e_cols, elem_size);
         if (b_cm) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, elem_size);
         if (c_cm) slicot_transpose_to_fortran(c, c_cm, c_rows, c_cols, elem_size);

         /* Fortran leading dimensions */
         lda_f = (l > 0) ? l : 1;
         lde_f = (l > 0) ? l : 1;
         ldb_f = (l > 0) ? l : 1;
         ldc_f = (p > 0) ? p : 1;

         /* Set pointers */
         a_ptr = a_cm; e_ptr = e_cm; b_ptr = b_cm; c_ptr = c_cm;

     } else {
         /* Column-major case - use original arrays */
         lda_f = lda; lde_f = lde; ldb_f = ldb; ldc_f = ldc;
         a_ptr = a; e_ptr = e; b_ptr = b; c_ptr = c;
     }

     /* --- Workspace Allocation --- */
     // No query needed, size is fixed
     dwork_size = MAX(1, 3 * (l + n));
     ldwork = dwork_size; // Use calculated size directly
     dwork = (double*)malloc((size_t)dwork_size * sizeof(double));
     CHECK_ALLOC(dwork);

     /* --- Call the computational routine --- */
     // LSCALE, RSCALE are 1D outputs, no transposition needed
     F77_FUNC(tg01ad, TG01AD)(&job_upper, &l, &n, &m, &p, &thresh,
                              a_ptr, &lda_f, e_ptr, &lde_f, b_ptr, &ldb_f, c_ptr, &ldc_f,
                              lscale, rscale, dwork, &info,
                              job_len);

     /* --- Copy results back to row-major format if needed --- */
     if (row_major && info == 0) {
         size_t a_rows = l; size_t a_cols = n; size_t a_size = a_rows * a_cols;
         size_t e_rows = l; size_t e_cols = n; size_t e_size = e_rows * e_cols;
         size_t b_rows = l; size_t b_cols = m; size_t b_size = b_rows * b_cols;
         size_t c_rows = p; size_t c_cols = n; size_t c_size = c_rows * c_cols;

         if (a_size > 0 && a_cm) slicot_transpose_to_c(a_cm, a, a_rows, a_cols, elem_size);
         if (e_size > 0 && e_cm) slicot_transpose_to_c(e_cm, e, e_rows, e_cols, elem_size);
         if (b_size > 0 && b_cm) slicot_transpose_to_c(b_cm, b, b_rows, b_cols, elem_size);
         if (c_size > 0 && c_cm) slicot_transpose_to_c(c_cm, c, c_rows, c_cols, elem_size);
         // LSCALE, RSCALE modified directly
     }
     // In column-major case, A, E, B, C, LSCALE, RSCALE modified in place.

  cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(a_cm); free(e_cm); free(b_cm); free(c_cm);

     return info;
 }