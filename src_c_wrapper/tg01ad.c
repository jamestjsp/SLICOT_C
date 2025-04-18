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
 #include "tg01ad.h"
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
 int slicot_tg01ad(char job, int l, int n, int m, int p, double thresh,
                   double* a, int lda, double* e, int lde,
                   double* b, int ldb, double* c, int ldc,
                   double* lscale, double* rscale, int row_major)
 {
     /* Local variables */
     int info = 0;
     double* dwork = NULL;
     int dwork_size = 0;
 
     const int job_len = 1;
     char job_upper = toupper(job);
 
     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *e_cm = NULL, *b_cm = NULL, *c_cm = NULL;
 
     /* --- Input Parameter Validation --- */
     if (l < 0) { info = -2; goto cleanup; }
     if (n < 0) { info = -3; goto cleanup; }
     if (m < 0) { info = -4; goto cleanup; }
     if (p < 0) { info = -5; goto cleanup; }
     if (thresh < 0.0) { info = -6; goto cleanup; }
     if (job_upper != 'A' && job_upper != 'B' && job_upper != 'C' && job_upper != 'N') {
         info = -1; goto cleanup;
     }
 
     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, l);
     int min_lde_f = MAX(1, l);
     int min_ldb_f = (m > 0) ? MAX(1, l) : 1;
     int min_ldc_f = MAX(1, p);
 
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
         if (ldb < min_ldb_f) { info = -12; goto cleanup; } // Check even if m=0? Fortran requires LDB>=1
         if (ldc < min_ldc_f) { info = -14; goto cleanup; }
     }
 
     /* --- Workspace Allocation --- */
     dwork_size = MAX(1, 3 * (l + n));
     dwork = (double*)malloc((size_t)dwork_size * sizeof(double));
     CHECK_ALLOC(dwork);
 
     /* --- Prepare Arrays and Call Fortran Routine --- */
     size_t elem_size = sizeof(double);
 
     // Determine sizes for potential copies
     size_t a_rows = l; size_t a_cols = n; size_t a_size = a_rows * a_cols;
     size_t e_rows = l; size_t e_cols = n; size_t e_size = e_rows * e_cols;
     size_t b_rows = l; size_t b_cols = m; size_t b_size = b_rows * b_cols;
     size_t c_rows = p; size_t c_cols = n; size_t c_size = c_rows * c_cols;
 
 
     if (row_major) {
         /* --- Row-Major Case --- */
 
         /* Allocate memory for column-major copies */
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (e_size > 0) { e_cm = (double*)malloc(e_size * elem_size); CHECK_ALLOC(e_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
 
         /* Transpose C inputs to Fortran copies */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         if (e_size > 0) slicot_transpose_to_fortran(e, e_cm, e_rows, e_cols, elem_size);
         if (b_size > 0) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, elem_size);
         if (c_size > 0) slicot_transpose_to_fortran(c, c_cm, c_rows, c_cols, elem_size);
 
         /* Fortran leading dimensions */
         int lda_f = (a_rows > 0) ? a_rows : 1;
         int lde_f = (e_rows > 0) ? e_rows : 1;
         int ldb_f = (b_rows > 0) ? b_rows : 1;
         int ldc_f = (c_rows > 0) ? c_rows : 1;
 
         /* Call the Fortran routine */
         F77_FUNC(tg01ad, TG01AD)(&job_upper, &l, &n, &m, &p, &thresh,
                                  a_cm, &lda_f, e_cm, &lde_f, b_cm, &ldb_f, c_cm, &ldc_f,
                                  lscale, rscale, dwork, &info,
                                  job_len);
 
         /* Copy back results from column-major temps to original row-major arrays */
         if (info == 0) {
             if (a_size > 0) slicot_transpose_to_c(a_cm, a, a_rows, a_cols, elem_size);
             if (e_size > 0) slicot_transpose_to_c(e_cm, e, e_rows, e_cols, elem_size);
             if (b_size > 0) slicot_transpose_to_c(b_cm, b, b_rows, b_cols, elem_size);
             if (c_size > 0) slicot_transpose_to_c(c_cm, c, c_rows, c_cols, elem_size);
             // LSCALE, RSCALE modified directly
         }
         /* Temps freed in cleanup */
 
     } else {
         /* --- Column-Major Case --- */
         /* Call the Fortran routine directly with user-provided arrays */
         F77_FUNC(tg01ad, TG01AD)(&job_upper, &l, &n, &m, &p, &thresh,
                                  a, &lda, e, &lde, b, &ldb, c, &ldc,
                                  lscale, rscale, dwork, &info,
                                  job_len);
         // A, E, B, C, LSCALE, RSCALE modified in place.
     }
 
  cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(a_cm);
     free(e_cm);
     free(b_cm);
     free(c_cm);
 
     return info;
 }
 