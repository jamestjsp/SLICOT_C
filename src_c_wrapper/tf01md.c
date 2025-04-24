/**
 * @file tf01md.c
 * @brief C wrapper implementation for SLICOT routine TF01MD
 *
 * This file provides a C wrapper implementation for the SLICOT routine TF01MD,
 * which computes the output response sequence of a linear time-invariant
 * discrete-time system given its state-space model (A,B,C,D) and
 * an input sequence.
 */

 #include <stdlib.h>
 #include <string.h> // For memcpy
 #include <stddef.h> // For size_t

 // Include the header file for this wrapper
 // #include "tf01md.h" // Assuming a header file exists
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, MIN, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  */
 extern void F77_FUNC(tf01md, TF01MD)(
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     const int* p,           // INTEGER P
     const int* ny,          // INTEGER NY
     const double* a,        // DOUBLE PRECISION A(LDA,*) (in)
     const int* lda,         // INTEGER LDA
     const double* b,        // DOUBLE PRECISION B(LDB,*) (in)
     const int* ldb,         // INTEGER LDB
     const double* c,        // DOUBLE PRECISION C(LDC,*) (in)
     const int* ldc,         // INTEGER LDC
     const double* d,        // DOUBLE PRECISION D(LDD,*) (in)
     const int* ldd,         // INTEGER LDD
     const double* u,        // DOUBLE PRECISION U(LDU,*) (in)
     const int* ldu,         // INTEGER LDU
     double* x,              // DOUBLE PRECISION X(*) (in/out)
     double* y,              // DOUBLE PRECISION Y(LDY,*) (output)
     const int* ldy,         // INTEGER LDY
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     int* info               // INTEGER INFO (output)
 );


 /* C wrapper function definition */
 SLICOT_EXPORT
 int slicot_tf01md(int n, int m, int p, int ny,
                   const double* a, int lda, const double* b, int ldb,
                   const double* c, int ldc, const double* d, int ldd,
                   const double* u, int ldu,
                   double* x, double* y, int ldy, int row_major)
 {
     /* Local variables */
     int info = 0;
     double* dwork = NULL;
     int dwork_size = 0;
     // No iwork needed

     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL, *u_cm = NULL;
     double *y_cm = NULL; // For output
     const double *a_ptr, *b_ptr, *c_ptr, *d_ptr, *u_ptr; // Input pointers
     double *y_ptr; // Output pointer
     int lda_f, ldb_f, ldc_f, ldd_f, ldu_f, ldy_f;

     /* --- Input Parameter Validation --- */
     if (n < 0) { info = -1; goto cleanup; }
     if (m < 0) { info = -2; goto cleanup; }
     if (p < 0) { info = -3; goto cleanup; }
     if (ny < 0) { info = -4; goto cleanup; }

     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, n);
     int min_ldb_f = MAX(1, n);
     int min_ldc_f = MAX(1, p);
     int min_ldd_f = MAX(1, p);
     int min_ldu_f = MAX(1, m);
     int min_ldy_f = MAX(1, p);

     if (row_major) {
         // For row-major C, LD is number of columns
         int min_lda_rm_cols = n;
         int min_ldb_rm_cols = m;
         int min_ldc_rm_cols = n;
         int min_ldd_rm_cols = m;
         int min_ldu_rm_cols = ny;
         int min_ldy_rm_cols = ny;
         if (lda < min_lda_rm_cols) { info = -6; goto cleanup; }
         if (ldb < min_ldb_rm_cols) { info = -8; goto cleanup; }
         if (ldc < min_ldc_rm_cols) { info = -10; goto cleanup; }
         if (ldd < min_ldd_rm_cols) { info = -12; goto cleanup; }
         if (ldu < min_ldu_rm_cols) { info = -14; goto cleanup; }
         if (ldy < min_ldy_rm_cols) { info = -17; goto cleanup; }
     } else {
         // For column-major C, LD is number of rows (Fortran style)
         if (lda < min_lda_f) { info = -6; goto cleanup; }
         if (ldb < min_ldb_f) { info = -8; goto cleanup; }
         if (ldc < min_ldc_f) { info = -10; goto cleanup; }
         if (ldd < min_ldd_f) { info = -12; goto cleanup; }
         if (ldu < min_ldu_f) { info = -14; goto cleanup; }
         if (ldy < min_ldy_f) { info = -17; goto cleanup; }
     }

     /* --- Prepare arrays for column-major format if using row-major --- */
     size_t elem_size = sizeof(double);
     if (row_major) {
         // Determine sizes for potential copies
         size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
         size_t b_rows = n; size_t b_cols = m; size_t b_size = b_rows * b_cols;
         size_t c_rows = p; size_t c_cols = n; size_t c_size = c_rows * c_cols;
         size_t d_rows = p; size_t d_cols = m; size_t d_size = d_rows * d_cols;
         size_t u_rows = m; size_t u_cols = ny; size_t u_size = u_rows * u_cols;
         size_t y_rows = p; size_t y_cols = ny; size_t y_size = y_rows * y_cols;

         /* Allocate memory for column-major copies */
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
         if (d_size > 0) { d_cm = (double*)malloc(d_size * elem_size); CHECK_ALLOC(d_cm); }
         if (u_size > 0) { u_cm = (double*)malloc(u_size * elem_size); CHECK_ALLOC(u_cm); }
         if (y_size > 0) { y_cm = (double*)malloc(y_size * elem_size); CHECK_ALLOC(y_cm); }

         /* Transpose C inputs to Fortran copies */
         if (a_cm) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         if (b_cm) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, elem_size);
         if (c_cm) slicot_transpose_to_fortran(c, c_cm, c_rows, c_cols, elem_size);
         if (d_cm) slicot_transpose_to_fortran(d, d_cm, d_rows, d_cols, elem_size);
         if (u_cm) slicot_transpose_to_fortran(u, u_cm, u_rows, u_cols, elem_size);

         /* Fortran leading dimensions */
         lda_f = (n > 0) ? n : 1;
         ldb_f = (n > 0) ? n : 1;
         ldc_f = (p > 0) ? p : 1;
         ldd_f = (p > 0) ? p : 1;
         ldu_f = (m > 0) ? m : 1;
         ldy_f = (p > 0) ? p : 1;

         /* Set pointers */
         a_ptr = a_cm; b_ptr = b_cm; c_ptr = c_cm; d_ptr = d_cm; u_ptr = u_cm;
         y_ptr = y_cm;

     } else {
         /* Column-major case - use original arrays */
         lda_f = lda; ldb_f = ldb; ldc_f = ldc; ldd_f = ldd; ldu_f = ldu; ldy_f = ldy;
         a_ptr = a; b_ptr = b; c_ptr = c; d_ptr = d; u_ptr = u;
         y_ptr = y;
     }

     /* --- Workspace Allocation --- */
     // No query needed, size is fixed
     dwork_size = n; // Size is N
     if (dwork_size < 1) dwork_size = 1; // Ensure at least 1
     dwork = (double*)malloc((size_t)dwork_size * sizeof(double));
     CHECK_ALLOC(dwork);

     /* --- Call the computational routine --- */
     // X is 1D in/out, no transposition needed
     F77_FUNC(tf01md, TF01MD)(&n, &m, &p, &ny,
                              a_ptr, &lda_f, b_ptr, &ldb_f, c_ptr, &ldc_f, d_ptr, &ldd_f,
                              u_ptr, &ldu_f, x, y_ptr, &ldy_f, dwork, &info);

     /* --- Copy results back to row-major format if needed --- */
     if (row_major && info == 0) {
         // X is modified directly (1D)
         size_t y_rows = p; size_t y_cols = ny; size_t y_size = y_rows * y_cols;
         if (y_size > 0 && y_cm) {
             slicot_transpose_to_c(y_cm, y, y_rows, y_cols, elem_size);
         }
     }
     // In column-major case, X, Y modified in place.

  cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(a_cm); free(b_cm); free(c_cm); free(d_cm); free(u_cm);
     free(y_cm);

     return info;
 }