/**
 * @file tf01rd.c
 * @brief C wrapper implementation for SLICOT routine TF01RD
 *
 * This file provides a C wrapper implementation for the SLICOT routine TF01RD,
 * which computes N Markov parameters M(k) = C*A^(k-1)*B for a
 * linear time-invariant system (A,B,C).
 */

 #include <stdlib.h>
 #include <string.h> // For memcpy
 #include <stddef.h> // For size_t

 // Include the header file for this wrapper
 // #include "tf01rd.h" // Assuming a header file exists
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, MIN, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  */
 extern void F77_FUNC(tf01rd, TF01RD)(
     const int* na,          // INTEGER NA
     const int* nb,          // INTEGER NB
     const int* nc,          // INTEGER NC
     const int* n,           // INTEGER N
     const double* a,        // DOUBLE PRECISION A(LDA,*) (in)
     const int* lda,         // INTEGER LDA
     const double* b,        // DOUBLE PRECISION B(LDB,*) (in)
     const int* ldb,         // INTEGER LDB
     const double* c,        // DOUBLE PRECISION C(LDC,*) (in)
     const int* ldc,         // INTEGER LDC
     double* h,              // DOUBLE PRECISION H(LDH,*) (output)
     const int* ldh,         // INTEGER LDH
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* info               // INTEGER INFO (output)
 );


 /* C wrapper function definition */
 SLICOT_EXPORT
 int slicot_tf01rd(int na, int nb, int nc, int n,
                   const double* a, int lda, const double* b, int ldb,
                   const double* c, int ldc, double* h, int ldh,
                   int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = 0; // Use calculated size
     double* dwork = NULL;
     int dwork_size = 0;
     // No iwork needed

     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL;
     double *h_cm = NULL; // For output
     const double *a_ptr, *b_ptr, *c_ptr; // Input pointers
     double *h_ptr; // Output pointer
     int lda_f, ldb_f, ldc_f, ldh_f;


     /* --- Input Parameter Validation --- */
     if (na < 0) { info = -1; goto cleanup; }
     if (nb < 0) { info = -2; goto cleanup; }
     if (nc < 0) { info = -3; goto cleanup; }
     if (n < 0) { info = -4; goto cleanup; }

     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, na);
     int min_ldb_f = MAX(1, na);
     int min_ldc_f = MAX(1, nc);
     int min_ldh_f = MAX(1, nc);

     if (row_major) {
         // For row-major C, LD is number of columns
         int min_lda_rm_cols = na;
         int min_ldb_rm_cols = nb;
         int min_ldc_rm_cols = na;
         int min_ldh_rm_cols = n * nb; // H is NC x (N*NB)
         if (lda < min_lda_rm_cols) { info = -6; goto cleanup; }
         if (ldb < min_ldb_rm_cols) { info = -8; goto cleanup; }
         if (ldc < min_ldc_rm_cols) { info = -10; goto cleanup; }
         if (ldh < min_ldh_rm_cols) { info = -12; goto cleanup; }
     } else {
         // For column-major C, LD is number of rows (Fortran style)
         if (lda < min_lda_f) { info = -6; goto cleanup; }
         if (ldb < min_ldb_f) { info = -8; goto cleanup; }
         if (ldc < min_ldc_f) { info = -10; goto cleanup; }
         if (ldh < min_ldh_f) { info = -12; goto cleanup; }
     }

     /* --- Prepare arrays for column-major format if using row-major --- */
     size_t elem_size = sizeof(double);
     if (row_major) {
         // Determine sizes for potential copies
         size_t a_rows = na; size_t a_cols = na; size_t a_size = a_rows * a_cols;
         size_t b_rows = na; size_t b_cols = nb; size_t b_size = b_rows * b_cols;
         size_t c_rows = nc; size_t c_cols = na; size_t c_size = c_rows * c_cols;
         size_t h_rows = nc; size_t h_cols = n * nb; size_t h_size = h_rows * h_cols;

         /* Allocate memory for column-major copies */
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
         if (h_size > 0) { h_cm = (double*)malloc(h_size * elem_size); CHECK_ALLOC(h_cm); }

         /* Transpose C inputs to Fortran copies */
         if (a_cm) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         if (b_cm) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, elem_size);
         if (c_cm) slicot_transpose_to_fortran(c, c_cm, c_rows, c_cols, elem_size);

         /* Fortran leading dimensions */
         lda_f = (na > 0) ? na : 1;
         ldb_f = (na > 0) ? na : 1;
         ldc_f = (nc > 0) ? nc : 1;
         ldh_f = (nc > 0) ? nc : 1;

         /* Set pointers */
         a_ptr = a_cm; b_ptr = b_cm; c_ptr = c_cm;
         h_ptr = h_cm;

     } else {
         /* Column-major case - use original arrays */
         lda_f = lda; ldb_f = ldb; ldc_f = ldc; ldh_f = ldh;
         a_ptr = a; b_ptr = b; c_ptr = c;
         h_ptr = h;
     }

     /* --- Workspace Allocation --- */
     // No query needed, size is fixed
     dwork_size = MAX(1, 2 * na * nc);
     ldwork = dwork_size; // Use calculated size directly
     dwork = (double*)malloc((size_t)dwork_size * sizeof(double));
     CHECK_ALLOC(dwork);

     /* --- Call the computational routine --- */
     F77_FUNC(tf01rd, TF01RD)(&na, &nb, &nc, &n,
                              a_ptr, &lda_f, b_ptr, &ldb_f, c_ptr, &ldc_f,
                              h_ptr, &ldh_f, dwork, &ldwork, &info);

     /* --- Copy results back to row-major format if needed --- */
     if (row_major && info == 0) {
         size_t h_rows = nc; size_t h_cols = n * nb; size_t h_size = h_rows * h_cols;
         if (h_size > 0 && h_cm) {
             slicot_transpose_to_c(h_cm, h, h_rows, h_cols, elem_size);
         }
     }
     // In column-major case, H modified in place.

  cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(a_cm); free(b_cm); free(c_cm);
     free(h_cm);

     return info;
 }