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
 #include "tf01rd.h"
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
 int slicot_tf01rd(int na, int nb, int nc, int n,
                   const double* a, int lda, const double* b, int ldb,
                   const double* c, int ldc, double* h, int ldh,
                   int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = -1; /* Use -1 for workspace query (although size is fixed) */
     double dwork_query; // Dummy for query call consistency
     double* dwork = NULL;
     int dwork_size = 0;
 
     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL;
     double *h_cm = NULL; // For output
 
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
 
     /* --- Workspace Allocation --- */
     // Calculate DWORK size based on documentation
     dwork_size = MAX(1, 2 * na * nc);
     ldwork = dwork_size; // Use calculated size directly
     dwork = (double*)malloc((size_t)dwork_size * sizeof(double));
     CHECK_ALLOC(dwork);
 
     /* --- Prepare Arrays and Call Fortran Routine --- */
     size_t elem_size = sizeof(double);
 
     // Determine sizes for potential copies
     size_t a_rows = na; size_t a_cols = na; size_t a_size = a_rows * a_cols;
     size_t b_rows = na; size_t b_cols = nb; size_t b_size = b_rows * b_cols;
     size_t c_rows = nc; size_t c_cols = na; size_t c_size = c_rows * c_cols;
     size_t h_rows = nc; size_t h_cols = n * nb; size_t h_size = h_rows * h_cols;
 
 
     if (row_major) {
         /* --- Row-Major Case --- */
 
         /* Allocate memory for column-major copies */
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
         if (h_size > 0) { h_cm = (double*)malloc(h_size * elem_size); CHECK_ALLOC(h_cm); }
 
         /* Transpose C inputs to Fortran copies */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         if (b_size > 0) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, elem_size);
         if (c_size > 0) slicot_transpose_to_fortran(c, c_cm, c_rows, c_cols, elem_size);
 
         /* Fortran leading dimensions */
         int lda_f = (a_rows > 0) ? a_rows : 1;
         int ldb_f = (b_rows > 0) ? b_rows : 1;
         int ldc_f = (c_rows > 0) ? c_rows : 1;
         int ldh_f = (h_rows > 0) ? h_rows : 1;
 
         /* Call the Fortran routine */
         F77_FUNC(tf01rd, TF01RD)(&na, &nb, &nc, &n,
                                  a_cm, &lda_f, b_cm, &ldb_f, c_cm, &ldc_f,
                                  h_cm, &ldh_f, dwork, &ldwork, &info);
 
         /* Copy back results from column-major temps to original row-major arrays */
         if (info == 0) {
             if (h_size > 0) slicot_transpose_to_c(h_cm, h, h_rows, h_cols, elem_size);
         }
         /* Temps freed in cleanup */
 
     } else {
         /* --- Column-Major Case --- */
         /* Call the Fortran routine directly with user-provided arrays */
         F77_FUNC(tf01rd, TF01RD)(&na, &nb, &nc, &n,
                                  a, &lda, b, &ldb, c, &ldc,
                                  h, &ldh, dwork, &ldwork, &info);
         // H modified in place.
     }
 
  cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(a_cm);
     free(b_cm);
     free(c_cm);
     free(h_cm);
 
     return info;
 }
 