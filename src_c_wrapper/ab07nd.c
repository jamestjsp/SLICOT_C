/**
 * @file ab07nd.c
 * @brief C wrapper implementation for SLICOT routine AB07ND
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB07ND,
 * which computes the inverse (Ai,Bi,Ci,Di) of a given system (A,B,C,D).
 */

 #include <stdlib.h>
 #include <stddef.h> // For size_t
 
 // Include the header file for this wrapper
 #include "ab07nd.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions
 
 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * A, B, C, D are input/output.
  */
 extern void F77_FUNC(ab07nd, AB07ND)(
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     double* a,              // DOUBLE PRECISION A(LDA,*) (in/out)
     const int* lda,         // INTEGER LDA
     double* b,              // DOUBLE PRECISION B(LDB,*) (in/out)
     const int* ldb,         // INTEGER LDB
     double* c,              // DOUBLE PRECISION C(LDC,*) (in/out)
     const int* ldc,         // INTEGER LDC
     double* d,              // DOUBLE PRECISION D(LDD,*) (in/out)
     const int* ldd,         // INTEGER LDD
     double* rcond,          // DOUBLE PRECISION RCOND (output)
     int* iwork,             // INTEGER IWORK(*)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* info               // INTEGER INFO (output)
 );
 
 
 /* C wrapper function definition */
 SLICOT_C_WRAPPER_API
 int slicot_ab07nd(int n, int m,
                   double* a, int lda, double* b, int ldb,
                   double* c, int ldc, double* d, int ldd,
                   double* rcond, int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = -1; /* Use -1 for workspace query */
     double dwork_query;
     double* dwork = NULL;
     int* iwork = NULL;
     int iwork_size = 0;
 
     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL;
 
     /* --- Input Parameter Validation --- */
 
     if (n < 0) { info = -1; goto cleanup; }
     if (m < 0) { info = -2; goto cleanup; }
 
     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, n);
     int min_ldb_f = MAX(1, n);
     int min_ldc_f = MAX(1, m);
     int min_ldd_f = MAX(1, m);
 
     if (row_major) {
         // For row-major C, LDA is the number of columns
         int min_lda_rm_cols = n;
         int min_ldb_rm_cols = m;
         int min_ldc_rm_cols = n;
         int min_ldd_rm_cols = m;
         if (lda < min_lda_rm_cols) { info = -4; goto cleanup; }
         if (ldb < min_ldb_rm_cols) { info = -6; goto cleanup; }
         if (ldc < min_ldc_rm_cols) { info = -8; goto cleanup; }
         if (ldd < min_ldd_rm_cols) { info = -10; goto cleanup; }
     } else {
         // For column-major C, LDA is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -4; goto cleanup; }
         if (ldb < min_ldb_f) { info = -6; goto cleanup; }
         if (ldc < min_ldc_f) { info = -8; goto cleanup; }
         if (ldd < min_ldd_f) { info = -10; goto cleanup; }
     }
 
     /* --- Workspace Allocation --- */
 
     // Allocate IWORK (size 2*M)
     iwork_size = MAX(1, 2 * m); // Ensure minimum size 1
     iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
     CHECK_ALLOC(iwork);
 
     // Allocate DWORK based on query
     ldwork = -1; // Query mode
     F77_FUNC(ab07nd, AB07ND)(&n, &m, a, &lda, b, &ldb, c, &ldc, d, &ldd,
                              rcond, iwork, &dwork_query, &ldwork, &info);
 
     if (info != 0 && info != (m + 1)) { // Allow INFO = M+1 during query
         // Query failed for reasons other than numerical singularity warning
         goto cleanup;
     }
     info = 0; // Reset info if it was M+1 from query
 
     // Get the required dwork size from query result
     ldwork = (int)dwork_query;
     // Check against minimum documented size: MAX(1, 4*M)
     ldwork = MAX(ldwork, MAX(1, 4 * m));
 
     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure
 
     /* --- Prepare Arrays and Call Fortran Routine --- */
 
     if (row_major) {
         /* --- Row-Major Case --- */
 
         /* Allocate memory for column-major copies */
         size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
         size_t b_rows = n; size_t b_cols = m; size_t b_size = b_rows * b_cols;
         size_t c_rows = m; size_t c_cols = n; size_t c_size = c_rows * c_cols; // Rows=M
         size_t d_rows = m; size_t d_cols = m; size_t d_size = d_rows * d_cols; // Rows=M
 
         if (a_size > 0) { a_cm = (double*)malloc(a_size * sizeof(double)); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * sizeof(double)); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * sizeof(double)); CHECK_ALLOC(c_cm); }
         if (d_size > 0) { d_cm = (double*)malloc(d_size * sizeof(double)); CHECK_ALLOC(d_cm); }
 
         /* Transpose C (row-major) inputs to Fortran (column-major) copies */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, sizeof(double));
         if (b_size > 0) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, sizeof(double));
         if (c_size > 0) slicot_transpose_to_fortran(c, c_cm, c_rows, c_cols, sizeof(double)); // Rows=M
         if (d_size > 0) slicot_transpose_to_fortran(d, d_cm, d_rows, d_cols, sizeof(double)); // Rows=M
 
         /* Fortran leading dimensions */
         int lda_f = (a_rows > 0) ? a_rows : 1;
         int ldb_f = (b_rows > 0) ? b_rows : 1;
         int ldc_f = (c_rows > 0) ? c_rows : 1; // Rows=M
         int ldd_f = (d_rows > 0) ? d_rows : 1; // Rows=M
 
         /* Call the Fortran routine */
         F77_FUNC(ab07nd, AB07ND)(&n, &m,
                                  a_cm, &lda_f, b_cm, &ldb_f, // Pass CM arrays
                                  c_cm, &ldc_f, d_cm, &ldd_f,
                                  rcond,                     // Pass output pointer
                                  iwork, dwork, &ldwork, &info); // Pass workspaces
 
         /* Copy back results from column-major temps to original row-major arrays */
         // Copy back even if info > 0 (singular), as calculation is completed
         if (info >= 0) {
             if (a_size > 0) slicot_transpose_to_c(a_cm, a, a_rows, a_cols, sizeof(double));
             if (b_size > 0) slicot_transpose_to_c(b_cm, b, b_rows, b_cols, sizeof(double));
             if (c_size > 0) slicot_transpose_to_c(c_cm, c, c_rows, c_cols, sizeof(double)); // Rows=M
             if (d_size > 0) slicot_transpose_to_c(d_cm, d, d_rows, d_cols, sizeof(double)); // Rows=M
         }
         /* Column-major copies will be freed in cleanup */
 
     } else {
         /* --- Column-Major Case --- */
 
         /* Call the Fortran routine directly with user-provided arrays */
         F77_FUNC(ab07nd, AB07ND)(&n, &m,
                                  a, &lda, b, &ldb, // Pass original arrays
                                  c, &ldc, d, &ldd,
                                  rcond,             // Pass output pointer
                                  iwork, dwork, &ldwork, &info); // Pass workspaces
         // A, B, C, D are modified in place by the Fortran call.
     }
 
 cleanup:
     /* --- Cleanup --- */
     // Free allocated memory (free(NULL) is safe)
     free(dwork);
     free(iwork);
     // Free column-major copies if they were allocated
     free(a_cm);
     free(b_cm);
     free(c_cm);
     free(d_cm);
 
     /* Return the info code from the Fortran routine or SLICOT_MEMORY_ERROR */
     return info;
 }
 