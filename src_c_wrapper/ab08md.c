/**
 * @file ab08md.c
 * @brief C wrapper implementation for SLICOT routine AB08MD
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB08MD.
 * This routine computes the normal rank of the system pencil corresponding
 * to the state space system (A, B, C, D). It can optionally perform scaling.
 * Refactored to align with ab01nd.c structure.
 */

 #include <stdlib.h>
 #include <math.h>   // For isnan, isinf (optional checks)
 #include <ctype.h>  // For toupper
 #include <stddef.h> // For size_t

 // Include the header file for this wrapper
 #include "ab08md.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

 /*
  * Declare the external FORTRAN routine using the F77_FUNC macro.
  * Explicitly includes the hidden 'int' argument for the length of 'equil'.
  */
 extern void F77_FUNC(ab08md, AB08MD)(
     const char* equil, const int* n, const int* m, const int* p,
     double* a, const int* lda, double* b, const int* ldb,
     double* c, const int* ldc, double* d, const int* ldd,
     int* rank, const double* tol, int* iwork, double* dwork,
     const int* ldwork, int* info,
     int equil_len /* Hidden length argument for equil */);

 /* C wrapper function definition */
 SLICOT_C_WRAPPER_API
 int slicot_ab08md(char equil, int n, int m, int p,
                   double* a, int lda,
                   double* b, int ldb,
                   double* c, int ldc,
                   double* d, int ldd,
                   int* rank, double tol, int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = -1; /* Use -1 for workspace query */
     double dwork_query;
     double *dwork = NULL;
     int *iwork = NULL;
     int iwork_size = 0;
     const int equil_len = 1; // Fortran expects 1-based length for strings

     char equil_upper = toupper(equil);

     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL;

     /* Pointers to pass to Fortran */
     double *a_ptr, *b_ptr, *c_ptr, *d_ptr;
     int lda_f, ldb_f, ldc_f, ldd_f;

     /* --- Input Parameter Validation --- */

     // Check non-negative dimensions
     if (n < 0) { info = -2; goto cleanup; }
     if (m < 0) { info = -3; goto cleanup; }
     if (p < 0) { info = -4; goto cleanup; }

     // Check EQUIL character
     if (equil_upper != 'S' && equil_upper != 'N') {
         info = -1; goto cleanup;
     }
     // Optional: Check TOL range if necessary (Fortran routine might handle it)
     // if (tol < 0.0) { info = -14; goto cleanup; }

     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, n);
     int min_ldb_f = MAX(1, n);
     int min_ldc_f = MAX(1, p);
     int min_ldd_f = MAX(1, p);

     if (row_major) {
         // For row-major C, LDA is the number of columns
         int min_lda_rm_cols = n;
         int min_ldb_rm_cols = m;
         int min_ldc_rm_cols = n;
         int min_ldd_rm_cols = m;
         if (n > 0 && lda < min_lda_rm_cols) { info = -6; goto cleanup; }
         if (n > 0 && ldb < min_ldb_rm_cols) { info = -8; goto cleanup; }
         if (p > 0 && ldc < min_ldc_rm_cols) { info = -10; goto cleanup; }
         if (p > 0 && ldd < min_ldd_rm_cols) { info = -12; goto cleanup; }
     } else {
         // For column-major C, LDA is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -6; goto cleanup; }
         if (ldb < min_ldb_f) { info = -8; goto cleanup; }
         if (ldc < min_ldc_f) { info = -10; goto cleanup; }
         if (ldd < min_ldd_f) { info = -12; goto cleanup; }
     }

     /* --- Workspace Allocation --- */

     // Determine iwork size from AB08MD documentation: 2*N + MAX(M,P) + 1
     iwork_size = 2 * n + MAX(m, p) + 1;
     iwork_size = MAX(1, iwork_size); // Ensure minimum size 1
     iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
     CHECK_ALLOC(iwork);

     // Call Fortran routine for DWORK workspace query
     ldwork = -1; // Query mode
     // Use dummy LDs for query if dimensions are 0
     int lda_q = row_major ? MAX(1, n) : lda;
     int ldb_q = row_major ? MAX(1, n) : ldb;
     int ldc_q = row_major ? MAX(1, p) : ldc;
     int ldd_q = row_major ? MAX(1, p) : ldd;

     F77_FUNC(ab08md, AB08MD)(&equil_upper, &n, &m, &p,
                              NULL, &lda_q, NULL, &ldb_q, // NULL arrays for query
                              NULL, &ldc_q, NULL, &ldd_q,
                              rank, &tol,             // Pass pointers/addresses
                              iwork, &dwork_query,    // Pass iwork, address for query result
                              &ldwork, &info,         // Pass address of ldwork (-1), address of info
                              equil_len);             // Pass hidden length

     if (info != 0) {
         // Query failed, likely due to invalid N, M, P, LDA, etc. passed to query
         goto cleanup; // Go to cleanup to free iwork
     }

     // Get the required dwork size from query result
     ldwork = (int)dwork_query;
     // Check against minimum documented size: MAX(1, N + MAX(N, M, P))
     int min_ldwork = MAX(1, n + MAX(n, MAX(m, p)));
     ldwork = MAX(ldwork, min_ldwork);

     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure

     /* --- Prepare Arrays and Call Fortran Routine --- */
     size_t elem_size = sizeof(double);

     if (row_major) {
         /* --- Row-Major Case --- */

         /* Allocate memory for column-major copies */
         size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
         size_t b_rows = n; size_t b_cols = m; size_t b_size = b_rows * b_cols;
         size_t c_rows = p; size_t c_cols = n; size_t c_size = c_rows * c_cols;
         size_t d_rows = p; size_t d_cols = m; size_t d_size = d_rows * d_cols;

         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
         if (d_size > 0) { d_cm = (double*)malloc(d_size * elem_size); CHECK_ALLOC(d_cm); }

         /* Transpose C (row-major) inputs to Fortran (column-major) copies */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         if (b_size > 0) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, elem_size);
         if (c_size > 0) slicot_transpose_to_fortran(c, c_cm, c_rows, c_cols, elem_size);
         if (d_size > 0) slicot_transpose_to_fortran(d, d_cm, d_rows, d_cols, elem_size);

         /* Fortran leading dimensions (number of rows) */
         lda_f = MAX(1, a_rows);
         ldb_f = MAX(1, b_rows);
         ldc_f = MAX(1, c_rows);
         ldd_f = MAX(1, d_rows);

         /* Set pointers for Fortran call */
         a_ptr = a_cm; b_ptr = b_cm; c_ptr = c_cm; d_ptr = d_cm;

     } else {
         /* --- Column-Major Case --- */
         lda_f = lda; ldb_f = ldb; ldc_f = ldc; ldd_f = ldd;
         a_ptr = a; b_ptr = b; c_ptr = c; d_ptr = d;
     }

     /* Call the computational routine */
     F77_FUNC(ab08md, AB08MD)(&equil_upper, &n, &m, &p,
                              a_ptr, &lda_f, b_ptr, &ldb_f,
                              c_ptr, &ldc_f, d_ptr, &ldd_f,
                              rank, &tol,
                              iwork, dwork, &ldwork, &info,
                              equil_len);

     /* Copy back results if row_major and scaling was done (EQUIL='S') */
     if (row_major && info == 0 && equil_upper == 'S') {
         size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
         size_t b_rows = n; size_t b_cols = m; size_t b_size = b_rows * b_cols;
         size_t c_rows = p; size_t c_cols = n; size_t c_size = c_rows * c_cols;
         size_t d_rows = p; size_t d_cols = m; size_t d_size = d_rows * d_cols;

         if (a_size > 0) slicot_transpose_to_c(a_cm, a, a_rows, a_cols, elem_size);
         if (b_size > 0) slicot_transpose_to_c(b_cm, b, b_rows, b_cols, elem_size);
         if (c_size > 0) slicot_transpose_to_c(c_cm, c, c_rows, c_cols, elem_size);
         if (d_size > 0) slicot_transpose_to_c(d_cm, d, d_rows, d_cols, elem_size);
     }

 cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(iwork);
     free(a_cm);
     free(b_cm);
     free(c_cm);
     free(d_cm);

     /* Return the info code from the Fortran routine or SLICOT_MEMORY_ERROR */
     return info;
 }
