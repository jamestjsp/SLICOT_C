/**
 * @file ab08md.c
 * @brief C wrapper implementation for SLICOT routine AB08MD
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB08MD.
 * This routine computes the normal rank of the system pencil corresponding
 * to the state space system (A, B, C, D). It can optionally perform scaling.
 *
 * Note: Explicitly passes Fortran hidden string length argument for 'equil'.
 */

 #include <stdlib.h>
 #include <math.h>
 #include <ctype.h> // For toupper
 
 #include "ab08md.h"         // Header for this wrapper function
 #include "slicot_utils.h"   // For MAX, SLICOT_MEMORY_ERROR, CHECK_ALLOC, transpose routines
 #include "slicot_f77.h"     // For F77_FUNC macro conventions
 
 /*
  * Declare the external FORTRAN routine using the F77_FUNC macro
  * - Handles potential compiler name mangling.
  * - Explicitly includes the hidden 'int' argument at the end
  * for the length of the 'equil' CHARACTER*1 argument.
  */
 extern void F77_FUNC(ab08md, AB08MD)(
     const char* equil, const int* n, const int* m, const int* p,
     double* a, const int* lda, double* b, const int* ldb,
     double* c, const int* ldc, double* d, const int* ldd,
     int* rank, const double* tol, int* iwork, double* dwork,
     const int* ldwork, int* info,
     int equil_len /* Hidden length argument for equil */);
 
 /* C wrapper function definition */
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
     int equil_len = 1; // Fortran expects 1-based length for strings
 
     char equil_upper = toupper(equil);
 
     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL;
 
     /* --- Input Parameter Validation --- */
 
     // Check non-negative dimensions
     if (n < 0) { info = -2; goto cleanup; }
     if (m < 0) { info = -3; goto cleanup; }
     if (p < 0) { info = -4; goto cleanup; }
 
     // Check EQUIL character
     if (equil_upper != 'S' && equil_upper != 'N') {
         info = -1; goto cleanup;
     }
 
     // Check leading dimensions based on storage order
     // Fortran minimum leading dimensions (number of rows)
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
         if (lda < min_lda_rm_cols) { info = -6; goto cleanup; }
         if (ldb < min_ldb_rm_cols) { info = -8; goto cleanup; }
         if (ldc < min_ldc_rm_cols) { info = -10; goto cleanup; }
         if (ldd < min_ldd_rm_cols) { info = -12; goto cleanup; }
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
     if (iwork_size < 1) iwork_size = 1; // Ensure minimum size 1
     iwork = (int*)malloc(iwork_size * sizeof(int));
     // Use CHECK_ALLOC macro (assumed defined in slicot_utils.h) for error handling
     // It sets info = SLICOT_MEMORY_ERROR and jumps to cleanup on failure.
     CHECK_ALLOC(iwork);
 
     // Call Fortran routine for DWORK workspace query
     ldwork = -1; // Query mode
     F77_FUNC(ab08md, AB08MD)(&equil_upper, &n, &m, &p,
                              NULL, &lda, NULL, &ldb, // NULL arrays for query
                              NULL, &ldc, NULL, &ldd,
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
     if (ldwork < 1) ldwork = 1; // Ensure minimum size 1
     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure
 
     /* --- Prepare Arrays and Call Fortran Routine --- */
 
     if (row_major) {
         /* --- Row-Major Case --- */
 
         /* Allocate memory for column-major copies */
         size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
         size_t b_rows = n; size_t b_cols = m; size_t b_size = b_rows * b_cols;
         size_t c_rows = p; size_t c_cols = n; size_t c_size = c_rows * c_cols;
         size_t d_rows = p; size_t d_cols = m; size_t d_size = d_rows * d_cols;
 
         if (a_size > 0) { a_cm = (double*)malloc(a_size * sizeof(double)); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * sizeof(double)); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * sizeof(double)); CHECK_ALLOC(c_cm); }
         if (d_size > 0) { d_cm = (double*)malloc(d_size * sizeof(double)); CHECK_ALLOC(d_cm); }
 
         /* Transpose C (row-major) inputs to Fortran (column-major) copies */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, sizeof(double));
         if (b_size > 0) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, sizeof(double));
         if (c_size > 0) slicot_transpose_to_fortran(c, c_cm, c_rows, c_cols, sizeof(double));
         if (d_size > 0) slicot_transpose_to_fortran(d, d_cm, d_rows, d_cols, sizeof(double));
 
         /* Fortran leading dimensions (number of rows in conceptual Fortran array) */
         int lda_f = (a_rows > 0) ? a_rows : 1;
         int ldb_f = (b_rows > 0) ? b_rows : 1;
         int ldc_f = (c_rows > 0) ? c_rows : 1;
         int ldd_f = (d_rows > 0) ? d_rows : 1;
 
         /* Call the Fortran routine with column-major copies */
         F77_FUNC(ab08md, AB08MD)(&equil_upper, &n, &m, &p,
                                  a_cm, &lda_f, b_cm, &ldb_f, // Pass CM arrays and Fortran LDs
                                  c_cm, &ldc_f, d_cm, &ldd_f,
                                  rank, &tol,                 // Pass original pointers/addresses
                                  iwork, dwork, &ldwork, &info, // Pass workspaces
                                  equil_len);                   // Pass hidden length
 
         /* Copy back results if scaling was done (EQUIL='S') */
         /* Note: AB08MD documentation does not explicitly state A,B,C,D are output */
         /* when EQUIL='S'. This copy-back is included based on the non-const */
         /* function signature allowing modification, for safety. */
         if (info == 0 && equil_upper == 'S') {
              if (a_size > 0) slicot_transpose_to_c(a_cm, a, a_rows, a_cols, sizeof(double));
              if (b_size > 0) slicot_transpose_to_c(b_cm, b, b_rows, b_cols, sizeof(double));
              if (c_size > 0) slicot_transpose_to_c(c_cm, c, c_rows, c_cols, sizeof(double));
              if (d_size > 0) slicot_transpose_to_c(d_cm, d, d_rows, d_cols, sizeof(double));
         }
         /* Column-major copies will be freed in cleanup */
 
     } else {
         /* --- Column-Major Case --- */
 
         /* Call the Fortran routine directly with user-provided arrays */
         F77_FUNC(ab08md, AB08MD)(&equil_upper, &n, &m, &p,
                                  a, &lda, b, &ldb, // Pass original arrays and LDs
                                  c, &ldc, d, &ldd,
                                  rank, &tol,                 // Pass original pointers/addresses
                                  iwork, dwork, &ldwork, &info, // Pass workspaces
                                  equil_len);                   // Pass hidden length
     }
 
 cleanup:
     /* --- Cleanup --- */
     // Free allocated workspace memory
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
 