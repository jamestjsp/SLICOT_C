/**
 * @file ab07nd.c
 * @brief C wrapper implementation for SLICOT routine AB07ND
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB07ND,
 * which computes the inverse (Ai,Bi,Ci,Di) of a given system (A,B,C,D).
 * Refactored to align with ab01nd.c structure.
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

     /* Pointers to pass to Fortran */
     double *a_ptr, *b_ptr, *c_ptr, *d_ptr;
     int lda_f, ldb_f, ldc_f, ldd_f;


     /* --- Input Parameter Validation --- */

     if (n < 0) { info = -1; goto cleanup; }
     if (m < 0) { info = -2; goto cleanup; }

     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, n);
     int min_ldb_f = MAX(1, n);
     int min_ldc_f = MAX(1, m); // C is M-by-N
     int min_ldd_f = MAX(1, m); // D is M-by-M

     if (row_major) {
         // For row-major C, LDA/LDB/LDC/LDD is the number of columns
         if (n > 0 && lda < n) { info = -4; goto cleanup; }
         if (n > 0 && ldb < m) { info = -6; goto cleanup; } // B is N-by-M
         if (m > 0 && ldc < n) { info = -8; goto cleanup; } // C is M-by-N
         if (m > 0 && ldd < m) { info = -10; goto cleanup; } // D is M-by-M
     } else {
         // For column-major C, LDA/LDB/LDC/LDD is the number of rows (Fortran style)
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

     // Need temporary pointers for query if row_major, even though data isn't used yet
     // Use dummy values for LDs if dimensions are 0
     int lda_q = row_major ? MAX(1, n) : lda;
     int ldb_q = row_major ? MAX(1, n) : ldb;
     int ldc_q = row_major ? MAX(1, m) : ldc;
     int ldd_q = row_major ? MAX(1, m) : ldd;

     // Perform workspace query.
     // Note: The input arrays (a,b,c,d) might be modified by the query call
     // if the Fortran routine doesn't strictly adhere to input/output roles
     // during query. This wrapper assumes they are treated as input only for query.
     // If row_major, we technically should query using temporary CM arrays,
     // but since the query only depends on dimensions, querying with original
     // pointers and Fortran-style LDs should be okay.
     F77_FUNC(ab07nd, AB07ND)(&n, &m, a, &lda_q, b, &ldb_q, c, &ldc_q, d, &ldd_q,
                              rcond, iwork, &dwork_query, &ldwork, &info);

     if (info != 0 && info != (m + 1)) { // Allow INFO = M+1 during query (D singular)
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
         size_t a_size = (size_t)n * n; if (n == 0) a_size = 0;
         size_t b_size = (size_t)n * m; if (n == 0 || m == 0) b_size = 0; // N rows, M cols
         size_t c_size = (size_t)m * n; if (m == 0 || n == 0) c_size = 0; // M rows, N cols
         size_t d_size = (size_t)m * m; if (m == 0) d_size = 0; // M rows, M cols

         if (a_size > 0) { a_cm = (double*)malloc(a_size * sizeof(double)); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * sizeof(double)); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * sizeof(double)); CHECK_ALLOC(c_cm); }
         if (d_size > 0) { d_cm = (double*)malloc(d_size * sizeof(double)); CHECK_ALLOC(d_cm); }

         /* Transpose C (row-major) inputs to Fortran (column-major) copies */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, n, n, sizeof(double));
         if (b_size > 0) slicot_transpose_to_fortran(b, b_cm, n, m, sizeof(double)); // N rows, M cols
         if (c_size > 0) slicot_transpose_to_fortran(c, c_cm, m, n, sizeof(double)); // M rows, N cols
         if (d_size > 0) slicot_transpose_to_fortran(d, d_cm, m, m, sizeof(double)); // M rows, M cols

         /* Fortran leading dimensions (number of rows) */
         lda_f = MAX(1, n);
         ldb_f = MAX(1, n);
         ldc_f = MAX(1, m);
         ldd_f = MAX(1, m);

         /* Set pointers for Fortran call */
         a_ptr = a_cm; b_ptr = b_cm; c_ptr = c_cm; d_ptr = d_cm;

     } else {
         /* --- Column-Major Case --- */
         lda_f = lda; ldb_f = ldb; ldc_f = ldc; ldd_f = ldd;
         a_ptr = a; b_ptr = b; c_ptr = c; d_ptr = d;
     }

     /* Call the Fortran routine */
     F77_FUNC(ab07nd, AB07ND)(&n, &m,
                              a_ptr, &lda_f, b_ptr, &ldb_f,
                              c_ptr, &ldc_f, d_ptr, &ldd_f,
                              rcond,                     // Pass output pointer
                              iwork, dwork, &ldwork, &info); // Pass workspaces

     /* Copy back results if row_major */
     // Copy back even if info > 0 (singular), as calculation is completed
     if (row_major && info >= 0) {
         size_t a_size = (size_t)n * n; if (n == 0) a_size = 0;
         size_t b_size = (size_t)n * m; if (n == 0 || m == 0) b_size = 0; // N rows, M cols
         size_t c_size = (size_t)m * n; if (m == 0 || n == 0) c_size = 0; // M rows, N cols
         size_t d_size = (size_t)m * m; if (m == 0) d_size = 0; // M rows, M cols

         if (a_size > 0) slicot_transpose_to_c(a_cm, a, n, n, sizeof(double));
         if (b_size > 0) slicot_transpose_to_c(b_cm, b, n, m, sizeof(double)); // N rows, M cols
         if (c_size > 0) slicot_transpose_to_c(c_cm, c, m, n, sizeof(double)); // M rows, N cols
         if (d_size > 0) slicot_transpose_to_c(d_cm, d, m, m, sizeof(double)); // M rows, M cols
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
