/**
 * @file sb04qd.c
 * @brief C wrapper implementation for SLICOT routine SB04QD
 *
 * This file provides a C wrapper implementation for the SLICOT routine SB04QD,
 * which solves the discrete-time Sylvester equation X + AXB = C
 * using the Hessenberg-Schur method.
 */

 #include <stdlib.h>
 #include <stddef.h> // For size_t
 
 // Include the header file for this wrapper
 #include "sb04qd.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions
 
 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * A, B, C are input/output. Z is output.
  */
 extern void F77_FUNC(sb04qd, SB04QD)(
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     double* a,              // DOUBLE PRECISION A(LDA,*) (in/out)
     const int* lda,         // INTEGER LDA
     double* b,              // DOUBLE PRECISION B(LDB,*) (in/out)
     const int* ldb,         // INTEGER LDB
     double* c,              // DOUBLE PRECISION C(LDC,*) (in/out -> X)
     const int* ldc,         // INTEGER LDC
     double* z,              // DOUBLE PRECISION Z(LDZ,*) (output)
     const int* ldz,         // INTEGER LDZ
     int* iwork,             // INTEGER IWORK(*)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* info               // INTEGER INFO (output)
 );
 
 
 /* C wrapper function definition */
 SLICOT_C_WRAPPER_API
 int slicot_sb04qd(int n, int m,
                   double* a, int lda, double* b, int ldb,
                   double* c, int ldc, double* z, int ldz,
                   int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = -1; /* Use -1 for workspace query */
     double dwork_query;
     double* dwork = NULL;
     int* iwork = NULL;
     int iwork_size = 0;
 
     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *z_cm = NULL;
 
     /* --- Input Parameter Validation --- */
 
     if (n < 0) { info = -1; goto cleanup; }
     if (m < 0) { info = -2; goto cleanup; }
 
     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, n);
     int min_ldb_f = MAX(1, m);
     int min_ldc_f = MAX(1, n);
     int min_ldz_f = MAX(1, m);
 
     if (row_major) {
         // For row-major C, LD is the number of columns
         int min_lda_rm_cols = n;
         int min_ldb_rm_cols = m;
         int min_ldc_rm_cols = m;
         int min_ldz_rm_cols = m;
         if (lda < min_lda_rm_cols) { info = -4; goto cleanup; }
         if (ldb < min_ldb_rm_cols) { info = -6; goto cleanup; }
         if (ldc < min_ldc_rm_cols) { info = -8; goto cleanup; }
         if (ldz < min_ldz_rm_cols) { info = -10; goto cleanup; }
     } else {
         // For column-major C, LD is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -4; goto cleanup; }
         if (ldb < min_ldb_f) { info = -6; goto cleanup; }
         if (ldc < min_ldc_f) { info = -8; goto cleanup; }
         if (ldz < min_ldz_f) { info = -10; goto cleanup; }
     }
 
     /* --- Workspace Allocation --- */
 
     // Allocate IWORK (size 4*N)
     iwork_size = MAX(1, 4 * n);
     iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
     CHECK_ALLOC(iwork);
 
     // Allocate DWORK based on query
     ldwork = -1; // Query mode
     F77_FUNC(sb04qd, SB04QD)(&n, &m, a, &lda, b, &ldb, c, &ldc, z, &ldz,
                              iwork, &dwork_query, &ldwork, &info);
 
     if (info < 0) { goto cleanup; } // Query failed due to invalid argument
     info = 0; // Reset info after query
 
     // Get the required dwork size from query result
     ldwork = (int)dwork_query;
     // Check against minimum documented size: MAX(1, 2*N*N + 9*N, 5*M, N + M)
     int min_ldwork = 1;
     min_ldwork = MAX(min_ldwork, 2 * n * n + 9 * n);
     min_ldwork = MAX(min_ldwork, 5 * m);
     min_ldwork = MAX(min_ldwork, n + m);
     ldwork = MAX(ldwork, min_ldwork);
 
     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure
 
     /* --- Prepare Arrays and Call Fortran Routine --- */
     size_t elem_size = sizeof(double);
 
     // Determine sizes for potential copies
     size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
     size_t b_rows = m; size_t b_cols = m; size_t b_size = b_rows * b_cols;
     size_t c_rows = n; size_t c_cols = m; size_t c_size = c_rows * c_cols;
     size_t z_rows = m; size_t z_cols = m; size_t z_size = z_rows * z_cols;
 
 
     if (row_major) {
         /* --- Row-Major Case --- */
 
         /* Allocate memory for column-major copies */
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
         if (z_size > 0) { z_cm = (double*)malloc(z_size * elem_size); CHECK_ALLOC(z_cm); }
 
         /* Transpose C inputs to Fortran copies */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         if (b_size > 0) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, elem_size);
         if (c_size > 0) slicot_transpose_to_fortran(c, c_cm, c_rows, c_cols, elem_size);
 
         /* Fortran leading dimensions */
         int lda_f = (a_rows > 0) ? a_rows : 1;
         int ldb_f = (b_rows > 0) ? b_rows : 1;
         int ldc_f = (c_rows > 0) ? c_rows : 1;
         int ldz_f = (z_rows > 0) ? z_rows : 1;
 
         /* Call the Fortran routine */
         F77_FUNC(sb04qd, SB04QD)(&n, &m,
                                  a_cm, &lda_f,           // Pass CM A (in/out)
                                  b_cm, &ldb_f,           // Pass CM B (in/out)
                                  c_cm, &ldc_f,           // Pass CM C (in/out -> X)
                                  z_cm, &ldz_f,           // Pass CM Z (out)
                                  iwork, dwork, &ldwork, &info);
 
         /* Copy back results from column-major temps to original row-major arrays */
         if (info == 0 || info > m) { // Copy back even if INFO > M (singular solve)
             if (a_size > 0) slicot_transpose_to_c(a_cm, a, a_rows, a_cols, elem_size); // Modified A
             if (b_size > 0) slicot_transpose_to_c(b_cm, b, b_rows, b_cols, elem_size); // Modified B
             if (c_size > 0) slicot_transpose_to_c(c_cm, c, c_rows, c_cols, elem_size); // Solution X
             if (z_size > 0) slicot_transpose_to_c(z_cm, z, z_rows, z_cols, elem_size); // Output Z
         }
         /* Temps freed in cleanup */
 
     } else {
         /* --- Column-Major Case --- */
 
         /* Call the Fortran routine directly with user-provided arrays */
         F77_FUNC(sb04qd, SB04QD)(&n, &m, a, &lda, b, &ldb, c, &ldc, z, &ldz,
                                  iwork, dwork, &ldwork, &info);
         // A, B, C, Z are modified in place.
     }
 
 cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(iwork);
     free(a_cm);
     free(b_cm);
     free(c_cm);
     free(z_cm);
 
     return info;
 }
 