/**
 * @file ab01nd.c
 * @brief C wrapper implementation for SLICOT routine AB01ND
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB01ND,
 * which finds a controllable realization for a linear time-invariant
 * multi-input system using orthogonal transformations.
 */

 #include <stdlib.h>
 #include <string.h> // For memset
 #include <ctype.h>  // For toupper
 
 // Include the header file for this wrapper
 #include "ab01nd.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines, set_identity
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions
 
 /* Helper function to set a matrix to identity (declaration assumed in slicot_utils.h) */
 /* If not declared there, uncomment or add declaration here: */
 /* static void set_identity(int n, double* mat, int ld, int row_major); */
 
 
 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * This handles potential name mangling issues between C and Fortran compilers.
  * All arguments are passed by reference (pointers).
  * The hidden string length argument for 'jobz' is explicitly included at the end.
  */
 extern void F77_FUNC(ab01nd, AB01ND)(
     const char* jobz,       // CHARACTER*1 JOBZ
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     double* a,              // DOUBLE PRECISION A(LDA,*) (in/out)
     const int* lda,         // INTEGER LDA
     double* b,              // DOUBLE PRECISION B(LDB,*) (in/out)
     const int* ldb,         // INTEGER LDB
     int* ncont,             // INTEGER NCONT (output)
     int* indcon,            // INTEGER INDCON (output)
     int* nblk,              // INTEGER NBLK(*) (output)
     double* z,              // DOUBLE PRECISION Z(LDZ,*) (output)
     const int* ldz,         // INTEGER LDZ
     double* tau,            // DOUBLE PRECISION TAU(*) (output)
     const double* tol,      // DOUBLE PRECISION TOL
     int* iwork,             // INTEGER IWORK(*)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* info,              // INTEGER INFO (output)
     int jobz_len            // Hidden length argument for jobz (integer)
 );
 
 
 /* C wrapper function definition */
 int slicot_ab01nd(char jobz, int n, int m,
                   double* a, int lda,
                   double* b, int ldb,
                   int* ncont, int* indcon, int* nblk,
                   double* z, int ldz,
                   double* tau, double tol,
                   int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = -1; /* Use -1 for workspace query */
     double dwork_query;
     double* dwork = NULL;
     int* iwork = NULL;
     int iwork_size = 0;
     int jobz_len = 1; // Fortran expects 1-based length for strings
 
     char jobz_upper = toupper(jobz);
 
     /* Pointers for column-major copies if needed */
     double *a_cm = NULL;
     double *b_cm = NULL;
     double *z_cm = NULL; // Needed if row_major and jobz != 'N'
 
     /* --- Input Parameter Validation --- */
 
     if (n < 0) { info = -2; goto cleanup; }
     if (m < 0) { info = -3; goto cleanup; }
 
     if (jobz_upper != 'N' && jobz_upper != 'F' && jobz_upper != 'I') {
         info = -1; goto cleanup;
     }
 
     // Check leading dimensions based on storage order and JOBZ
     int min_lda_f = MAX(1, n);
     int min_ldb_f = MAX(1, n);
     int min_ldz_f = (jobz_upper == 'N') ? 1 : MAX(1, n);
 
     if (row_major) {
         // For row-major C, LDA/LDB/LDZ is the number of columns
         int min_lda_rm_cols = n;
         int min_ldb_rm_cols = m;
         int min_ldz_rm_cols = (jobz_upper == 'N') ? 1 : n; // Need n columns for Z
         if (lda < min_lda_rm_cols) { info = -5; goto cleanup; }
         if (ldb < min_ldb_rm_cols) { info = -7; goto cleanup; }
         if (jobz_upper != 'N' && ldz < min_ldz_rm_cols) { info = -12; goto cleanup; }
     } else {
         // For column-major C, LDA/LDB/LDZ is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -5; goto cleanup; }
         if (ldb < min_ldb_f) { info = -7; goto cleanup; }
         if (ldz < min_ldz_f) { info = -12; goto cleanup; } // Check even if JOBZ='N' as per docs
     }
 
     /* --- Workspace Allocation --- */
 
     // Allocate IWORK (size M)
     iwork_size = MAX(1, m); // Ensure minimum size 1
     iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
     CHECK_ALLOC(iwork);
 
     // Allocate DWORK based on query
     ldwork = -1; // Query mode
     // Note: Passing actual arrays (or NULL if known safe) to query
     F77_FUNC(ab01nd, AB01ND)(&jobz_upper, &n, &m, a, &lda, b, &ldb,
                              ncont, indcon, nblk, z, &ldz, tau, &tol,
                              iwork, &dwork_query, &ldwork, &info,
                              jobz_len);
 
     if (info != 0) {
         // Query failed, likely due to invalid N, M, LDA, LDB, LDZ etc. passed to query
         goto cleanup;
     }
 
     // Get the required dwork size from query result
     ldwork = (int)dwork_query;
     // Check against minimum documented size: MAX(1, N, 3*M)
     ldwork = MAX(ldwork, MAX(1, MAX(n, 3 * m)));
 
     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure
 
     /* --- Prepare Arrays and Call Fortran Routine --- */
 
     if (row_major) {
         /* --- Row-Major Case --- */
 
         /* Allocate memory for column-major copies */
         size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
         size_t b_rows = n; size_t b_cols = m; size_t b_size = b_rows * b_cols;
         size_t z_rows = n; size_t z_cols = n; size_t z_size = z_rows * z_cols;
         int ldz_f = 1; // Default for JOBZ='N'
 
         if (a_size > 0) { a_cm = (double*)malloc(a_size * sizeof(double)); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * sizeof(double)); CHECK_ALLOC(b_cm); }
 
         if (jobz_upper != 'N') {
             ldz_f = (z_rows > 0) ? z_rows : 1; // Fortran LDZ is number of rows
             if (z_size > 0) { z_cm = (double*)malloc(z_size * sizeof(double)); CHECK_ALLOC(z_cm); }
             // Initialize Z to identity if requested
             if (jobz_upper == 'I') {
                 // Create identity in z_cm (column-major)
                 set_identity(n, z_cm, ldz_f, 0); // 0 for column-major
             }
         }
 
         /* Transpose C (row-major) inputs to Fortran (column-major) copies */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, sizeof(double));
         if (b_size > 0) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, sizeof(double));
 
         /* Fortran leading dimensions */
         int lda_f = (a_rows > 0) ? a_rows : 1;
         int ldb_f = (b_rows > 0) ? b_rows : 1;
         // ldz_f calculated above
 
         /* Call the Fortran routine */
         F77_FUNC(ab01nd, AB01ND)(&jobz_upper, &n, &m,
                                  a_cm, &lda_f,           // Pass CM A
                                  b_cm, &ldb_f,           // Pass CM B
                                  ncont, indcon, nblk,    // Pass output pointers
                                  (jobz_upper == 'N' ? NULL : z_cm), // Pass CM Z or NULL
                                  &ldz_f,                 // Pass Fortran LDZ
                                  tau,                    // Pass original TAU (1D output)
                                  &tol,                   // Pass address of tol
                                  iwork, dwork, &ldwork, &info, // Pass workspaces
                                  jobz_len);              // Pass hidden length
 
         /* Copy back results */
         if (info == 0) {
             // Transpose modified a_cm back to original row-major a
             if (a_size > 0) slicot_transpose_to_c(a_cm, a, a_rows, a_cols, sizeof(double));
             // Transpose modified b_cm back to original row-major b
             if (b_size > 0) slicot_transpose_to_c(b_cm, b, b_rows, b_cols, sizeof(double));
             // Transpose computed z_cm back to original row-major z (if computed)
             if (jobz_upper != 'N' && z_size > 0) {
                 slicot_transpose_to_c(z_cm, z, z_rows, z_cols, sizeof(double));
             }
             // NBLK and TAU were modified in place (passed directly), no copy needed.
         }
         /* Column-major copies a_cm, b_cm, z_cm will be freed in cleanup */
 
     } else {
         /* --- Column-Major Case --- */
 
         /* Initialize Z to identity if requested */
         if (jobz_upper == 'I') {
             set_identity(n, z, ldz, 0); // 0 for column-major
         }
 
         /* Call the Fortran routine directly with user-provided arrays */
         F77_FUNC(ab01nd, AB01ND)(&jobz_upper, &n, &m,
                                  a, &lda,                // Pass original A
                                  b, &ldb,                // Pass original B
                                  ncont, indcon, nblk,    // Pass output pointers
                                  (jobz_upper == 'N' ? NULL : z), // Pass original Z or NULL
                                  &ldz,                   // Pass original LDZ
                                  tau,                    // Pass original TAU
                                  &tol,                   // Pass address of tol
                                  iwork, dwork, &ldwork, &info, // Pass workspaces
                                  jobz_len);              // Pass hidden length
         // A, B, NBLK, Z, TAU are modified in place by the Fortran call.
     }
 
 cleanup:
     /* --- Cleanup --- */
     // Free allocated workspace memory
     free(dwork);
     free(iwork);
     // Free column-major copies if they were allocated
     free(a_cm);
     free(b_cm);
     free(z_cm);
 
     /* Return the info code from the Fortran routine or SLICOT_MEMORY_ERROR */
     return info;
 }
 