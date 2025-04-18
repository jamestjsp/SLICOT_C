/**
 * @file ab01md.c
 * @brief C wrapper implementation for SLICOT routine AB01MD
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB01MD,
 * which finds a controllable realization for a linear time-invariant
 * single-input system using orthogonal transformations.
 */

 #include <stdlib.h>
 #include <string.h> // For memset
 #include <ctype.h>  // For toupper
 
 // Include the header file for this wrapper
 #include "ab01md.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions
 
 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * This handles potential name mangling issues between C and Fortran compilers.
  * All arguments are passed by reference (pointers).
  * The hidden string length argument for 'jobz' is explicitly included at the end.
  */
 extern void F77_FUNC(ab01md, AB01MD)(
     const char* jobz,       // CHARACTER*1 JOBZ
     const int* n,           // INTEGER N
     double* a,              // DOUBLE PRECISION A(LDA,*) (in/out)
     const int* lda,         // INTEGER LDA
     double* b,              // DOUBLE PRECISION B(*) (in/out)
     int* ncont,             // INTEGER NCONT (output)
     double* z,              // DOUBLE PRECISION Z(LDZ,*) (output)
     const int* ldz,         // INTEGER LDZ
     double* tau,            // DOUBLE PRECISION TAU(*) (output)
     const double* tol,      // DOUBLE PRECISION TOL
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* info,              // INTEGER INFO (output)
     int jobz_len            // Hidden length argument for jobz (integer)
 );
 
 
 /* C wrapper function definition */
 int slicot_ab01md(char jobz, int n,
                   double* a, int lda,
                   double* b,
                   int* ncont,
                   double* z, int ldz,
                   double* tau, double tol, int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = -1; /* Use -1 for workspace query */
     double dwork_query;
     double* dwork = NULL;
     int jobz_len = 1; // Fortran expects 1-based length for strings
 
     char jobz_upper = toupper(jobz);
 
     /* Pointers for column-major copies if needed */
     double *a_cm = NULL;
     double *z_cm = NULL; // Needed if row_major and jobz != 'N'
 
     /* --- Input Parameter Validation --- */
 
     if (n < 0) { info = -2; goto cleanup; }
 
     if (jobz_upper != 'N' && jobz_upper != 'F' && jobz_upper != 'I') {
         info = -1; goto cleanup;
     }
 
     // Check leading dimensions based on storage order and JOBZ
     int min_lda_f = MAX(1, n);
     int min_ldz_f = (jobz_upper == 'N') ? 1 : MAX(1, n);
 
     if (row_major) {
         // For row-major C, LDA/LDZ is the number of columns
         int min_lda_rm_cols = n;
         int min_ldz_rm_cols = (jobz_upper == 'N') ? 1 : n; // Need n columns for Z
         if (lda < min_lda_rm_cols) { info = -4; goto cleanup; }
         if (jobz_upper != 'N' && ldz < min_ldz_rm_cols) { info = -8; goto cleanup; }
     } else {
         // For column-major C, LDA/LDZ is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -4; goto cleanup; }
         if (ldz < min_ldz_f) { info = -8; goto cleanup; } // Check even if JOBZ='N' as per docs
     }
 
     /* --- Workspace Query --- */
 
     // Allocate DWORK based on query
     ldwork = -1; // Query mode
     // Note: Passing actual arrays (or NULL if known safe) to query
     // Here we pass them as the routine might need N based on LDA/LDZ.
     F77_FUNC(ab01md, AB01MD)(&jobz_upper, &n, a, &lda, b, ncont, z, &ldz, tau,
                              &tol, &dwork_query, &ldwork, &info,
                              jobz_len);
 
     if (info != 0) {
         // Query failed, likely due to invalid N, LDA, LDZ etc. passed to query
         goto cleanup;
     }
 
     // Get the required dwork size from query result
     ldwork = (int)dwork_query;
     // Check against minimum documented size: MAX(1,N)
     ldwork = MAX(ldwork, MAX(1, n));
 
     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure
 
     /* --- Prepare Arrays and Call Fortran Routine --- */
 
     if (row_major) {
         /* --- Row-Major Case --- */
 
         /* Allocate memory for column-major copy of A */
         size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
         if (a_size > 0) { a_cm = (double*)malloc(a_size * sizeof(double)); CHECK_ALLOC(a_cm); }
 
         /* Allocate memory for column-major copy of Z if needed */
         size_t z_rows = n; size_t z_cols = n; size_t z_size = z_rows * z_cols;
         int ldz_f = 1; // Default for JOBZ='N'
         if (jobz_upper != 'N') {
             ldz_f = (z_rows > 0) ? z_rows : 1; // Fortran LDZ is number of rows
             if (z_size > 0) { z_cm = (double*)malloc(z_size * sizeof(double)); CHECK_ALLOC(z_cm); }
             // Initialize Z to identity if requested
             if (jobz_upper == 'I') {
                 // Create identity in z_cm (column-major)
                 set_identity(n, z_cm, ldz_f, 0); // 0 for column-major
             }
         }
 
         /* Transpose C (row-major) input A to Fortran (column-major) copy */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, sizeof(double));
 
         /* Fortran leading dimension for A */
         int lda_f = (a_rows > 0) ? a_rows : 1;
 
         /* Call the Fortran routine */
         // Pass a_cm, original b, ncont, z_cm (or NULL if JOBZ='N'), ldz_f, original tau
         F77_FUNC(ab01md, AB01MD)(&jobz_upper, &n,
                                  a_cm, &lda_f,           // Pass CM A
                                  b,                      // Pass original B (1D)
                                  ncont,                  // Pass NCONT pointer
                                  (jobz_upper == 'N' ? NULL : z_cm), // Pass CM Z or NULL
                                  &ldz_f,                 // Pass Fortran LDZ
                                  tau,                    // Pass original TAU (1D output)
                                  &tol,                   // Pass address of tol
                                  dwork, &ldwork, &info,  // Pass workspace
                                  jobz_len);              // Pass hidden length
 
         /* Copy back results */
         if (info == 0) {
             // Transpose modified a_cm back to original row-major a
             if (a_size > 0) slicot_transpose_to_c(a_cm, a, a_rows, a_cols, sizeof(double));
 
             // Transpose computed z_cm back to original row-major z (if computed)
             if (jobz_upper != 'N' && z_size > 0) {
                 slicot_transpose_to_c(z_cm, z, z_rows, z_cols, sizeof(double));
             }
             // B and TAU were modified in place (passed directly), no copy needed.
         }
         /* Column-major copies a_cm, z_cm will be freed in cleanup */
 
     } else {
         /* --- Column-Major Case --- */
 
         /* Initialize Z to identity if requested */
         if (jobz_upper == 'I') {
             set_identity(n, z, ldz, 0); // 0 for column-major
         }
 
         /* Call the Fortran routine directly with user-provided arrays */
         F77_FUNC(ab01md, AB01MD)(&jobz_upper, &n,
                                  a, &lda,                // Pass original A
                                  b,                      // Pass original B
                                  ncont,                  // Pass NCONT pointer
                                  (jobz_upper == 'N' ? NULL : z), // Pass original Z or NULL
                                  &ldz,                   // Pass original LDZ
                                  tau,                    // Pass original TAU
                                  &tol,                   // Pass address of tol
                                  dwork, &ldwork, &info,  // Pass workspace
                                  jobz_len);              // Pass hidden length
         // A, B, Z, TAU are modified in place by the Fortran call.
     }
 
 cleanup:
     /* --- Cleanup --- */
     // Free allocated workspace memory
     free(dwork);
     // Free column-major copies if they were allocated
     free(a_cm);
     free(z_cm);
 
     /* Return the info code from the Fortran routine or SLICOT_MEMORY_ERROR */
     return info;
 }
 