/**
 * @file ab13ed.c
 * @brief C wrapper implementation for SLICOT routine AB13ED
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB13ED,
 * which estimates the distance from a real matrix A to the nearest
 * complex matrix with an eigenvalue on the imaginary axis.
 */

 #include <stdlib.h>
 #include <math.h>   // For sqrt, potentially needed if TOL check involves EPS
 #include <stddef.h> // For size_t
 
 // Include the header file for this wrapper
 #include "ab13ed.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions
 
 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * Note A is input only (const). LOW, HIGH are output.
  */
 extern void F77_FUNC(ab13ed, AB13ED)(
     const int* n,           // INTEGER N
     const double* a,        // DOUBLE PRECISION A(LDA,*)
     const int* lda,         // INTEGER LDA
     double* low,            // DOUBLE PRECISION LOW (output)
     double* high,           // DOUBLE PRECISION HIGH (output)
     const double* tol,      // DOUBLE PRECISION TOL
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* info               // INTEGER INFO (output)
 );
 
 
 /* C wrapper function definition */
 int slicot_ab13ed(int n, const double* a, int lda,
                   double* low, double* high, double tol,
                   int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = -1; /* Use -1 for workspace query */
     double dwork_query;
     double* dwork = NULL;
 
     /* Pointers for column-major copies if needed */
     double *a_cm = NULL;
 
     /* --- Input Parameter Validation --- */
 
     if (n < 0) { info = -1; goto cleanup; }
     // TOL check: Fortran routine uses sqrt(eps) if tol < sqrt(eps).
     // No explicit check needed here unless we want to enforce tol >= 0.
 
     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, n);
 
     if (row_major) {
         // For row-major C, LDA is the number of columns
         int min_lda_rm_cols = n;
         if (lda < min_lda_rm_cols) { info = -3; goto cleanup; }
     } else {
         // For column-major C, LDA is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -3; goto cleanup; }
     }
 
     /* --- Workspace Allocation --- */
 
     // Allocate DWORK based on query
     ldwork = -1; // Query mode
     F77_FUNC(ab13ed, AB13ED)(&n, a, &lda, low, high, &tol,
                              &dwork_query, &ldwork, &info);
 
     if (info < 0) { goto cleanup; } // Query failed due to invalid argument
     info = 0; // Reset info after query
 
     // Get the required dwork size from query result
     ldwork = (int)dwork_query;
     // Check against minimum documented size: MAX(1, 3*N*(N+1))
     int min_ldwork = MAX(1, 3 * n * (n + 1));
     ldwork = MAX(ldwork, min_ldwork);
 
     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure
 
     /* --- Prepare Arrays and Call Fortran Routine --- */
     size_t elem_size = sizeof(double);
 
     if (row_major) {
         /* --- Row-Major Case --- */
 
         /* Allocate memory for column-major copy of A */
         size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
 
         /* Transpose C (row-major) input A to Fortran (column-major) copy */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
 
         /* Fortran leading dimension */
         int lda_f = (a_rows > 0) ? a_rows : 1;
 
         /* Call the Fortran routine */
         F77_FUNC(ab13ed, AB13ED)(&n,
                                  a_cm, &lda_f,           // Pass CM A
                                  low, high, &tol,       // Pass output pointers, address of tol
                                  dwork, &ldwork, &info); // Pass workspace
 
         /* No copy-back needed for A as it's input only */
         /* LOW and HIGH are filled directly */
         /* Temp array a_cm will be freed in cleanup */
 
     } else {
         /* --- Column-Major Case --- */
 
         /* Call the Fortran routine directly with user-provided arrays */
         F77_FUNC(ab13ed, AB13ED)(&n,
                                  a, &lda,                // Pass original A
                                  low, high, &tol,       // Pass output pointers, address of tol
                                  dwork, &ldwork, &info); // Pass workspace
         // LOW and HIGH are filled directly.
     }
 
 cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(a_cm); // Safe even if NULL
 
     return info;
 }
 