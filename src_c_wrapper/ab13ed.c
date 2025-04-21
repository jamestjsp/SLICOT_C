/**
 * @file ab13ed.c
 * @brief C wrapper implementation for SLICOT routine AB13ED
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB13ED,
 * which estimates the distance from a real matrix A to the nearest
 * complex matrix with an eigenvalue on the imaginary axis.
 * Refactored to align with ab01nd.c structure.
 */

 #include <stdlib.h>
 #include <math.h>   // For isnan, isinf (optional checks)
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
 SLICOT_C_WRAPPER_API
 int slicot_ab13ed(int n, const double* a, int lda,
                   double* low, double* high, double tol,
                   int row_major)
 {
     /* Local variables */
     int info = 0;
     double* dwork = NULL;

     /* Pointers for column-major copies if needed */
     double *a_cm = NULL;

     /* Pointers to pass to Fortran */
     const double *a_ptr;
     int lda_f;

     /* --- Input Parameter Validation --- */

     if (n < 0) { info = -1; goto cleanup; }
     // TOL check: Fortran routine uses sqrt(eps) if tol < sqrt(eps).
     // No explicit check needed here unless we want to enforce tol >= 0.
     // if (tol < 0.0) { info = -6; goto cleanup; }

     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, n);

     if (row_major) {
         // For row-major C, LDA is the number of columns
         int min_lda_rm_cols = n;
         if (n > 0 && lda < min_lda_rm_cols) { info = -3; goto cleanup; }
     } else {
         // For column-major C, LDA is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -3; goto cleanup; }
     }

     /* --- Workspace Allocation --- */

     // Calculate the minimum required workspace size: MAX(1, 3*N*(N+1))
     int ldwork = MAX(1, 3 * n * (n + 1));

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
         if (a_cm) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);

         /* Fortran leading dimension */
         lda_f = MAX(1, a_rows);

         /* Set pointers for Fortran call */
         a_ptr = a_cm;

     } else {
         /* --- Column-Major Case --- */
         lda_f = lda;
         a_ptr = a;
     }

     /* Call the computational routine */
     F77_FUNC(ab13ed, AB13ED)(&n,
                              a_ptr, &lda_f,           // Pass A ptr
                              low, high, &tol,       // Pass output pointers, address of tol
                              dwork, &ldwork, &info); // Pass workspace

     /* No copy-back needed for A as it's input only */
     /* LOW and HIGH are filled directly */

 cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(a_cm); // Safe even if NULL

     return info;
 }
