/**
 * @file sb10yd.c
 * @brief C wrapper implementation for SLICOT routine SB10YD
 *
 * This file provides a C wrapper implementation for the SLICOT routine SB10YD,
 * which fits frequency response data with a stable, minimum phase
 * SISO system.
 */

 #include <stdlib.h>
 #include <stddef.h> // For size_t
 #include <complex.h> // For C99 complex types if used by slicot_utils.h
 
 // Include the header file for this wrapper
 #include "sb10yd.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, MIN, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines, slicot_complex_double, SLICOT_COMPLEX_REAL
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions
 
 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * N is input/output. A, B, C, D are output.
  * Uses COMPLEX*16 workspace ZWORK.
  */
 extern void F77_FUNC(sb10yd, SB10YD)(
     const int* discfl,      // INTEGER DISCFL
     const int* flag,        // INTEGER FLAG
     const int* lendat,      // INTEGER LENDAT
     const double* rfrdat,   // DOUBLE PRECISION RFRDAT(*)
     const double* ifrdat,   // DOUBLE PRECISION IFRDAT(*)
     const double* omega,    // DOUBLE PRECISION OMEGA(*)
     int* n,                 // INTEGER N (in/out)
     double* a,              // DOUBLE PRECISION A(LDA,*) (output)
     const int* lda,         // INTEGER LDA
     double* b,              // DOUBLE PRECISION B(*) (output)
     double* c,              // DOUBLE PRECISION C(*) (output)
     double* d,              // DOUBLE PRECISION D(*) (output)
     const double* tol,      // DOUBLE PRECISION TOL
     int* iwork,             // INTEGER IWORK(*)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     slicot_complex_double* zwork, // COMPLEX*16 ZWORK(*)
     const int* lzwork,      // INTEGER LZWORK
     int* info               // INTEGER INFO (output)
 );
 
 
 /* C wrapper function definition */
 int slicot_sb10yd(int discfl, int flag, int lendat,
                   const double* rfrdat, const double* ifrdat, const double* omega,
                   int* n, double* a, int lda, double* b, double* c, double* d,
                   double tol, int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = -1; /* Use -1 for workspace query */
     int lzwork = -1; /* Use -1 for workspace query */
     double dwork_query[2]; // DWORK(1)=opt LDWORK, DWORK(2)=opt LZWORK
     double* dwork = NULL;
     slicot_complex_double* zwork = NULL;
     int* iwork = NULL;
     int iwork_size = 0;
     int n_in = (n != NULL) ? *n : -1; // Get initial N for validation/workspace
 
     /* Pointers for column-major copies if needed */
     double *a_cm = NULL;
 
     /* --- Input Parameter Validation --- */
 
     if (lendat < 2) { info = -3; goto cleanup; }
     if (n == NULL || a == NULL || b == NULL || c == NULL || d == NULL) {
          info = -99; goto cleanup; // Check essential output pointers
     }
     if (n_in < 0 || n_in > lendat - 1) { info = -7; goto cleanup; }
     if (discfl != 0 && discfl != 1) { info = -1; goto cleanup; }
     if (flag != 0 && flag != 1) { info = -2; goto cleanup; }
 
     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, n_in);
 
     if (row_major) {
         // For row-major C, LDA is the number of columns
         int min_lda_rm_cols = n_in;
         if (lda < min_lda_rm_cols) { info = -9; goto cleanup; }
     } else {
         // For column-major C, LDA is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -9; goto cleanup; }
     }
 
     /* --- Workspace Allocation --- */
 
     // Allocate IWORK (size max(2, 2*N+1))
     iwork_size = MAX(2, 2 * n_in + 1);
     iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
     CHECK_ALLOC(iwork);
 
     // Query DWORK and ZWORK sizes
     ldwork = -1; // Query mode
     lzwork = -1; // Query mode
     F77_FUNC(sb10yd, SB10YD)(&discfl, &flag, &lendat, rfrdat, ifrdat, omega,
                              &n_in, a, &lda, b, c, d, &tol,
                              iwork, dwork_query, &ldwork,
                              NULL, &lzwork, &info); // Pass NULL for zwork in query
 
     if (info < 0) { goto cleanup; } // Query failed due to invalid argument
     info = 0; // Reset info after query
 
     // Get the required workspace sizes from query results
     ldwork = (int)dwork_query[0];
     lzwork = (int)dwork_query[1]; // Optimal LZWORK is in DWORK(2)
 
     // Check against minimum documented sizes
     int min_ldwork = 2;
     int hnpts = 2048; // As per docs
     int lw1 = 2 * lendat + 4 * hnpts;
     int lw2 = lendat + 6 * hnpts;
     int mn = MIN(2 * lendat, 2 * n_in + 1);
     int lw3 = 0;
     if (n_in > 0) lw3 = 2*lendat*(2*n_in+1) + MAX(2*lendat, 2*n_in+1) + MAX(mn + 6*n_in + 4, 2*mn + 1);
     else          lw3 = 4*lendat + 5;
     int lw4 = 0;
     if (flag == 1) lw4 = MAX(n_in*n_in + 5*n_in, 6*n_in + 1 + MIN(1, n_in));
 
     min_ldwork = MAX(min_ldwork, MAX(lw1, MAX(lw2, MAX(lw3, lw4))));
     ldwork = MAX(ldwork, min_ldwork);
 
     int min_lzwork = (n_in > 0) ? MAX(1, lendat * (2 * n_in + 3)) : MAX(1, lendat);
     lzwork = MAX(lzwork, min_lzwork);
 
 
     // Allocate workspaces
     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork);
     zwork = (slicot_complex_double*)malloc((size_t)lzwork * sizeof(slicot_complex_double));
     CHECK_ALLOC(zwork);
 
     /* --- Prepare Arrays and Call Fortran Routine --- */
     size_t elem_size = sizeof(double);
 
     // Determine sizes for potential copies
     // Use n_in for allocation, actual N might change
     size_t a_rows = n_in; size_t a_cols = n_in; size_t a_size = a_rows * a_cols;
 
     if (row_major) {
         /* --- Row-Major Case --- */
 
         /* Allocate memory for column-major copy of A */
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
 
         /* Fortran leading dimension */
         int lda_f = (a_rows > 0) ? a_rows : 1;
 
         /* Call the Fortran routine */
         F77_FUNC(sb10yd, SB10YD)(&discfl, &flag, &lendat, rfrdat, ifrdat, omega,
                                  n,                      // Pass pointer to N (in/out)
                                  a_cm, &lda_f,           // Pass CM A (out)
                                  b, c, d, &tol,          // Pass 1D outputs, tol
                                  iwork, dwork, &ldwork, zwork, &lzwork, &info);
 
         /* Copy back results from column-major temps to original row-major arrays */
         if (info == 0) {
             int n_out = *n; // Use final N for copy back dimension
             if (n_out > 0 && a_cm) {
                  slicot_transpose_to_c(a_cm, a, n_out, n_out, elem_size);
             }
             // N, B, C, D modified directly
         }
         /* Temps freed in cleanup */
 
     } else {
         /* --- Column-Major Case --- */
 
         /* Call the Fortran routine directly with user-provided arrays */
         F77_FUNC(sb10yd, SB10YD)(&discfl, &flag, &lendat, rfrdat, ifrdat, omega,
                                  n,                      // Pass pointer to N (in/out)
                                  a, &lda,                // Pass original A
                                  b, c, d, &tol,          // Pass 1D outputs, tol
                                  iwork, dwork, &ldwork, zwork, &lzwork, &info);
         // N, A, B, C, D modified in place.
     }
 
 cleanup:
     /* --- Cleanup --- */
     free(zwork);
     free(dwork);
     free(iwork);
     free(a_cm); // Safe even if NULL
 
     return info;
 }
 