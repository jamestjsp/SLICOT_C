/**
 * @file sb01bd.c
 * @brief C wrapper implementation for SLICOT routine SB01BD
 *
 * This file provides a C wrapper implementation for the SLICOT routine SB01BD,
 * which determines the state feedback matrix F for a given system (A,B)
 * such that the closed-loop state matrix A+B*F has specified eigenvalues.
 */

 #include <stdlib.h>
 #include <ctype.h>
 #include <stddef.h> // For size_t
 
 // Include the header file for this wrapper
 #include "sb01bd.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions
 
 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * A, WR, WI are input/output. B is input.
  * NFP, NAP, NUP, F, Z, IWARN are output.
  */
 extern void F77_FUNC(sb01bd, SB01BD)(
     const char* dico,       // CHARACTER*1 DICO
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     const int* np,          // INTEGER NP
     const double* alpha,    // DOUBLE PRECISION ALPHA
     double* a,              // DOUBLE PRECISION A(LDA,*) (in/out)
     const int* lda,         // INTEGER LDA
     const double* b,        // DOUBLE PRECISION B(LDB,*)
     const int* ldb,         // INTEGER LDB
     double* wr,             // DOUBLE PRECISION WR(*) (in/out)
     double* wi,             // DOUBLE PRECISION WI(*) (in/out)
     int* nfp,               // INTEGER NFP (output)
     int* nap,               // INTEGER NAP (output)
     int* nup,               // INTEGER NUP (output)
     double* f,              // DOUBLE PRECISION F(LDF,*) (output)
     const int* ldf,         // INTEGER LDF
     double* z,              // DOUBLE PRECISION Z(LDZ,*) (output)
     const int* ldz,         // INTEGER LDZ
     const double* tol,      // DOUBLE PRECISION TOL
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* iwarn,             // INTEGER IWARN (output)
     int* info,              // INTEGER INFO (output)
     int dico_len            // Hidden length
 );
 
 
 /* C wrapper function definition */
 SLICOT_C_WRAPPER_API
 int slicot_sb01bd(char dico, int n, int m, int np, double alpha,
                   double* a, int lda, const double* b, int ldb,
                   double* wr, double* wi,
                   int* nfp, int* nap, int* nup,
                   double* f, int ldf, double* z, int ldz,
                   double tol, int* iwarn, int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = -1; /* Use -1 for workspace query */
     double dwork_query;
     double* dwork = NULL;
     const int dico_len = 1;
 
     char dico_upper = toupper(dico);
 
     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *b_cm = NULL, *f_cm = NULL, *z_cm = NULL;
 
     /* --- Input Parameter Validation --- */
 
     if (n < 0) { info = -2; goto cleanup; }
     if (m < 0) { info = -3; goto cleanup; }
     if (np < 0) { info = -4; goto cleanup; } // NP can be 0
     if (dico_upper != 'C' && dico_upper != 'D') { info = -1; goto cleanup; }
     // ALPHA check depends on DICO
     if (dico_upper == 'D' && alpha < 0.0) { info = -5; goto cleanup; }
 
     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, n);
     int min_ldb_f = MAX(1, n);
     int min_ldf_f = MAX(1, m);
     int min_ldz_f = MAX(1, n);
 
     if (row_major) {
         // For row-major C, LDA is the number of columns
         int min_lda_rm_cols = n;
         int min_ldb_rm_cols = m;
         int min_ldf_rm_cols = n;
         int min_ldz_rm_cols = n;
         if (lda < min_lda_rm_cols) { info = -7; goto cleanup; }
         if (ldb < min_ldb_rm_cols) { info = -9; goto cleanup; }
         if (ldf < min_ldf_rm_cols) { info = -16; goto cleanup; }
         if (ldz < min_ldz_rm_cols) { info = -18; goto cleanup; }
     } else {
         // For column-major C, LDA is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -7; goto cleanup; }
         if (ldb < min_ldb_f) { info = -9; goto cleanup; }
         if (ldf < min_ldf_f) { info = -16; goto cleanup; }
         if (ldz < min_ldz_f) { info = -18; goto cleanup; }
     }
 
     /* --- Workspace Allocation --- */
 
     // Allocate DWORK based on query
     ldwork = -1; // Query mode
     F77_FUNC(sb01bd, SB01BD)(&dico_upper, &n, &m, &np, &alpha, a, &lda, b, &ldb,
                              wr, wi, nfp, nap, nup, f, &ldf, z, &ldz,
                              &tol, &dwork_query, &ldwork, iwarn, &info,
                              dico_len);
 
     if (info < 0) { goto cleanup; } // Query failed due to invalid argument
     info = 0; // Reset info after query
 
     // Get the required dwork size from query result
     ldwork = (int)dwork_query;
     // Check against minimum documented size: MAX( 1,5*M,5*N,2*N+4*M )
     int min_ldwork = 1;
     min_ldwork = MAX(min_ldwork, 5 * m);
     min_ldwork = MAX(min_ldwork, 5 * n);
     min_ldwork = MAX(min_ldwork, 2 * n + 4 * m);
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
         size_t f_rows = m; size_t f_cols = n; size_t f_size = f_rows * f_cols;
         size_t z_rows = n; size_t z_cols = n; size_t z_size = z_rows * z_cols;
 
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (f_size > 0) { f_cm = (double*)malloc(f_size * elem_size); CHECK_ALLOC(f_cm); }
         if (z_size > 0) { z_cm = (double*)malloc(z_size * elem_size); CHECK_ALLOC(z_cm); }
 
         /* Transpose C (row-major) inputs to Fortran (column-major) copies */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         if (b_size > 0) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, elem_size);
         // WR, WI are 1D, no transpose needed
 
         /* Fortran leading dimensions */
         int lda_f = (a_rows > 0) ? a_rows : 1;
         int ldb_f = (b_rows > 0) ? b_rows : 1;
         int ldf_f = (f_rows > 0) ? f_rows : 1;
         int ldz_f = (z_rows > 0) ? z_rows : 1;
 
         /* Call the Fortran routine */
         F77_FUNC(sb01bd, SB01BD)(&dico_upper, &n, &m, &np, &alpha,
                                  a_cm, &lda_f,           // Pass CM A (in/out)
                                  b_cm, &ldb_f,           // Pass CM B (in)
                                  wr, wi,                 // Pass WR, WI (in/out)
                                  nfp, nap, nup,          // Pass output pointers
                                  f_cm, &ldf_f,           // Pass CM F (out)
                                  z_cm, &ldz_f,           // Pass CM Z (out)
                                  &tol, dwork, &ldwork, iwarn, &info,
                                  dico_len);
 
         /* Copy back results from column-major temps to original row-major arrays */
         if (info == 0 || info == 3 || info == 4) { // Copy back even on INFO 3, 4
             if (a_size > 0) slicot_transpose_to_c(a_cm, a, a_rows, a_cols, elem_size); // Modified A
             if (f_size > 0) slicot_transpose_to_c(f_cm, f, f_rows, f_cols, elem_size); // Output F
             if (z_size > 0) slicot_transpose_to_c(z_cm, z, z_rows, z_cols, elem_size); // Output Z
             // WR, WI, NFP, NAP, NUP, IWARN modified directly
         }
         /* Column-major copies will be freed in cleanup */
 
     } else {
         /* --- Column-Major Case --- */
 
         /* Call the Fortran routine directly with user-provided arrays */
         F77_FUNC(sb01bd, SB01BD)(&dico_upper, &n, &m, &np, &alpha,
                                  a, &lda,                // Pass original A
                                  b, &ldb,                // Pass original B
                                  wr, wi,                 // Pass original WR, WI
                                  nfp, nap, nup,          // Pass output pointers
                                  f, &ldf,                // Pass original F
                                  z, &ldz,                // Pass original Z
                                  &tol, dwork, &ldwork, iwarn, &info,
                                  dico_len);
         // A, WR, WI, NFP, NAP, NUP, F, Z, IWARN are modified in place.
     }
 
 cleanup:
     /* --- Cleanup --- */
     free(dwork);
     // No iwork for this routine
     free(a_cm);
     free(b_cm);
     free(f_cm);
     free(z_cm);
 
     return info;
 }
 