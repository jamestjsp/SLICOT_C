/**
 * @file sb10jd.c
 * @brief C wrapper implementation for SLICOT routine SB10JD
 *
 * This file provides a C wrapper implementation for the SLICOT routine SB10JD,
 * which converts a descriptor state-space system into regular
 * state-space form.
 */

 #include <stdlib.h>
 #include <stddef.h> // For size_t
 
 // Include the header file for this wrapper
 #include "sb10jd.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions
 
 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * A, B, C, D, E are input/output. NSYS is output.
  */
 extern void F77_FUNC(sb10jd, SB10JD)(
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     const int* np,          // INTEGER NP
     double* a,              // DOUBLE PRECISION A(LDA,*) (in/out -> Ad)
     const int* lda,         // INTEGER LDA
     double* b,              // DOUBLE PRECISION B(LDB,*) (in/out -> Bd)
     const int* ldb,         // INTEGER LDB
     double* c,              // DOUBLE PRECISION C(LDC,*) (in/out -> Cd)
     const int* ldc,         // INTEGER LDC
     double* d,              // DOUBLE PRECISION D(LDD,*) (in/out -> Dd)
     const int* ldd,         // INTEGER LDD
     double* e,              // DOUBLE PRECISION E(LDE,*) (in/out -> destroyed)
     const int* lde,         // INTEGER LDE
     int* nsys,              // INTEGER NSYS (output)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* info               // INTEGER INFO (output)
 );
 
 
 /* C wrapper function definition */
 int slicot_sb10jd(int n, int m, int np,
                   double* a, int lda, double* b, int ldb,
                   double* c, int ldc, double* d, int ldd,
                   double* e, int lde, int* nsys,
                   int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = -1; /* Use -1 for workspace query */
     double dwork_query;
     double* dwork = NULL;
 
     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL, *e_cm = NULL;
 
     /* --- Input Parameter Validation --- */
 
     if (n < 0) { info = -1; goto cleanup; }
     if (m < 0) { info = -2; goto cleanup; }
     if (np < 0) { info = -3; goto cleanup; }
 
     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, n);
     int min_ldb_f = MAX(1, n);
     int min_ldc_f = MAX(1, np);
     int min_ldd_f = MAX(1, np);
     int min_lde_f = MAX(1, n);
 
     if (row_major) {
         // For row-major C, LD is the number of columns
         int min_lda_rm_cols = n;
         int min_ldb_rm_cols = m;
         int min_ldc_rm_cols = n;
         int min_ldd_rm_cols = m;
         int min_lde_rm_cols = n;
         if (lda < min_lda_rm_cols) { info = -5; goto cleanup; }
         if (ldb < min_ldb_rm_cols) { info = -7; goto cleanup; }
         if (ldc < min_ldc_rm_cols) { info = -9; goto cleanup; }
         if (ldd < min_ldd_rm_cols) { info = -11; goto cleanup; }
         if (lde < min_lde_rm_cols) { info = -13; goto cleanup; }
     } else {
         // For column-major C, LD is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -5; goto cleanup; }
         if (ldb < min_ldb_f) { info = -7; goto cleanup; }
         if (ldc < min_ldc_f) { info = -9; goto cleanup; }
         if (ldd < min_ldd_f) { info = -11; goto cleanup; }
         if (lde < min_lde_f) { info = -13; goto cleanup; }
     }
 
     /* --- Workspace Allocation --- */
 
     // Allocate DWORK based on query
     ldwork = -1; // Query mode
     F77_FUNC(sb10jd, SB10JD)(&n, &m, &np, a, &lda, b, &ldb, c, &ldc, d, &ldd,
                              e, &lde, nsys, &dwork_query, &ldwork, &info);
 
     if (info < 0 && info != -16) { info = info; goto cleanup; } // Query failed due to invalid argument
     info = 0; // Reset info after query
 
     // Get the required dwork size from query result
     ldwork = (int)dwork_query;
     // Check against minimum documented size: max( 1, 2*N*N + 2*N + N*MAX( 5, N + M + NP ) )
     int min_ldwork = 1;
     min_ldwork = MAX(min_ldwork, 2*n*n + 2*n + n*MAX(5, n + m + np));
     ldwork = MAX(ldwork, min_ldwork);
 
     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure
 
     /* --- Prepare Arrays and Call Fortran Routine --- */
     size_t elem_size = sizeof(double);
 
     // Determine sizes for potential copies
     size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
     size_t b_rows = n; size_t b_cols = m; size_t b_size = b_rows * b_cols;
     size_t c_rows = np; size_t c_cols = n; size_t c_size = c_rows * c_cols;
     size_t d_rows = np; size_t d_cols = m; size_t d_size = d_rows * d_cols;
     size_t e_rows = n; size_t e_cols = n; size_t e_size = e_rows * e_cols;
 
 
     if (row_major) {
         /* --- Row-Major Case --- */
 
         /* Allocate memory for column-major copies */
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
         if (d_size > 0) { d_cm = (double*)malloc(d_size * elem_size); CHECK_ALLOC(d_cm); }
         if (e_size > 0) { e_cm = (double*)malloc(e_size * elem_size); CHECK_ALLOC(e_cm); }
 
         /* Transpose C inputs to Fortran copies */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         if (b_size > 0) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, elem_size);
         if (c_size > 0) slicot_transpose_to_fortran(c, c_cm, c_rows, c_cols, elem_size);
         if (d_size > 0) slicot_transpose_to_fortran(d, d_cm, d_rows, d_cols, elem_size);
         if (e_size > 0) slicot_transpose_to_fortran(e, e_cm, e_rows, e_cols, elem_size);
 
         /* Fortran leading dimensions */
         int lda_f = (a_rows > 0) ? a_rows : 1;
         int ldb_f = (b_rows > 0) ? b_rows : 1;
         int ldc_f = (c_rows > 0) ? c_rows : 1;
         int ldd_f = (d_rows > 0) ? d_rows : 1;
         int lde_f = (e_rows > 0) ? e_rows : 1;
 
         /* Call the Fortran routine */
         F77_FUNC(sb10jd, SB10JD)(&n, &m, &np,
                                  a_cm, &lda_f, b_cm, &ldb_f, c_cm, &ldc_f, d_cm, &ldd_f,
                                  e_cm, &lde_f, nsys, dwork, &ldwork, &info);
 
         /* Copy back results from column-major temps to original row-major arrays */
         if (info == 0) {
             int nsys_val = *nsys; // Get the order of the converted system
             if (a_size > 0) slicot_transpose_to_c(a_cm, a, nsys_val, nsys_val, elem_size); // Ad
             if (b_size > 0) slicot_transpose_to_c(b_cm, b, nsys_val, m, elem_size);       // Bd
             if (c_size > 0) slicot_transpose_to_c(c_cm, c, np, nsys_val, elem_size);       // Cd
             if (d_size > 0) slicot_transpose_to_c(d_cm, d, np, m, elem_size);       // Dd
             // NSYS modified directly
         }
         /* Temps freed in cleanup */
 
     } else {
         /* --- Column-Major Case --- */
 
         /* Call the Fortran routine directly with user-provided arrays */
         F77_FUNC(sb10jd, SB10JD)(&n, &m, &np, a, &lda, b, &ldb, c, &ldc, d, &ldd,
                                  e, &lde, nsys, dwork, &ldwork, &info);
         // A, B, C, D, E, NSYS modified in place.
     }
 
 cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(a_cm);
     free(b_cm);
     free(c_cm);
     free(d_cm);
     free(e_cm);
 
     return info;
 }
 