/**
 * @file ab09nd.c
 * @brief C wrapper implementation for SLICOT routine AB09ND
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB09ND,
 * which computes a reduced order model (Ar,Br,Cr,Dr) for the ALPHA-stable
 * part of an original state-space representation (A,B,C,D) using
 * Singular Perturbation Approximation (SPA) methods.
 * Refactored to align with ab01nd.c structure.
 */

 #include <stdlib.h>
 #include <ctype.h>  // For toupper
 #include <stddef.h> // For size_t

 // Include the header file for this wrapper
 #include "ab09nd.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * Note A, B, C, D are input/output. NR is input/output. NS is output.
  */
 extern void F77_FUNC(ab09nd, AB09ND)(
     const char* dico,       // CHARACTER*1 DICO
     const char* job,        // CHARACTER*1 JOB
     const char* equil,      // CHARACTER*1 EQUIL
     const char* ordsel,     // CHARACTER*1 ORDSEL
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     const int* p,           // INTEGER P
     int* nr,                // INTEGER NR (in/out)
     const double* alpha,    // DOUBLE PRECISION ALPHA
     double* a,              // DOUBLE PRECISION A(LDA,*) (in/out)
     const int* lda,         // INTEGER LDA
     double* b,              // DOUBLE PRECISION B(LDB,*) (in/out)
     const int* ldb,         // INTEGER LDB
     double* c,              // DOUBLE PRECISION C(LDC,*) (in/out)
     const int* ldc,         // INTEGER LDC
     double* d,              // DOUBLE PRECISION D(LDD,*) (in/out)
     const int* ldd,         // INTEGER LDD
     int* ns,                // INTEGER NS (output)
     double* hsv,            // DOUBLE PRECISION HSV(*) (output)
     const double* tol1,     // DOUBLE PRECISION TOL1
     const double* tol2,     // DOUBLE PRECISION TOL2
     int* iwork,             // INTEGER IWORK(*)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* iwarn,             // INTEGER IWARN (output)
     int* info,              // INTEGER INFO (output)
     int dico_len,           // Hidden length
     int job_len,            // Hidden length
     int equil_len,          // Hidden length
     int ordsel_len          // Hidden length
 );


 /* C wrapper function definition */
 SLICOT_C_WRAPPER_API
 int slicot_ab09nd(char dico, char job, char equil, char ordsel,
                   int n, int m, int p, int* nr, double alpha,
                   double* a, int lda,
                   double* b, int ldb,
                   double* c, int ldc,
                   double* d, int ldd,
                   int* ns, double* hsv, double tol1, double tol2,
                   int* iwarn, int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork;
     double* dwork = NULL;
     int* iwork = NULL;
     int iwork_size;
     const int dico_len = 1, job_len = 1, equil_len = 1, ordsel_len = 1;

     char dico_upper = toupper(dico);
     char job_upper = toupper(job);
     char equil_upper = toupper(equil);
     char ordsel_upper = toupper(ordsel);

     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL;

     /* Pointers to pass to Fortran */
     double *a_ptr, *b_ptr, *c_ptr, *d_ptr;
     int lda_f, ldb_f, ldc_f, ldd_f;

     /* --- Input Parameter Validation --- */

     if (n < 0) { info = -5; goto cleanup; }
     if (m < 0) { info = -6; goto cleanup; }
     if (p < 0) { info = -7; goto cleanup; }
     if (dico_upper != 'C' && dico_upper != 'D') { info = -1; goto cleanup; }
     if (job_upper != 'B' && job_upper != 'N') { info = -2; goto cleanup; }
     if (equil_upper != 'S' && equil_upper != 'N') { info = -3; goto cleanup; }
     if (ordsel_upper != 'F' && ordsel_upper != 'A') { info = -4; goto cleanup; }
     if (ordsel_upper == 'F' && (*nr < 0 || *nr > n) ) { info = -8; goto cleanup; }
     // ALPHA check depends on DICO
     if (dico_upper == 'C' && alpha > 0.0) { info = -9; goto cleanup; }
     if (dico_upper == 'D' && (alpha < 0.0 || alpha > 1.0)) { info = -9; goto cleanup; }
     // Check TOLs
     if (tol1 < 0.0) { info = -18; goto cleanup; }
     if (tol2 < 0.0) { info = -19; goto cleanup; }
     if (tol1 > 0.0 && tol2 > 0.0 && tol2 > tol1) { info = -19; goto cleanup; }

     // Check leading dimensions based on storage order
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
         if (n > 0 && lda < min_lda_rm_cols) { info = -11; goto cleanup; }
         if (n > 0 && ldb < min_ldb_rm_cols) { info = -13; goto cleanup; }
         if (p > 0 && ldc < min_ldc_rm_cols) { info = -15; goto cleanup; }
         if (p > 0 && ldd < min_ldd_rm_cols) { info = -17; goto cleanup; }
     } else {
         // For column-major C, LDA is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -11; goto cleanup; }
         if (ldb < min_ldb_f) { info = -13; goto cleanup; }
         if (ldc < min_ldc_f) { info = -15; goto cleanup; }
         if (ldd < min_ldd_f) { info = -17; goto cleanup; }
     }

     /* --- Prepare Arrays --- */
     size_t elem_size = sizeof(double);
     
     /* Set pointers to original arrays for workspace query */
     if (row_major) {
         /* --- Need to create temporary column-major copies --- */
         size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
         size_t b_rows = n; size_t b_cols = m; size_t b_size = b_rows * b_cols;
         size_t c_rows = p; size_t c_cols = n; size_t c_size = c_rows * c_cols;
         size_t d_rows = p; size_t d_cols = m; size_t d_size = d_rows * d_cols;

         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
         if (d_size > 0) { d_cm = (double*)malloc(d_size * elem_size); CHECK_ALLOC(d_cm); }

         /* Convert to column-major */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         if (b_size > 0) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, elem_size);
         if (c_size > 0) slicot_transpose_to_fortran(c, c_cm, c_rows, c_cols, elem_size);
         if (d_size > 0) slicot_transpose_to_fortran(d, d_cm, d_rows, d_cols, elem_size);

         /* Set pointers for Fortran and leading dimensions */
         a_ptr = a_cm; b_ptr = b_cm; c_ptr = c_cm; d_ptr = d_cm;
         lda_f = MAX(1, a_rows); 
         ldb_f = MAX(1, b_rows);
         ldc_f = MAX(1, c_rows); 
         ldd_f = MAX(1, d_rows);
     } else {
         /* --- Column-Major Case: Use original arrays --- */
         a_ptr = a; b_ptr = b; c_ptr = c; d_ptr = d;
         lda_f = lda; ldb_f = ldb; ldc_f = ldc; ldd_f = ldd;
     }

     /* --- Workspace Allocation --- */
     // Set IWORK size based on JOB parameter
     iwork_size = (job_upper == 'B') ? MAX(1, 0) : MAX(1, 2 * n);
     iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
     CHECK_ALLOC(iwork);

     // Calculate LDWORK directly using the formula
     int max_nmp = MAX(n, MAX(m, p));
     ldwork = MAX(1, n * (2 * n + max_nmp + 5) + n * (n + 1) / 2);

     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork);

     /* --- Call the computational routine --- */
     F77_FUNC(ab09nd, AB09ND)(&dico_upper, &job_upper, &equil_upper, &ordsel_upper,
                              &n, &m, &p, nr, &alpha,
                              a_ptr, &lda_f,
                              b_ptr, &ldb_f,
                              c_ptr, &ldc_f,
                              d_ptr, &ldd_f,
                              ns, hsv, &tol1, &tol2, iwork, dwork, &ldwork,
                              iwarn, &info,
                              dico_len, job_len, equil_len, ordsel_len);

     /* Copy back results from column-major temps to original row-major arrays */
     if (row_major && info == 0) {
         int nr_val = *nr; // Get the final reduced order
         // Copy back only the reduced portions for A, B, C
         if (nr_val > 0) { // Check nr_val is valid before using as dimension
             // Use slicot_transpose_to_c_with_ld to properly handle leading dimensions
             slicot_transpose_to_c_with_ld(a_cm, a, nr_val, nr_val, lda_f, lda, elem_size);
             slicot_transpose_to_c_with_ld(b_cm, b, nr_val, m, ldb_f, ldb, elem_size);
             slicot_transpose_to_c_with_ld(c_cm, c, p, nr_val, ldc_f, ldc, elem_size);
         }
         // Copy back the potentially modified D matrix
         slicot_transpose_to_c_with_ld(d_cm, d, p, m, ldd_f, ldd, elem_size);
     }

 cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(iwork);
     free(a_cm);
     free(b_cm);
     free(c_cm);
     free(d_cm);

     return info;
 }
