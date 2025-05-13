/**
 * @file ab09md.c
 * @brief C wrapper implementation for SLICOT routine AB09MD
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB09MD,
 * which computes a reduced order model (Ar,Br,Cr) for the ALPHA-stable
 * part of an original state-space representation (A,B,C) using
 * Balance & Truncate methods.
 * Refactored to align with ab01nd.c structure.
 */

 #include <stdlib.h>
 #include <ctype.h>  // For toupper
 #include <stddef.h> // For size_t

 // Include the header file for this wrapper
 #include "ab09md.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * Note A, B, C are input/output. NR is input/output. NS is output.
  */
 extern void F77_FUNC(ab09md, AB09MD)(
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
     int* ns,                // INTEGER NS (output)
     double* hsv,            // DOUBLE PRECISION HSV(*) (output)
     const double* tol,      // DOUBLE PRECISION TOL
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
 SLICOT_EXPORT
 int slicot_ab09md(char dico, char job, char equil, char ordsel,
                   int n, int m, int p, int* nr, double alpha,
                   double* a, int lda,
                   double* b, int ldb,
                   double* c, int ldc,
                   int* ns, double* hsv, double tol, int* iwarn,
                   int row_major)
 {
     /* Local variables */
     int info = 0;
     int* iwork = NULL;
     double* dwork = NULL;
     int iwork_size = 0;
     int ldwork = -1; /* Use -1 for workspace query */
     const int dico_len = 1, job_len = 1, equil_len = 1, ordsel_len = 1;

     char dico_upper = toupper(dico);
     char job_upper = toupper(job);
     char equil_upper = toupper(equil);
     char ordsel_upper = toupper(ordsel);

     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL;

     /* Pointers to pass to Fortran */
     double *a_ptr, *b_ptr, *c_ptr;
     int lda_f, ldb_f, ldc_f;


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
     // Optional: Check TOL range if necessary
     // if (tol < 0.0) { info = -18; goto cleanup; }

     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, n);
     int min_ldb_f = MAX(1, n);
     int min_ldc_f = MAX(1, p);

     if (row_major) {
         // For row-major C, LDA is the number of columns
         if (n > 0 && lda < n) { info = -11; goto cleanup; }
         if (n > 0 && m > 0 && ldb < m) { info = -13; goto cleanup; }
         if (p > 0 && n > 0 && ldc < n) { info = -15; goto cleanup; }
     } else {
         // For column-major C, LDA is the number of rows (Fortran style)
         // Only check if dimensions are positive - critical for zero-dimension case
         if (n > 0 && lda < min_lda_f) { info = -11; goto cleanup; }
         if (n > 0 && ldb < min_ldb_f) { info = -13; goto cleanup; }
         if (p > 0 && n > 0 && ldc < min_ldc_f) { info = -15; goto cleanup; }
     }

     /* --- Workspace Allocation --- */

     // Allocate IWORK (size N if JOB='N', 0 otherwise)
     if (job_upper == 'N') {
         iwork_size = MAX(1, n);
         if (n > 0) {  // Only allocate if n > 0
             iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
             CHECK_ALLOC(iwork);
         } else {
             iwork = NULL;  // For n=0, pass NULL to Fortran
         }
     } else {
         iwork = NULL; // Pass NULL for size 0
     }

     // Call Fortran routine for DWORK workspace query
     ldwork = -1; // Query mode
     double dwork_query[1]; /* Use array for workspace query */
     
     // Use dummy LDs for query
     int lda_q = row_major ? MAX(1, n) : lda;
     int ldb_q = row_major ? MAX(1, n) : ldb;
     int ldc_q = row_major ? MAX(1, p) : ldc;
     int nr_q = *nr;  // Save original value
     int ns_q = 0;    // Dummy output variable
     int iwarn_q = 0; // Dummy output variable
     
     F77_FUNC(ab09md, AB09MD)(&dico_upper, &job_upper, &equil_upper, &ordsel_upper,
                              &n, &m, &p, &nr_q, &alpha,
                              NULL, &lda_q, NULL, &ldb_q, NULL, &ldc_q,
                              &ns_q, NULL, &tol, iwork, dwork_query, &ldwork,
                              &iwarn_q, &info,
                              dico_len, job_len, equil_len, ordsel_len);

     if (info != 0) {
         // Query failed
         goto cleanup;
     }
     
     // Get the required dwork size from query result
     ldwork = (int)dwork_query[0];
     
     // Ensure minimum size based on documentation formula
     int min_ldwork_doc = MAX(1, n*(2*n+MAX(n,MAX(m,p))+5) + n*(n+1)/2);
     ldwork = MAX(ldwork, min_ldwork_doc);

     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork);

     /* --- Prepare Arrays and Call Fortran Routine --- */
     size_t elem_size = sizeof(double);

     if (row_major) {
         /* --- Row-Major Case --- */

         /* Allocate memory for column-major copies */
         size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
         size_t b_rows = n; size_t b_cols = m; size_t b_size = b_rows * b_cols;
         size_t c_rows = p; size_t c_cols = n; size_t c_size = c_rows * c_cols;

         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }

         /* Transpose C (row-major) inputs to Fortran (column-major) copies */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         if (b_size > 0) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, elem_size);
         if (c_size > 0) slicot_transpose_to_fortran(c, c_cm, c_rows, c_cols, elem_size);

         /* Fortran leading dimensions */
         lda_f = MAX(1, a_rows);
         ldb_f = MAX(1, b_rows);
         ldc_f = MAX(1, c_rows);

         /* Set pointers for Fortran call */
         a_ptr = a_size > 0 ? a_cm : NULL;
         b_ptr = b_size > 0 ? b_cm : NULL;
         c_ptr = c_size > 0 ? c_cm : NULL;

     } else {
         /* --- Column-Major Case --- */
         lda_f = lda; ldb_f = ldb; ldc_f = ldc;
         
         // Set pointers - pass NULL if dimensions are zero
         size_t a_size = (size_t)n * n;
         size_t b_size = (size_t)n * m;
         size_t c_size = (size_t)p * n;
         
         a_ptr = (a_size > 0) ? a : NULL;
         b_ptr = (b_size > 0) ? b : NULL;
         c_ptr = (c_size > 0) ? c : NULL;
     }

     /* Call the computational routine */
     F77_FUNC(ab09md, AB09MD)(&dico_upper, &job_upper, &equil_upper, &ordsel_upper,
                              &n, &m, &p, nr, &alpha,
                              a_ptr, &lda_f,           // Pass A ptr
                              b_ptr, &ldb_f,           // Pass B ptr
                              c_ptr, &ldc_f,           // Pass C ptr
                              ns, hsv, &tol, iwork, dwork, &ldwork,
                              iwarn, &info,
                              dico_len, job_len, equil_len, ordsel_len);

     /* Copy back results from column-major temps to original row-major arrays */
     if (row_major && info == 0) {
         int nr_val = *nr; // Get the final reduced order
         // Copy back only the reduced portions
         if (nr_val > 0) { // Only copy if we have a positive reduced order
             size_t ar_rows = nr_val; size_t ar_cols = nr_val;
             size_t br_rows = nr_val; size_t br_cols = m;
             size_t cr_rows = p; size_t cr_cols = nr_val;

             // Make sure the destination array has enough space
             if (ar_rows * ar_cols > 0) { 
                 if (lda < nr_val) {
                     info = -11; // Insufficient columns in row-major A
                     goto cleanup;
                 }
                 slicot_transpose_to_c_with_ld(a_cm, a, ar_rows, ar_cols, ar_rows, lda, elem_size);
             }
             
             if (br_rows * br_cols > 0) { 
                 if (ldb < m) {
                     info = -13; // Insufficient columns in row-major B
                     goto cleanup;
                 }
                 slicot_transpose_to_c_with_ld(b_cm, b, br_rows, br_cols, br_rows, ldb, elem_size);
             }
             
             if (cr_rows * cr_cols > 0) { 
                 if (ldc < nr_val) {
                     info = -15; // Insufficient columns in row-major C
                     goto cleanup;
                 }
                 slicot_transpose_to_c_with_ld(c_cm, c, cr_rows, cr_cols, cr_rows, ldc, elem_size);
             }
         }
         // NS and HSV were filled directly.
     }

 cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(iwork); // Safe even if NULL
     free(a_cm);
     free(b_cm);
     free(c_cm);

     return info;
 }
