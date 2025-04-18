/**
 * @file ab09md.c
 * @brief C wrapper implementation for SLICOT routine AB09MD
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB09MD,
 * which computes a reduced order model (Ar,Br,Cr) for the ALPHA-stable
 * part of an original state-space representation (A,B,C) using
 * Balance & Truncate methods.
 */

 #include <stdlib.h>
 #include <ctype.h>
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
 SLICOT_C_WRAPPER_API
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
     int ldwork = -1; /* Use -1 for workspace query */
     double dwork_query;
     double* dwork = NULL;
     int* iwork = NULL;
     int iwork_size = 0;
     const int dico_len = 1, job_len = 1, equil_len = 1, ordsel_len = 1;
 
     char dico_upper = toupper(dico);
     char job_upper = toupper(job);
     char equil_upper = toupper(equil);
     char ordsel_upper = toupper(ordsel);
 
     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL;
 
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
 
 
     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, n);
     int min_ldb_f = MAX(1, n);
     int min_ldc_f = MAX(1, p);
 
     if (row_major) {
         // For row-major C, LDA is the number of columns
         int min_lda_rm_cols = n;
         int min_ldb_rm_cols = m;
         int min_ldc_rm_cols = n;
         if (lda < min_lda_rm_cols) { info = -11; goto cleanup; }
         if (ldb < min_ldb_rm_cols) { info = -13; goto cleanup; }
         if (ldc < min_ldc_rm_cols) { info = -15; goto cleanup; }
     } else {
         // For column-major C, LDA is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -11; goto cleanup; }
         if (ldb < min_ldb_f) { info = -13; goto cleanup; }
         if (ldc < min_ldc_f) { info = -15; goto cleanup; }
     }
 
     /* --- Workspace Allocation --- */
 
     // Allocate IWORK (size N if JOB='N', 0 otherwise)
     if (job_upper == 'N') {
         iwork_size = MAX(1, n);
         iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
         CHECK_ALLOC(iwork);
     } else {
         iwork = NULL;
     }
 
     // Allocate DWORK based on query
     ldwork = -1; // Query mode
     F77_FUNC(ab09md, AB09MD)(&dico_upper, &job_upper, &equil_upper, &ordsel_upper,
                              &n, &m, &p, nr, &alpha, a, &lda, b, &ldb, c, &ldc,
                              ns, hsv, &tol, iwork, &dwork_query, &ldwork,
                              iwarn, &info,
                              dico_len, job_len, equil_len, ordsel_len);
 
     if (info != 0) { goto cleanup; } // Query failed
 
     // Get the required dwork size from query result
     ldwork = (int)dwork_query;
     // Check against minimum documented size: MAX(1,N*(2*N+MAX(N,M,P)+5)+N*(N+1)/2)
     int min_ldwork = 1;
     int max_nmp = MAX(n, MAX(m,p));
     min_ldwork = MAX(min_ldwork, n * (2 * n + max_nmp + 5) + n * (n + 1) / 2);
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
         size_t c_rows = p; size_t c_cols = n; size_t c_size = c_rows * c_cols;
 
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
 
         /* Transpose C (row-major) inputs to Fortran (column-major) copies */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         if (b_size > 0) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, elem_size);
         if (c_size > 0) slicot_transpose_to_fortran(c, c_cm, c_rows, c_cols, elem_size);
 
         /* Fortran leading dimensions */
         int lda_f = (a_rows > 0) ? a_rows : 1;
         int ldb_f = (b_rows > 0) ? b_rows : 1;
         int ldc_f = (c_rows > 0) ? c_rows : 1;
 
         /* Call the Fortran routine */
         F77_FUNC(ab09md, AB09MD)(&dico_upper, &job_upper, &equil_upper, &ordsel_upper,
                                  &n, &m, &p, nr, &alpha,
                                  a_cm, &lda_f,           // Pass CM A
                                  b_cm, &ldb_f,           // Pass CM B
                                  c_cm, &ldc_f,           // Pass CM C
                                  ns, hsv, &tol, iwork, dwork, &ldwork,
                                  iwarn, &info,
                                  dico_len, job_len, equil_len, ordsel_len);
 
         /* Copy back results from column-major temps to original row-major arrays */
         if (info == 0) {
             int nr_val = *nr; // Get the final reduced order
             // Copy back only the reduced portions
             if (nr_val > 0) {
                 if (a_size > 0) slicot_transpose_to_c(a_cm, a, nr_val, nr_val, elem_size);
                 if (b_size > 0) slicot_transpose_to_c(b_cm, b, nr_val, m, elem_size);
                 if (c_size > 0) slicot_transpose_to_c(c_cm, c, p, nr_val, elem_size);
             }
             // NS and HSV were filled directly.
         }
         /* Column-major copies will be freed in cleanup */
 
     } else {
         /* --- Column-Major Case --- */
 
         /* Call the Fortran routine directly with user-provided arrays */
         F77_FUNC(ab09md, AB09MD)(&dico_upper, &job_upper, &equil_upper, &ordsel_upper,
                                  &n, &m, &p, nr, &alpha,
                                  a, &lda,                // Pass original A
                                  b, &ldb,                // Pass original B
                                  c, &ldc,                // Pass original C
                                  ns, hsv, &tol, iwork, dwork, &ldwork,
                                  iwarn, &info,
                                  dico_len, job_len, equil_len, ordsel_len);
         // A, B, C, NS, HSV are modified in place by the Fortran call.
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
 