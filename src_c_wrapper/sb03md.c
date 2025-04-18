/**
 * @file sb03md.c
 * @brief C wrapper implementation for SLICOT routine SB03MD
 *
 * This file provides a C wrapper implementation for the SLICOT routine SB03MD,
 * which solves continuous- or discrete-time Lyapunov equations and
 * optionally estimates the separation condition number.
 */

 #include <stdlib.h>
 #include <ctype.h>
 #include <stddef.h> // For size_t
 #include <string.h> // For memcpy
 
 // Include the header file for this wrapper
 #include "sb03md.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines, slicot_copy_symmetric_part
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions
 
 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * A, U, C are input/output depending on flags.
  * SCALE, SEP, FERR, WR, WI are output.
  */
 extern void F77_FUNC(sb03md, SB03MD)(
     const char* dico,       // CHARACTER*1 DICO
     const char* job,        // CHARACTER*1 JOB
     const char* fact,       // CHARACTER*1 FACT
     const char* trana,      // CHARACTER*1 TRANA
     const int* n,           // INTEGER N
     double* a,              // DOUBLE PRECISION A(LDA,*) (in/out)
     const int* lda,         // INTEGER LDA
     double* u,              // DOUBLE PRECISION U(LDU,*) (in/out)
     const int* ldu,         // INTEGER LDU
     double* c,              // DOUBLE PRECISION C(LDC,*) (in/out)
     const int* ldc,         // INTEGER LDC
     double* scale,          // DOUBLE PRECISION SCALE (output)
     double* sep,            // DOUBLE PRECISION SEP (output)
     double* ferr,           // DOUBLE PRECISION FERR (output)
     double* wr,             // DOUBLE PRECISION WR(*) (output)
     double* wi,             // DOUBLE PRECISION WI(*) (output)
     int* iwork,             // INTEGER IWORK(*)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* info,              // INTEGER INFO (output)
     int dico_len,           // Hidden length
     int job_len,            // Hidden length
     int fact_len,           // Hidden length
     int trana_len           // Hidden length
 );
 
 
 /* C wrapper function definition */
 SLICOT_C_WRAPPER_API
 int slicot_sb03md(char dico, char job, char fact, char trana, int n,
                   double* a, int lda, double* u, int ldu,
                   double* c, int ldc, double* scale, double* sep, double* ferr,
                   double* wr, double* wi, int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = -1; /* Use -1 for workspace query */
     double dwork_query;
     double* dwork = NULL;
     int* iwork = NULL; // Needed only if JOB = 'S' or 'B'
     int iwork_size = 0;
 
     const int dico_len = 1, job_len = 1, fact_len = 1, trana_len = 1;
 
     char dico_upper = toupper(dico);
     char job_upper = toupper(job);
     char fact_upper = toupper(fact);
     char trana_upper = toupper(trana);
 
     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *u_cm = NULL, *c_cm = NULL;
 
     /* --- Input Parameter Validation --- */
 
     if (n < 0) { info = -5; goto cleanup; }
     if (dico_upper != 'C' && dico_upper != 'D') { info = -1; goto cleanup; }
     if (job_upper != 'X' && job_upper != 'S' && job_upper != 'B') { info = -2; goto cleanup; }
     if (fact_upper != 'F' && fact_upper != 'N') { info = -3; goto cleanup; }
     if (trana_upper != 'N' && trana_upper != 'T' && trana_upper != 'C') { info = -4; goto cleanup; }
 
     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, n);
     int min_ldu_f = MAX(1, n);
     int min_ldc_f = (job_upper == 'S') ? 1 : MAX(1, n);
 
     if (row_major) {
         // For row-major C, LD is the number of columns
         int min_lda_rm_cols = n;
         int min_ldu_rm_cols = n;
         int min_ldc_rm_cols = (job_upper == 'S') ? 1 : n;
         if (lda < min_lda_rm_cols) { info = -7; goto cleanup; }
         if (ldu < min_ldu_rm_cols) { info = -9; goto cleanup; }
         if (ldc < min_ldc_rm_cols) { info = -11; goto cleanup; }
     } else {
         // For column-major C, LD is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -7; goto cleanup; }
         if (ldu < min_ldu_f) { info = -9; goto cleanup; }
         if (ldc < min_ldc_f) { info = -11; goto cleanup; }
     }
 
     /* --- Workspace Allocation --- */
 
     // Allocate IWORK only if needed
     if (job_upper == 'S' || job_upper == 'B') {
         iwork_size = MAX(1, n * n);
         iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
         CHECK_ALLOC(iwork);
     } else {
         iwork = NULL; // Not referenced
     }
 
     // Allocate DWORK based on query
     ldwork = -1; // Query mode
     F77_FUNC(sb03md, SB03MD)(&dico_upper, &job_upper, &fact_upper, &trana_upper, &n,
                              a, &lda, u, &ldu, c, &ldc, scale, sep, ferr, wr, wi,
                              iwork, &dwork_query, &ldwork, &info,
                              dico_len, job_len, fact_len, trana_len);
 
     if (info < 0) { goto cleanup; } // Query failed due to invalid argument
     info = 0; // Reset info after query
 
     // Get the required dwork size from query result
     ldwork = (int)dwork_query;
     // Check against minimum documented size
     int min_ldwork = 1;
     if (job_upper == 'X') {
         if (fact_upper == 'F') min_ldwork = (dico_upper == 'C') ? MAX(1, n*n) : MAX(1, MAX(n*n, 2*n));
         else                   min_ldwork = MAX(1, MAX(n*n, 3*n));
     } else { // JOB = 'S' or 'B'
         if (fact_upper == 'F') min_ldwork = (dico_upper == 'C') ? MAX(1, 2*n*n) : MAX(1, 2*n*n + 2*n);
         else                   min_ldwork = (dico_upper == 'C') ? MAX(1, MAX(2*n*n, 3*n)) : MAX(1, 2*n*n + 2*n);
     }
     ldwork = MAX(ldwork, min_ldwork);
 
     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure
 
     /* --- Prepare Arrays and Call Fortran Routine --- */
     size_t elem_size = sizeof(double);
 
     // Determine sizes for potential copies
     size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
     size_t u_rows = n; size_t u_cols = n; size_t u_size = u_rows * u_cols;
     size_t c_rows = n; size_t c_cols = n; size_t c_size = c_rows * c_cols;
 
     if (row_major) {
         /* --- Row-Major Case --- */
 
         /* Allocate memory for column-major copies */
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (u_size > 0) { u_cm = (double*)malloc(u_size * elem_size); CHECK_ALLOC(u_cm); }
         if (c_size > 0 && (job_upper == 'X' || job_upper == 'B')) {
              c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm);
         }
 
         /* Transpose C inputs to Fortran copies */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         if (fact_upper == 'F' && u_size > 0) { // U is input only if FACT='F'
             slicot_transpose_to_fortran(u, u_cm, u_rows, u_cols, elem_size);
         }
         if (c_size > 0 && (job_upper == 'X' || job_upper == 'B')) { // C is input only if JOB='X' or 'B'
             slicot_transpose_symmetric_to_fortran(c, c_cm, n, 'U', elem_size); // Assuming UPLO='U' for C (docs say symmetric)
         }
 
         /* Fortran leading dimensions */
         int lda_f = (a_rows > 0) ? a_rows : 1;
         int ldu_f = (u_rows > 0) ? u_rows : 1;
         int ldc_f = (c_rows > 0) ? c_rows : 1;
 
         /* Call the Fortran routine */
         F77_FUNC(sb03md, SB03MD)(&dico_upper, &job_upper, &fact_upper, &trana_upper, &n,
                                  a_cm, &lda_f,           // Pass CM A (in/out if FACT='N')
                                  u_cm, &ldu_f,           // Pass CM U (in if FACT='F', out if FACT='N')
                                  (job_upper == 'S' ? NULL : c_cm), &ldc_f, // Pass CM C or NULL
                                  scale, sep, ferr, wr, wi,
                                  iwork, dwork, &ldwork, &info,
                                  dico_len, job_len, fact_len, trana_len);
 
         /* Copy back results from column-major temps to original row-major arrays */
         if (info == 0 || info == (n + 1)) {
             if (fact_upper == 'N' && a_size > 0) slicot_transpose_to_c(a_cm, a, a_rows, a_cols, elem_size); // Schur factor S
             if (fact_upper == 'N' && u_size > 0) slicot_transpose_to_c(u_cm, u, u_rows, u_cols, elem_size); // Schur vectors U
             if ((job_upper == 'X' || job_upper == 'B') && c_size > 0 && c_cm) {
                 slicot_transpose_symmetric_to_c(c_cm, c, n, 'U', elem_size); // Solution X
             }
             // SCALE, SEP, FERR, WR, WI modified directly
         }
         /* Temps freed in cleanup */
 
     } else {
         /* --- Column-Major Case --- */
         // Need to copy symmetric input C to ensure full matrix is passed if needed
         if (c_size > 0 && (job_upper == 'X' || job_upper == 'B')) {
             c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm);
             slicot_copy_symmetric_part(c, c_cm, n, 'U', ldc, elem_size); // Assuming UPLO='U' for C
         }
 
         /* Call the Fortran routine directly with user-provided arrays (or copy of C) */
         F77_FUNC(sb03md, SB03MD)(&dico_upper, &job_upper, &fact_upper, &trana_upper, &n,
                                  a, &lda,                // Pass original A
                                  u, &ldu,                // Pass original U
                                  (job_upper == 'S' ? NULL : (c_cm ? c_cm : c)), &ldc, // Pass copy of C or NULL
                                  scale, sep, ferr, wr, wi,
                                  iwork, dwork, &ldwork, &info,
                                  dico_len, job_len, fact_len, trana_len);
 
         // Copy back solution X from c_cm to c if column major copy was made
         if (info == 0 || info == (n + 1)) {
              if ((job_upper == 'X' || job_upper == 'B') && c_size > 0 && c_cm) {
                  slicot_copy_symmetric_part(c_cm, c, n, 'U', ldc, elem_size);
              }
         }
         // A, U, SCALE, SEP, FERR, WR, WI modified directly
     }
 
 cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(iwork); // Safe even if NULL
     free(a_cm);
     free(u_cm);
     free(c_cm);
 
     return info;
 }
 