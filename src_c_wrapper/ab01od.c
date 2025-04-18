/**
 * @file ab01od.c
 * @brief C wrapper implementation for SLICOT routine AB01OD
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB01OD,
 * which reduces a state-space system (A,B) to staircase form using
 * orthogonal state and input transformations.
 */

 #include <stdlib.h>
 #include <string.h> // For memset
 #include <ctype.h>  // For toupper
 
 // Include the header file for this wrapper
 #include "ab01od.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines, set_identity
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions
 
 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * This handles potential name mangling issues between C and Fortran compilers.
  * All arguments are passed by reference (pointers).
  * Hidden string length arguments for character arguments are explicitly included.
  */
 extern void F77_FUNC(ab01od, AB01OD)(
     const char* stages,     // CHARACTER*1 STAGES
     const char* jobu,       // CHARACTER*1 JOBU
     const char* jobv,       // CHARACTER*1 JOBV
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     double* a,              // DOUBLE PRECISION A(LDA,*) (in/out)
     const int* lda,         // INTEGER LDA
     double* b,              // DOUBLE PRECISION B(LDB,*) (in/out)
     const int* ldb,         // INTEGER LDB
     double* u,              // DOUBLE PRECISION U(LDU,*) (in/out)
     const int* ldu,         // INTEGER LDU
     double* v,              // DOUBLE PRECISION V(LDV,*) (output)
     const int* ldv,         // INTEGER LDV
     int* ncont,             // INTEGER NCONT (in/out)
     int* indcon,            // INTEGER INDCON (in/out)
     int* kstair,            // INTEGER KSTAIR(*) (in/out)
     const double* tol,      // DOUBLE PRECISION TOL
     int* iwork,             // INTEGER IWORK(*)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* info,              // INTEGER INFO (output)
     int stages_len,         // Hidden length for stages
     int jobu_len,           // Hidden length for jobu
     int jobv_len            // Hidden length for jobv
 );
 
 
 /* C wrapper function definition */
 SLICOT_C_WRAPPER_API
 int slicot_ab01od(char stages, char jobu, char jobv,
                   int n, int m,
                   double* a, int lda,
                   double* b, int ldb,
                   double* u, int ldu,
                   double* v, int ldv,
                   int* ncont, int* indcon, int* kstair,
                   double tol, int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = -1; /* Use -1 for workspace query */
     double dwork_query;
     double* dwork = NULL;
     int* iwork = NULL;
     int iwork_size = 0;
     const int stages_len = 1;
     const int jobu_len = 1;
     const int jobv_len = 1;
 
     char stages_upper = toupper(stages);
     char jobu_upper = toupper(jobu);
     char jobv_upper = toupper(jobv);
 
     /* Pointers for column-major copies if needed */
     double *a_cm = NULL;
     double *b_cm = NULL;
     double *u_cm = NULL; // Needed if row_major and jobu == 'I'
     double *v_cm = NULL; // Needed if row_major and jobv == 'I' and stages != 'F'
 
     /* --- Input Parameter Validation --- */
 
     if (n < 0) { info = -4; goto cleanup; }
     if (m < 0) { info = -5; goto cleanup; }
 
     if (stages_upper != 'F' && stages_upper != 'B' && stages_upper != 'A') {
         info = -1; goto cleanup;
     }
     if (jobu_upper != 'N' && jobu_upper != 'I') {
         info = -2; goto cleanup;
     }
     if (jobv_upper != 'N' && jobv_upper != 'I') {
         info = -3; goto cleanup;
     }
 
     // Check leading dimensions based on storage order and flags
     int min_lda_f = MAX(1, n);
     int min_ldb_f = MAX(1, n);
     int min_ldu_f = (jobu_upper == 'N') ? 1 : MAX(1, n);
     int min_ldv_f = (jobv_upper == 'N' || stages_upper == 'F') ? 1 : MAX(1, m);
 
     if (row_major) {
         // For row-major C, LDA/LDB/LDU/LDV is the number of columns
         int min_lda_rm_cols = n;
         int min_ldb_rm_cols = m;
         int min_ldu_rm_cols = (jobu_upper == 'N') ? 1 : n;
         int min_ldv_rm_cols = (jobv_upper == 'N' || stages_upper == 'F') ? 1 : m;
         if (lda < min_lda_rm_cols) { info = -7; goto cleanup; }
         if (ldb < min_ldb_rm_cols) { info = -9; goto cleanup; }
         if (jobu_upper == 'I' && ldu < min_ldu_rm_cols) { info = -11; goto cleanup; }
         if (jobv_upper == 'I' && stages_upper != 'F' && ldv < min_ldv_rm_cols) { info = -13; goto cleanup; }
     } else {
         // For column-major C, LDA/LDB/LDU/LDV is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -7; goto cleanup; }
         if (ldb < min_ldb_f) { info = -9; goto cleanup; }
         if (ldu < min_ldu_f) { info = -11; goto cleanup; } // Check even if JOBU='N'
         if (ldv < min_ldv_f) { info = -13; goto cleanup; } // Check even if JOBV='N' or STAGES='F'
     }
 
     /* --- Workspace Allocation --- */
 
     // Allocate IWORK (size M) only if needed
     if (stages_upper != 'B') {
         iwork_size = MAX(1, m); // Ensure minimum size 1
         iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
         CHECK_ALLOC(iwork);
     } else {
         iwork = NULL; // Not referenced if STAGES = 'B'
     }
 
 
     // Allocate DWORK based on query
     ldwork = -1; // Query mode
     F77_FUNC(ab01od, AB01OD)(&stages_upper, &jobu_upper, &jobv_upper,
                              &n, &m, a, &lda, b, &ldb, u, &ldu, v, &ldv,
                              ncont, indcon, kstair, &tol,
                              iwork, // Pass NULL if STAGES == 'B'
                              &dwork_query, &ldwork, &info,
                              stages_len, jobu_len, jobv_len);
 
     if (info != 0) {
         // Query failed
         goto cleanup;
     }
 
     // Get the required dwork size from query result
     ldwork = (int)dwork_query;
     // Check against minimum documented size based on STAGES
     int min_ldwork;
     if (stages_upper != 'B') {
         min_ldwork = MAX(1, n + MAX(n, 3 * m));
     } else {
         min_ldwork = MAX(1, m + MAX(n, m));
     }
     ldwork = MAX(ldwork, min_ldwork);
 
     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure
 
     /* --- Prepare Arrays and Call Fortran Routine --- */
 
     if (row_major) {
         /* --- Row-Major Case --- */
 
         /* Allocate memory for column-major copies */
         size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
         size_t b_rows = n; size_t b_cols = m; size_t b_size = b_rows * b_cols;
         size_t u_rows = n; size_t u_cols = n; size_t u_size = u_rows * u_cols;
         size_t v_rows = m; size_t v_cols = m; size_t v_size = v_rows * v_cols;
         int ldu_f = 1;
         int ldv_f = 1;
 
         if (a_size > 0) { a_cm = (double*)malloc(a_size * sizeof(double)); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * sizeof(double)); CHECK_ALLOC(b_cm); }
 
         if (jobu_upper == 'I') {
             ldu_f = (u_rows > 0) ? u_rows : 1;
             if (u_size > 0) { u_cm = (double*)malloc(u_size * sizeof(double)); CHECK_ALLOC(u_cm); }
             if (stages_upper == 'B') {
                 // U is input, transpose it
                 if (u_size > 0) slicot_transpose_to_fortran(u, u_cm, u_rows, u_cols, sizeof(double));
             } else {
                 // U is output, initialize to identity
                  if (u_size > 0) set_identity(n, u_cm, ldu_f, 0); // 0 for column-major
             }
         }
 
         if (jobv_upper == 'I' && stages_upper != 'F') {
             ldv_f = (v_rows > 0) ? v_rows : 1;
             if (v_size > 0) { v_cm = (double*)malloc(v_size * sizeof(double)); CHECK_ALLOC(v_cm); }
             // V is output, initialize to identity
             if (v_size > 0) set_identity(m, v_cm, ldv_f, 0); // 0 for column-major
         }
 
         /* Transpose C (row-major) inputs A and B to Fortran (column-major) copies */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, sizeof(double));
         if (b_size > 0) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, sizeof(double));
 
         /* Fortran leading dimensions */
         int lda_f = (a_rows > 0) ? a_rows : 1;
         int ldb_f = (b_rows > 0) ? b_rows : 1;
         // ldu_f, ldv_f calculated above
 
         /* Call the Fortran routine */
         F77_FUNC(ab01od, AB01OD)(&stages_upper, &jobu_upper, &jobv_upper,
                                  &n, &m,
                                  a_cm, &lda_f,           // Pass CM A
                                  b_cm, &ldb_f,           // Pass CM B
                                  (jobu_upper == 'N' ? NULL : u_cm), // Pass CM U or NULL
                                  &ldu_f,                 // Pass Fortran LDU
                                  (jobv_upper == 'N' || stages_upper == 'F' ? NULL : v_cm), // Pass CM V or NULL
                                  &ldv_f,                 // Pass Fortran LDV
                                  ncont, indcon, kstair,  // Pass output pointers
                                  &tol,                   // Pass address of tol
                                  iwork, dwork, &ldwork, &info, // Pass workspaces
                                  stages_len, jobu_len, jobv_len); // Pass hidden lengths
 
         /* Copy back results */
         if (info == 0) {
             // Transpose modified a_cm back to original row-major a
             if (a_size > 0) slicot_transpose_to_c(a_cm, a, a_rows, a_cols, sizeof(double));
             // Transpose modified b_cm back to original row-major b
             if (b_size > 0) slicot_transpose_to_c(b_cm, b, b_rows, b_cols, sizeof(double));
             // Transpose computed/updated u_cm back to original row-major u
             if (jobu_upper == 'I' && u_size > 0) {
                 slicot_transpose_to_c(u_cm, u, u_rows, u_cols, sizeof(double));
             }
             // Transpose computed v_cm back to original row-major v
             if (jobv_upper == 'I' && stages_upper != 'F' && v_size > 0) {
                 slicot_transpose_to_c(v_cm, v, v_rows, v_cols, sizeof(double));
             }
             // NCONT, INDCON, KSTAIR were modified in place (passed directly).
         }
         /* Column-major copies will be freed in cleanup */
 
     } else {
         /* --- Column-Major Case --- */
 
         /* Initialize U and V if requested */
         if (jobu_upper == 'I' && stages_upper != 'B') {
             set_identity(n, u, ldu, 0); // 0 for column-major
         }
         if (jobv_upper == 'I' && stages_upper != 'F') {
             set_identity(m, v, ldv, 0); // 0 for column-major
         }
 
         /* Call the Fortran routine directly with user-provided arrays */
         F77_FUNC(ab01od, AB01OD)(&stages_upper, &jobu_upper, &jobv_upper,
                                  &n, &m,
                                  a, &lda,                // Pass original A
                                  b, &ldb,                // Pass original B
                                  (jobu_upper == 'N' ? NULL : u), // Pass original U or NULL
                                  &ldu,                   // Pass original LDU
                                  (jobv_upper == 'N' || stages_upper == 'F' ? NULL : v), // Pass original V or NULL
                                  &ldv,                   // Pass original LDV
                                  ncont, indcon, kstair,  // Pass output pointers
                                  &tol,                   // Pass address of tol
                                  iwork, dwork, &ldwork, &info, // Pass workspaces
                                  stages_len, jobu_len, jobv_len); // Pass hidden lengths
         // A, B, U, V, NCONT, INDCON, KSTAIR are modified in place by the Fortran call.
     }
 
 cleanup:
     /* --- Cleanup --- */
     // Free allocated memory (free(NULL) is safe)
     free(dwork);
     free(iwork); // iwork might be NULL if STAGES == 'B'
     // Free column-major copies if they were allocated
     free(a_cm);
     free(b_cm);
     free(u_cm);
     free(v_cm);
 
     /* Return the info code from the Fortran routine or SLICOT_MEMORY_ERROR */
     return info;
 }
 