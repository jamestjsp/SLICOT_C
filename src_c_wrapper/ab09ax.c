/**
 * @file ab09ax.c
 * @brief C wrapper implementation for SLICOT routine AB09AX
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB09AX,
 * which computes a reduced order model (Ar,Br,Cr) for a stable system
 * (A,B,C) using Balance & Truncate methods, where A is already in
 * real Schur form. It also returns the truncation matrices T and TI.
 */

 #include <stdlib.h>
 #include <ctype.h>
 #include <stddef.h> // For size_t
 
 // Include the header file for this wrapper
 #include "ab09ax.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions
 
 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * Note A, B, C are input/output. NR is input/output. T, TI are output.
  */
 extern void F77_FUNC(ab09ax, AB09AX)(
     const char* dico,       // CHARACTER*1 DICO
     const char* job,        // CHARACTER*1 JOB
     const char* ordsel,     // CHARACTER*1 ORDSEL
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     const int* p,           // INTEGER P
     int* nr,                // INTEGER NR (in/out)
     double* a,              // DOUBLE PRECISION A(LDA,*) (in/out)
     const int* lda,         // INTEGER LDA
     double* b,              // DOUBLE PRECISION B(LDB,*) (in/out)
     const int* ldb,         // INTEGER LDB
     double* c,              // DOUBLE PRECISION C(LDC,*) (in/out)
     const int* ldc,         // INTEGER LDC
     double* hsv,            // DOUBLE PRECISION HSV(*) (output)
     double* t,              // DOUBLE PRECISION T(LDT,*) (output)
     const int* ldt,         // INTEGER LDT
     double* ti,             // DOUBLE PRECISION TI(LDTI,*) (output)
     const int* ldti,        // INTEGER LDTI
     const double* tol,      // DOUBLE PRECISION TOL
     int* iwork,             // INTEGER IWORK(*)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* iwarn,             // INTEGER IWARN (output)
     int* info,              // INTEGER INFO (output)
     int dico_len,           // Hidden length
     int job_len,            // Hidden length
     int ordsel_len          // Hidden length
 );
 
 
 /* C wrapper function definition */
 SLICOT_C_WRAPPER_API
 int slicot_ab09ax(char dico, char job, char ordsel,
                   int n, int m, int p, int* nr,
                   double* a, int lda,
                   double* b, int ldb,
                   double* c, int ldc,
                   double* hsv,
                   double* t, int ldt,
                   double* ti, int ldti,
                   double tol, int* iwarn,
                   int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = -1; /* Use -1 for workspace query */
     double dwork_query;
     double* dwork = NULL;
     int* iwork = NULL;
     int iwork_size = 0;
     const int dico_len = 1, job_len = 1, ordsel_len = 1;
 
     char dico_upper = toupper(dico);
     char job_upper = toupper(job);
     char ordsel_upper = toupper(ordsel);
 
     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL;
     double *t_cm = NULL, *ti_cm = NULL; // For output matrices T, TI
 
     /* --- Input Parameter Validation --- */
 
     if (n < 0) { info = -4; goto cleanup; }
     if (m < 0) { info = -5; goto cleanup; }
     if (p < 0) { info = -6; goto cleanup; }
     if (dico_upper != 'C' && dico_upper != 'D') { info = -1; goto cleanup; }
     if (job_upper != 'B' && job_upper != 'N') { info = -2; goto cleanup; }
     if (ordsel_upper != 'F' && ordsel_upper != 'A') { info = -3; goto cleanup; }
     if (ordsel_upper == 'F' && (*nr < 0 || *nr > n) ) { info = -7; goto cleanup; }
 
     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, n);
     int min_ldb_f = MAX(1, n);
     int min_ldc_f = MAX(1, p);
     int min_ldt_f = MAX(1, n);
     // LDTI depends on NR, which might not be known yet if ORDSEL='A'.
     // Fortran checks LDTI >= MAX(1,NR). We check >= 1 initially.
     // A check after NR is computed might be needed if ORDSEL='A'.
     int min_ldti_f = 1; // Minimum check, Fortran does more specific check
 
     if (row_major) {
         // For row-major C, LDA is the number of columns
         int min_lda_rm_cols = n;
         int min_ldb_rm_cols = m;
         int min_ldc_rm_cols = n;
         int min_ldt_rm_cols = n; // T is N x NR -> needs N cols in C row-major
         int min_ldti_rm_cols = n; // TI is NR x N -> needs N cols in C row-major
         if (lda < min_lda_rm_cols) { info = -9; goto cleanup; }
         if (ldb < min_ldb_rm_cols) { info = -11; goto cleanup; }
         if (ldc < min_ldc_rm_cols) { info = -13; goto cleanup; }
         if (ldt < min_ldt_rm_cols) { info = -16; goto cleanup; }
         if (ldti < min_ldti_rm_cols) { info = -18; goto cleanup; }
     } else {
         // For column-major C, LDA is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -9; goto cleanup; }
         if (ldb < min_ldb_f) { info = -11; goto cleanup; }
         if (ldc < min_ldc_f) { info = -13; goto cleanup; }
         if (ldt < min_ldt_f) { info = -16; goto cleanup; }
         if (ldti < min_ldti_f) { info = -18; goto cleanup; }
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
     F77_FUNC(ab09ax, AB09AX)(&dico_upper, &job_upper, &ordsel_upper,
                              &n, &m, &p, nr, a, &lda, b, &ldb, c, &ldc,
                              hsv, t, &ldt, ti, &ldti, &tol,
                              iwork, &dwork_query, &ldwork,
                              iwarn, &info,
                              dico_len, job_len, ordsel_len);
 
     if (info != 0) { goto cleanup; } // Query failed
 
     // Get the required dwork size from query result
     ldwork = (int)dwork_query;
     // Check against minimum documented size: MAX(1,N*(MAX(N,M,P)+5) + N*(N+1)/2)
     int min_ldwork = 1;
     min_ldwork = MAX(min_ldwork, n * (MAX(n, MAX(m, p)) + 5) + n * (n + 1) / 2);
     ldwork = MAX(ldwork, min_ldwork);
 
     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure
 
     /* --- Prepare Arrays and Call Fortran Routine --- */
     size_t elem_size = sizeof(double);
     int nr_in = *nr; // Save input NR if needed for validation later
 
     if (row_major) {
         /* --- Row-Major Case --- */
 
         /* Allocate memory for column-major copies */
         size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
         size_t b_rows = n; size_t b_cols = m; size_t b_size = b_rows * b_cols;
         size_t c_rows = p; size_t c_cols = n; size_t c_size = c_rows * c_cols;
         // T(N, NR), TI(NR, N) - NR is output, use N for allocation size initially
         size_t t_rows = n; size_t t_cols = n; size_t t_size = t_rows * t_cols;
         size_t ti_rows = n; size_t ti_cols = n; size_t ti_size = ti_rows * ti_cols;
 
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
         // Allocate T and TI based on N initially, actual content is smaller
         if (t_size > 0) { t_cm = (double*)malloc(t_size * elem_size); CHECK_ALLOC(t_cm); }
         if (ti_size > 0) { ti_cm = (double*)malloc(ti_size * elem_size); CHECK_ALLOC(ti_cm); }
 
 
         /* Transpose C (row-major) inputs to Fortran (column-major) copies */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         if (b_size > 0) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, elem_size);
         if (c_size > 0) slicot_transpose_to_fortran(c, c_cm, c_rows, c_cols, elem_size);
 
         /* Fortran leading dimensions */
         int lda_f = (a_rows > 0) ? a_rows : 1;
         int ldb_f = (b_rows > 0) ? b_rows : 1;
         int ldc_f = (c_rows > 0) ? c_rows : 1;
         int ldt_f = (t_rows > 0) ? t_rows : 1; // LDT >= N
         int ldti_f = (ti_rows > 0) ? ti_rows: 1; // LDTI >= N (use N for allocation, Fortran requires >= NR)
 
         /* Call the Fortran routine */
         F77_FUNC(ab09ax, AB09AX)(&dico_upper, &job_upper, &ordsel_upper,
                                  &n, &m, &p, nr, // Pass pointer to nr
                                  a_cm, &lda_f,           // Pass CM A
                                  b_cm, &ldb_f,           // Pass CM B
                                  c_cm, &ldc_f,           // Pass CM C
                                  hsv,                    // Pass hsv pointer
                                  t_cm, &ldt_f,           // Pass CM T
                                  ti_cm, &ldti_f,         // Pass CM TI
                                  &tol, iwork, dwork, &ldwork,
                                  iwarn, &info,
                                  dico_len, job_len, ordsel_len);
 
         /* Copy back results from column-major temps to original row-major arrays */
         if (info == 0) {
             int nr_val = *nr; // Get the final reduced order
             // Check LDTI again now that NR is known
             if (ldti < MAX(1, nr_val)) {
                  info = -18; // LDTI is too small for computed NR
                  goto cleanup;
             }
 
             // Copy back reduced A, B, C and computed T, TI
             if (nr_val > 0) {
                 if (a_size > 0) slicot_transpose_to_c(a_cm, a, nr_val, nr_val, elem_size);
                 if (b_size > 0) slicot_transpose_to_c(b_cm, b, nr_val, m, elem_size);
                 if (c_size > 0) slicot_transpose_to_c(c_cm, c, p, nr_val, elem_size);
                 if (t_size > 0) slicot_transpose_to_c(t_cm, t, n, nr_val, elem_size); // T is N x NR
                 if (ti_size > 0) slicot_transpose_to_c(ti_cm, ti, nr_val, n, elem_size); // TI is NR x N
             }
             // HSV was filled directly.
         }
         /* Column-major copies will be freed in cleanup */
 
     } else {
         /* --- Column-Major Case --- */
 
         /* Call the Fortran routine directly with user-provided arrays */
         F77_FUNC(ab09ax, AB09AX)(&dico_upper, &job_upper, &ordsel_upper,
                                  &n, &m, &p, nr, // Pass pointer to nr
                                  a, &lda,                // Pass original A
                                  b, &ldb,                // Pass original B
                                  c, &ldc,                // Pass original C
                                  hsv,                    // Pass hsv pointer
                                  t, &ldt,                // Pass original T
                                  ti, &ldti,              // Pass original TI
                                  &tol, iwork, dwork, &ldwork,
                                  iwarn, &info,
                                  dico_len, job_len, ordsel_len);
         // A, B, C, HSV, T, TI are modified in place by the Fortran call.
     }
 
 cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(iwork); // Safe even if NULL
     free(a_cm);
     free(b_cm);
     free(c_cm);
     free(t_cm);
     free(ti_cm);
 
     return info;
 }
 