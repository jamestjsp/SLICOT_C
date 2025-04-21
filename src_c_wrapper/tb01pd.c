/**
 * @file tb01pd.c
 * @brief C wrapper implementation for SLICOT routine TB01PD
 *
 * This file provides a C wrapper implementation for the SLICOT routine TB01PD,
 * which finds a reduced (minimal, controllable, or observable)
 * state-space representation (Ar,Br,Cr) in upper block Hessenberg form
 * for a given system (A,B,C).
 */

 #include <stdlib.h>
 #include <ctype.h>
 #include <stddef.h> // For size_t
 #include <string.h> // For memcpy if needed

 // Include the header file for this wrapper
 #include "tb01pd.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * A, B, C are input/output. NR is output. IWORK contains output info.
  */
 extern void F77_FUNC(tb01pd, TB01PD)(
     const char* job,        // CHARACTER*1 JOB
     const char* equil,      // CHARACTER*1 EQUIL
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     const int* p,           // INTEGER P
     double* a,              // DOUBLE PRECISION A(LDA,*) (in/out)
     const int* lda,         // INTEGER LDA
     double* b,              // DOUBLE PRECISION B(LDB,*) (in/out)
     const int* ldb,         // INTEGER LDB
     double* c,              // DOUBLE PRECISION C(LDC,*) (in/out)
     const int* ldc,         // INTEGER LDC
     int* nr,                // INTEGER NR (output)
     const double* tol,      // DOUBLE PRECISION TOL (in)
     int* iwork,             // INTEGER IWORK(*) (output info/workspace)
     double* dwork,          // DOUBLE PRECISION DWORK(*) (workspace)
     const int* ldwork,      // INTEGER LDWORK
     int* info,              // INTEGER INFO (output)
     int job_len,            // Hidden length
     int equil_len           // Hidden length
 );


 /* C wrapper function definition */
 SLICOT_C_WRAPPER_API
 int slicot_tb01pd(char job, char equil, int n, int m, int p,
                   double* a, int lda, double* b, int ldb,
                   double* c, int ldc, int* nr, double tol,
                   int* iwork, /* Now treated as workspace/output */
                   int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = -1; /* Use -1 for workspace query */
     double dwork_query;
     double* dwork = NULL;
     // IWORK is now managed externally as workspace/output. No allocation here.
     const int job_len = 1, equil_len = 1;

     char job_upper = toupper(job);
     char equil_upper = toupper(equil);

     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL;
     double *a_ptr, *b_ptr, *c_ptr;
     int lda_f, ldb_f, ldc_f;

     /* --- Input Parameter Validation --- */

     if (job_upper != 'M' && job_upper != 'C' && job_upper != 'O') { info = -1; goto cleanup; }
     if (equil_upper != 'S' && equil_upper != 'N') { info = -2; goto cleanup; }
     if (n < 0) { info = -3; goto cleanup; }
     if (m < 0) { info = -4; goto cleanup; }
     if (p < 0) { info = -5; goto cleanup; }
     // TOL check done by Fortran

     // Check leading dimensions based on storage order and JOB
     int max_mp = MAX(m, p);
     int min_lda_f = MAX(1, n);
     int min_ldb_f = MAX(1, n);
     // LDC needs space for C and potentially workspace related to B if JOB='C'
     int c_fort_rows = (job_upper == 'C') ? max_mp : p;
     int min_ldc_f = MAX(1, c_fort_rows);

     // B needs space for B and potentially workspace related to C if JOB='O'
     int b_fort_cols = (job_upper == 'O') ? max_mp : m;

     if (row_major) {
         // For row-major C, LD is the number of columns
         int min_lda_rm_cols = n;
         int min_ldb_rm_cols = b_fort_cols;
         int min_ldc_rm_cols = n;
         if (lda < min_lda_rm_cols) { info = -7; goto cleanup; }
         if (ldb < min_ldb_rm_cols) { info = -9; goto cleanup; }
         if (ldc < min_ldc_rm_cols) { info = -11; goto cleanup; }
     } else {
         // For column-major C, LD is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -7; goto cleanup; }
         if (ldb < min_ldb_f) { info = -9; goto cleanup; }
         if (ldc < min_ldc_f) { info = -11; goto cleanup; }
     }

     /* --- Prepare arrays for column-major format if using row-major --- */
     size_t elem_size = sizeof(double);
     if (row_major) {
         /* Allocate memory for column-major copies */
         size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
         // Allocate B copy based on required Fortran columns
         size_t b_rows = n; size_t b_work_cols = b_fort_cols; size_t b_work_size = b_rows * b_work_cols;
         // Allocate C copy based on required Fortran rows
         size_t c_work_rows = c_fort_rows; size_t c_cols = n; size_t c_work_size = c_work_rows * c_cols;

         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_work_size > 0) { b_cm = (double*)malloc(b_work_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_work_size > 0) { c_cm = (double*)malloc(c_work_size * elem_size); CHECK_ALLOC(c_cm); }

         /* Transpose C inputs to Fortran copies */
         if (a_cm) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         // Copy only actual N x M part of B into potentially larger workspace
         if (b_cm && n > 0 && m > 0) slicot_transpose_to_fortran(b, b_cm, n, m, elem_size);
         // Copy only actual P x N part of C into potentially larger workspace
         if (c_cm && p > 0 && n > 0) slicot_transpose_to_fortran(c, c_cm, p, n, elem_size);

         /* Fortran leading dimensions */
         lda_f = (a_rows > 0) ? a_rows : 1;
         ldb_f = (b_rows > 0) ? b_rows : 1;
         ldc_f = (c_work_rows > 0) ? c_work_rows : 1; // Use workspace rows for LDC_F

         /* Set pointers */
         a_ptr = a_cm;
         b_ptr = b_cm;
         c_ptr = c_cm;

     } else {
         /* Column-major case - use original arrays */
         lda_f = lda;
         ldb_f = ldb;
         ldc_f = ldc;
         a_ptr = a;
         b_ptr = b;
         c_ptr = c;
     }

     /* --- Workspace Allocation --- */

     // Perform workspace query for DWORK
     ldwork = -1; // Query mode
     F77_FUNC(tb01pd, TB01PD)(&job_upper, &equil_upper, &n, &m, &p,
                              a_ptr, &lda_f, b_ptr, &ldb_f, c_ptr, &ldc_f, nr, &tol,
                              iwork, &dwork_query, &ldwork, &info,
                              job_len, equil_len);

     if (info < 0 && info != -16) { goto cleanup; } // Query failed due to invalid argument (allow INFO=-16)
     info = 0; // Reset info after query

     // Get the required dwork size from query result
     ldwork = (int)dwork_query;
     // Check against minimum documented size: MAX(1, N + MAX(N, 3*M, 3*P))
     int min_ldwork = 1;
     // Properly nest MAX for 3 arguments
     min_ldwork = MAX(min_ldwork, n + MAX(n, MAX(3 * m, 3 * p)));
     ldwork = MAX(ldwork, min_ldwork);

     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure


     /* --- Call the computational routine --- */
     F77_FUNC(tb01pd, TB01PD)(&job_upper, &equil_upper, &n, &m, &p,
                              a_ptr, &lda_f,           // Pass A ptr
                              b_ptr, &ldb_f,           // Pass B ptr (original or workspace)
                              c_ptr, &ldc_f,           // Pass C ptr (original or workspace)
                              nr, &tol, iwork, dwork, &ldwork, &info,
                              job_len, equil_len);

     /* --- Copy results back to row-major format if needed --- */
     if (row_major && info == 0) {
         int nr_val = *nr;
         size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
         size_t b_rows = n; size_t b_cols = m; size_t b_size = b_rows * b_cols; // Actual B size
         size_t c_rows = p; size_t c_cols = n; size_t c_size = c_rows * c_cols; // Actual C size

         if (nr_val >= 0) { // Allow nr = 0
             if (a_cm) slicot_transpose_to_c(a_cm, a, nr_val, nr_val, elem_size); // Ar
             if (b_cm && nr_val > 0 && m > 0) slicot_transpose_to_c(b_cm, b, nr_val, m, elem_size);       // Br
             if (c_cm && p > 0 && nr_val > 0) slicot_transpose_to_c(c_cm, c, p, nr_val, elem_size);       // Cr
         }
         // NR, IWORK modified directly
     }
     // In column-major case, A, B, C, NR, IWORK modified in place.

 cleanup:
     /* --- Cleanup --- */
     free(dwork);
     // IWORK is external workspace/output, not freed here
     free(a_cm);
     free(b_cm);
     free(c_cm);

     return info;
 }