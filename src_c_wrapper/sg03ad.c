/**
 * @file sg03ad.c
 * @brief C wrapper implementation for SLICOT routine SG03AD
 *
 * This file provides a C wrapper implementation for the SLICOT routine SG03AD,
 * which solves continuous- or discrete-time generalized Lyapunov
 * equations and optionally estimates the separation condition number.
 */

 #include <stdlib.h>
 #include <ctype.h>
 #include <stddef.h> // For size_t
 #include <string.h> // For memcpy
 
 // Include the header file for this wrapper
 #include "sg03ad.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines, slicot_copy_symmetric_part
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions
 
 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * A, E, Q, Z, X are input/output depending on flags.
  * SCALE, SEP, FERR, ALPHAR, ALPHAI, BETA are output.
  */
 extern void F77_FUNC(sg03ad, SG03AD)(
     const char* dico,       // CHARACTER*1 DICO
     const char* job,        // CHARACTER*1 JOB
     const char* fact,       // CHARACTER*1 FACT
     const char* trans,      // CHARACTER*1 TRANS
     const char* uplo,       // CHARACTER*1 UPLO
     const int* n,           // INTEGER N
     double* a,              // DOUBLE PRECISION A(LDA,*) (in/out)
     const int* lda,         // INTEGER LDA
     double* e,              // DOUBLE PRECISION E(LDE,*) (in/out)
     const int* lde,         // INTEGER LDE
     double* q,              // DOUBLE PRECISION Q(LDQ,*) (in/out)
     const int* ldq,         // INTEGER LDQ
     double* z,              // DOUBLE PRECISION Z(LDZ,*) (in/out)
     const int* ldz,         // INTEGER LDZ
     double* x,              // DOUBLE PRECISION X(LDX,*) (in/out)
     const int* ldx,         // INTEGER LDX
     double* scale,          // DOUBLE PRECISION SCALE (output)
     double* sep,            // DOUBLE PRECISION SEP (output)
     double* ferr,           // DOUBLE PRECISION FERR (output)
     double* alphar,         // DOUBLE PRECISION ALPHAR(*) (output)
     double* alphai,         // DOUBLE PRECISION ALPHAI(*) (output)
     double* beta,           // DOUBLE PRECISION BETA(*) (output)
     int* iwork,             // INTEGER IWORK(*)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* info,              // INTEGER INFO (output)
     int dico_len,           // Hidden length
     int job_len,            // Hidden length
     int fact_len,           // Hidden length
     int trans_len,          // Hidden length
     int uplo_len            // Hidden length
 );
 
 
 /* C wrapper function definition */
 int slicot_sg03ad(char dico, char job, char fact, char trans, char uplo, int n,
                   double* a, int lda, double* e, int lde,
                   double* q, int ldq, double* z, int ldz,
                   double* x, int ldx, double* scale, double* sep, double* ferr,
                   double* alphar, double* alphai, double* beta, int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = -1; /* Use -1 for workspace query */
     double dwork_query;
     double* dwork = NULL;
     int* iwork = NULL; // Needed only if JOB = 'S' or 'B'
     int iwork_size = 0;
 
     const int dico_len = 1, job_len = 1, fact_len = 1, trans_len = 1, uplo_len = 1;
 
     char dico_upper = toupper(dico);
     char job_upper = toupper(job);
     char fact_upper = toupper(fact);
     char trans_upper = toupper(trans);
     char uplo_upper = toupper(uplo);
 
     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *e_cm = NULL, *q_cm = NULL, *z_cm = NULL, *x_cm = NULL;
 
     /* --- Input Parameter Validation --- */
 
     if (n < 0) { info = -6; goto cleanup; }
     if (dico_upper != 'C' && dico_upper != 'D') { info = -1; goto cleanup; }
     if (job_upper != 'X' && job_upper != 'S' && job_upper != 'B') { info = -2; goto cleanup; }
     if (fact_upper != 'F' && fact_upper != 'N') { info = -3; goto cleanup; }
     if (trans_upper != 'N' && trans_upper != 'T') { info = -4; goto cleanup; }
     if (uplo_upper != 'U' && uplo_upper != 'L') { info = -5; goto cleanup; }
 
     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, n);
     int min_lde_f = MAX(1, n);
     int min_ldq_f = MAX(1, n);
     int min_ldz_f = MAX(1, n);
     int min_ldx_f = (job_upper == 'S') ? 1 : MAX(1, n);
 
     if (row_major) {
         // For row-major C, LD is the number of columns
         int min_lda_rm_cols = n;
         int min_lde_rm_cols = n;
         int min_ldq_rm_cols = n;
         int min_ldz_rm_cols = n;
         int min_ldx_rm_cols = (job_upper == 'S') ? 1 : n;
         if (lda < min_lda_rm_cols) { info = -8; goto cleanup; }
         if (lde < min_lde_rm_cols) { info = -10; goto cleanup; }
         if (ldq < min_ldq_rm_cols) { info = -12; goto cleanup; }
         if (ldz < min_ldz_rm_cols) { info = -14; goto cleanup; }
         if (ldx < min_ldx_rm_cols) { info = -16; goto cleanup; }
     } else {
         // For column-major C, LD is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -8; goto cleanup; }
         if (lde < min_lde_f) { info = -10; goto cleanup; }
         if (ldq < min_ldq_f) { info = -12; goto cleanup; }
         if (ldz < min_ldz_f) { info = -14; goto cleanup; }
         if (ldx < min_ldx_f) { info = -16; goto cleanup; }
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
     F77_FUNC(sg03ad, SG03AD)(&dico_upper, &job_upper, &fact_upper, &trans_upper, &uplo_upper, &n,
                              a, &lda, e, &lde, q, &ldq, z, &ldz, x, &ldx,
                              scale, sep, ferr, alphar, alphai, beta,
                              iwork, &dwork_query, &ldwork, &info,
                              dico_len, job_len, fact_len, trans_len, uplo_len);
 
     if (info < 0 && info != -24) { info = info; goto cleanup; } // Query failed due to invalid argument
     info = 0; // Reset info after query
 
     // Get the required dwork size from query result
     ldwork = (int)dwork_query;
     // Check against minimum documented size
     int min_ldwork = 1;
     if (job_upper == 'X') {
         min_ldwork = (fact_upper == 'F') ? MAX(1, n) : MAX(1, 4 * n);
     } else { // JOB = 'S' or 'B'
         min_ldwork = (fact_upper == 'F') ? MAX(1, 2*n*n) : MAX(1, MAX(2*n*n, 4*n));
     }
     ldwork = MAX(ldwork, min_ldwork);
 
     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure
 
     /* --- Prepare Arrays and Call Fortran Routine --- */
     size_t elem_size = sizeof(double);
 
     // Determine sizes for potential copies
     size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
     size_t e_rows = n; size_t e_cols = n; size_t e_size = e_rows * e_cols;
     size_t q_rows = n; size_t q_cols = n; size_t q_size = q_rows * q_cols;
     size_t z_rows = n; size_t z_cols = n; size_t z_size = z_rows * z_cols;
     size_t x_rows = n; size_t x_cols = n; size_t x_size = x_rows * x_cols;
 
     if (row_major) {
         /* --- Row-Major Case --- */
 
         /* Allocate memory for column-major copies */
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (e_size > 0) { e_cm = (double*)malloc(e_size * elem_size); CHECK_ALLOC(e_cm); }
         if (q_size > 0) { q_cm = (double*)malloc(q_size * elem_size); CHECK_ALLOC(q_cm); }
         if (z_size > 0) { z_cm = (double*)malloc(z_size * elem_size); CHECK_ALLOC(z_cm); }
         if (x_size > 0 && (job_upper == 'X' || job_upper == 'B')) {
              x_cm = (double*)malloc(x_size * elem_size); CHECK_ALLOC(x_cm);
         }
 
         /* Transpose C inputs to Fortran copies */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         if (e_size > 0) slicot_transpose_to_fortran(e, e_cm, e_rows, e_cols, elem_size);
         if (fact_upper == 'F' && q_size > 0) slicot_transpose_to_fortran(q, q_cm, q_rows, q_cols, elem_size);
         if (fact_upper == 'F' && z_size > 0) slicot_transpose_to_fortran(z, z_cm, z_rows, z_cols, elem_size);
         if (x_size > 0 && (job_upper == 'X' || job_upper == 'B')) {
             slicot_transpose_symmetric_to_fortran(x, x_cm, n, uplo_upper, elem_size); // Y is symmetric
         }
 
         /* Fortran leading dimensions */
         int lda_f = (a_rows > 0) ? a_rows : 1;
         int lde_f = (e_rows > 0) ? e_rows : 1;
         int ldq_f = (q_rows > 0) ? q_rows : 1;
         int ldz_f = (z_rows > 0) ? z_rows : 1;
         int ldx_f = (x_rows > 0) ? x_rows : 1;
 
         /* Call the Fortran routine */
         F77_FUNC(sg03ad, SG03AD)(&dico_upper, &job_upper, &fact_upper, &trans_upper, &uplo_upper, &n,
                                  a_cm, &lda_f, e_cm, &lde_f, q_cm, &ldq_f, z_cm, &ldz_f,
                                  (job_upper == 'S' ? NULL : x_cm), &ldx_f,
                                  scale, sep, ferr, alphar, alphai, beta,
                                  iwork, dwork, &ldwork, &info,
                                  dico_len, job_len, fact_len, trans_len, uplo_len);
 
         /* Copy back results from column-major temps to original row-major arrays */
         if (info == 0 || info == 3 || info == 4) {
             if (fact_upper == 'N' && a_size > 0) slicot_transpose_to_c(a_cm, a, a_rows, a_cols, elem_size); // A_s
             if (fact_upper == 'N' && e_size > 0) slicot_transpose_to_c(e_cm, e, e_rows, e_cols, elem_size); // E_s
             if (fact_upper == 'N' && q_size > 0) slicot_transpose_to_c(q_cm, q, q_rows, q_cols, elem_size); // Q vectors
             if (fact_upper == 'N' && z_size > 0) slicot_transpose_to_c(z_cm, z, z_rows, z_cols, elem_size); // Z vectors
             if ((job_upper == 'X' || job_upper == 'B') && x_size > 0 && x_cm) {
                 slicot_transpose_symmetric_to_c(x_cm, x, n, uplo_upper, elem_size); // Solution X
             }
             // SCALE, SEP, FERR, ALPHAR, ALPHAI, BETA modified directly
         }
         /* Temps freed in cleanup */
 
     } else {
         /* --- Column-Major Case --- */
         // Need to copy symmetric input X (used as Y) to ensure full matrix is passed if needed
         if (x_size > 0 && (job_upper == 'X' || job_upper == 'B')) {
             x_cm = (double*)malloc(x_size * elem_size); CHECK_ALLOC(x_cm);
             slicot_copy_symmetric_part(x, x_cm, n, uplo_upper, ldx, elem_size); // Copy Y
         }
 
         /* Call the Fortran routine directly with user-provided arrays (or copy of X) */
         F77_FUNC(sg03ad, SG03AD)(&dico_upper, &job_upper, &fact_upper, &trans_upper, &uplo_upper, &n,
                                  a, &lda, e, &lde, q, &ldq, z, &ldz,
                                  (job_upper == 'S' ? NULL : (x_cm ? x_cm : x)), &ldx, // Pass copy of Y or NULL
                                  scale, sep, ferr, alphar, alphai, beta,
                                  iwork, dwork, &ldwork, &info,
                                  dico_len, job_len, fact_len, trans_len, uplo_len);
 
         // Copy back solution X from x_cm to x if column major copy was made
         if (info == 0 || info == 3 || info == 4) {
              if ((job_upper == 'X' || job_upper == 'B') && x_size > 0 && x_cm) {
                  slicot_copy_symmetric_part(x_cm, x, n, uplo_upper, ldx, elem_size);
              }
         }
         // A, E, Q, Z, SCALE, SEP, FERR, ALPHAR, ALPHAI, BETA modified directly
     }
 
 cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(iwork); // Safe even if NULL
     free(a_cm);
     free(e_cm);
     free(q_cm);
     free(z_cm);
     free(x_cm); // Used for Y input copy in col-major, or X output copy in row-major
 
     return info;
 }
 