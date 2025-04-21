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
 #include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines, slicot_copy_symmetric_part, slicot_transpose_symmetric_to_fortran, slicot_transpose_symmetric_to_c
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
     double* a,              // DOUBLE PRECISION A(LDA,*) (in/out if FACT='N')
     const int* lda,         // INTEGER LDA
     double* e,              // DOUBLE PRECISION E(LDE,*) (in/out if FACT='N')
     const int* lde,         // INTEGER LDE
     double* q,              // DOUBLE PRECISION Q(LDQ,*) (in/out if FACT='N')
     const int* ldq,         // INTEGER LDQ
     double* z,              // DOUBLE PRECISION Z(LDZ,*) (in/out if FACT='N')
     const int* ldz,         // INTEGER LDZ
     double* x,              // DOUBLE PRECISION X(LDX,*) (in/out if JOB='X' or 'B')
     const int* ldx,         // INTEGER LDX
     double* scale,          // DOUBLE PRECISION SCALE (output)
     double* sep,            // DOUBLE PRECISION SEP (output)
     double* ferr,           // DOUBLE PRECISION FERR (output)
     double* alphar,         // DOUBLE PRECISION ALPHAR(*) (output if FACT='N')
     double* alphai,         // DOUBLE PRECISION ALPHAI(*) (output if FACT='N')
     double* beta,           // DOUBLE PRECISION BETA(*) (output if FACT='N')
     int* iwork,             // INTEGER IWORK(*) (workspace if JOB='S' or 'B')
     double* dwork,          // DOUBLE PRECISION DWORK(*) (workspace)
     const int* ldwork,      // INTEGER LDWORK
     int* info,              // INTEGER INFO (output)
     int dico_len,           // Hidden length
     int job_len,            // Hidden length
     int fact_len,           // Hidden length
     int trans_len,          // Hidden length
     int uplo_len            // Hidden length
 );


 /* C wrapper function definition */
 SLICOT_C_WRAPPER_API
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
     double *a_ptr, *e_ptr, *q_ptr, *z_ptr, *x_ptr;
     int lda_f, lde_f, ldq_f, ldz_f, ldx_f;


     /* --- Input Parameter Validation --- */

     if (dico_upper != 'C' && dico_upper != 'D') { info = -1; goto cleanup; }
     if (job_upper != 'X' && job_upper != 'S' && job_upper != 'B') { info = -2; goto cleanup; }
     if (fact_upper != 'F' && fact_upper != 'N') { info = -3; goto cleanup; }
     if (trans_upper != 'N' && trans_upper != 'T') { info = -4; goto cleanup; }
     if (uplo_upper != 'U' && uplo_upper != 'L') { info = -5; goto cleanup; }
     if (n < 0) { info = -6; goto cleanup; }

     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, n);
     int min_lde_f = MAX(1, n);
     int min_ldq_f = MAX(1, n);
     int min_ldz_f = MAX(1, n);
     int min_ldx_f = (job_upper == 'S') ? 1 : MAX(1, n); // X not referenced if JOB='S'

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

     /* --- Prepare arrays for column-major format if using row-major --- */
     size_t elem_size = sizeof(double);
     if (row_major) {
         /* Allocate memory for column-major copies */
         size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
         size_t e_rows = n; size_t e_cols = n; size_t e_size = e_rows * e_cols;
         size_t q_rows = n; size_t q_cols = n; size_t q_size = q_rows * q_cols;
         size_t z_rows = n; size_t z_cols = n; size_t z_size = z_rows * z_cols;
         size_t x_rows = n; size_t x_cols = n; size_t x_size = (job_upper != 'S') ? x_rows * x_cols : 0; // Size 0 if not referenced

         // Allocate based on usage (FACT, JOB flags)
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); } // A always in/out? Docs say in/out if FACT='N'
         if (e_size > 0) { e_cm = (double*)malloc(e_size * elem_size); CHECK_ALLOC(e_cm); } // E always in/out? Docs say in/out if FACT='N'
         if (fact_upper == 'F' && q_size > 0) { q_cm = (double*)malloc(q_size * elem_size); CHECK_ALLOC(q_cm); } // Q input if FACT='F', output if 'N'
         if (fact_upper == 'F' && z_size > 0) { z_cm = (double*)malloc(z_size * elem_size); CHECK_ALLOC(z_cm); } // Z input if FACT='F', output if 'N'
         if (x_size > 0) { x_cm = (double*)malloc(x_size * elem_size); CHECK_ALLOC(x_cm); } // X input if JOB!='S', output if JOB!='S'

         /* Transpose C inputs to Fortran copies */
         if (a_cm) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         if (e_cm) slicot_transpose_to_fortran(e, e_cm, e_rows, e_cols, elem_size);
         if (fact_upper == 'F' && q_cm) slicot_transpose_to_fortran(q, q_cm, q_rows, q_cols, elem_size);
         if (fact_upper == 'F' && z_cm) slicot_transpose_to_fortran(z, z_cm, z_rows, z_cols, elem_size);
         if (x_cm) slicot_transpose_symmetric_to_fortran(x, x_cm, n, uplo_upper, elem_size); // X (Y) is symmetric input

         /* Fortran leading dimensions */
         lda_f = (a_rows > 0) ? a_rows : 1;
         lde_f = (e_rows > 0) ? e_rows : 1;
         ldq_f = (q_rows > 0) ? q_rows : 1;
         ldz_f = (z_rows > 0) ? z_rows : 1;
         ldx_f = (x_rows > 0) ? x_rows : 1;

         /* Set pointers */
         a_ptr = a_cm;
         e_ptr = e_cm;
         q_ptr = (fact_upper == 'F') ? q_cm : q; // Use original Q if output (FACT='N')
         z_ptr = (fact_upper == 'F') ? z_cm : z; // Use original Z if output (FACT='N')
         x_ptr = (job_upper == 'S') ? NULL : x_cm; // Use copy of X (Y), or NULL if not referenced


     } else {
         /* Column-major case - use original arrays, potentially copy symmetric X (Y) */
         lda_f = lda;
         lde_f = lde;
         ldq_f = ldq;
         ldz_f = ldz;
         ldx_f = ldx;
         a_ptr = a;
         e_ptr = e;
         q_ptr = q;
         z_ptr = z;
         x_ptr = x;

         // Need to copy symmetric input X (used as Y) if JOB != 'S'
         size_t x_rows = n; size_t x_cols = n; size_t x_size = (job_upper != 'S') ? x_rows * x_cols : 0;
         if (x_size > 0) {
             x_cm = (double*)malloc(x_size * elem_size); CHECK_ALLOC(x_cm);
             slicot_copy_symmetric_part(x, x_cm, n, uplo_upper, ldx, elem_size); // Copy Y
             x_ptr = x_cm; // Use the temporary full copy
             ldx_f = n;    // Adjust LD for the full copy
         } else {
             x_ptr = NULL; // Pass NULL if JOB='S'
         }
     }


     /* --- Workspace Allocation --- */

     // Allocate IWORK only if needed
     if (job_upper == 'S' || job_upper == 'B') {
         iwork_size = MAX(1, n * n); // Ensure size >= 1
         iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
         CHECK_ALLOC(iwork);
     } else {
         iwork = NULL; // Not referenced
     }

     // Perform workspace query with small temporary array
     double dwork_temp[2]; // Small array to receive query result
     ldwork = -1; // Query mode
     
     F77_FUNC(sg03ad, SG03AD)(&dico_upper, &job_upper, &fact_upper, &trans_upper, &uplo_upper, &n,
                              a_ptr, &lda_f, e_ptr, &lde_f, q_ptr, &ldq_f, z_ptr, &ldz_f, x_ptr, &ldx_f,
                              scale, sep, ferr, alphar, alphai, beta,
                              iwork, dwork_temp, &ldwork, &info,
                              dico_len, job_len, fact_len, trans_len, uplo_len);
     
     if (info < 0 && info != -24) { goto cleanup; } // Query failed due to invalid argument (allow INFO=-24 from query)
     info = 0; // Reset info after query
     
     // Calculate minimum workspace based on requirements
     int min_ldwork;
     if (job_upper == 'X') {
         if (fact_upper == 'F') {
             min_ldwork = MAX(1, n);
         } else { // FACT = 'N'
             min_ldwork = MAX(1, 4 * n);
         }
     } else { // JOB = 'S' or 'B'
         if (fact_upper == 'F') {
             min_ldwork = MAX(1, 2 * n * n);
         } else { // FACT = 'N'
             min_ldwork = MAX(1, MAX(2 * n * n, 4 * n));
         }
     }
     
     // Use the larger of the queried optimal size and the minimum required size
     ldwork = (int)dwork_temp[0]; // First element contains optimal ldwork
     ldwork = MAX(ldwork, min_ldwork);

     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure


     /* --- Call the computational routine --- */
     F77_FUNC(sg03ad, SG03AD)(&dico_upper, &job_upper, &fact_upper, &trans_upper, &uplo_upper, &n,
                              a_ptr, &lda_f, e_ptr, &lde_f, q_ptr, &ldq_f, z_ptr, &ldz_f, x_ptr, &ldx_f,
                              scale, sep, ferr, alphar, alphai, beta,
                              iwork, dwork, &ldwork, &info,
                              dico_len, job_len, fact_len, trans_len, uplo_len);

     /* --- Copy results back to row-major format if needed --- */
     if (row_major && (info == 0 || info == 3 || info == 4)) { // Copy back based on success/certain errors
         size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
         size_t e_rows = n; size_t e_cols = n; size_t e_size = e_rows * e_cols;
         size_t q_rows = n; size_t q_cols = n; size_t q_size = q_rows * q_cols;
         size_t z_rows = n; size_t z_cols = n; size_t z_size = z_rows * z_cols;
         size_t x_rows = n; size_t x_cols = n; size_t x_size = (job_upper != 'S') ? x_rows * x_cols : 0;

         if (fact_upper == 'N' && a_cm) slicot_transpose_to_c(a_cm, a, a_rows, a_cols, elem_size); // A_s
         if (fact_upper == 'N' && e_cm) slicot_transpose_to_c(e_cm, e, e_rows, e_cols, elem_size); // E_s
         if (fact_upper == 'N' && q_size > 0) slicot_transpose_to_c(q_ptr, q, q_rows, q_cols, elem_size); // Q vectors (q_ptr points to original q)
         if (fact_upper == 'N' && z_size > 0) slicot_transpose_to_c(z_ptr, z, z_rows, z_cols, elem_size); // Z vectors (z_ptr points to original z)

         if (x_cm) { // Copy solution X (stored in x_cm) back to x
             slicot_transpose_symmetric_to_c(x_cm, x, n, uplo_upper, elem_size); // Solution X
         }
         // SCALE, SEP, FERR, ALPHAR, ALPHAI, BETA modified directly
     } else if (!row_major && (info == 0 || info == 3 || info == 4)) {
         // Copy back solution X from temporary full x_cm to original symmetric x
         size_t x_rows = n; size_t x_cols = n; size_t x_size = (job_upper != 'S') ? x_rows * x_cols : 0;
         if (x_cm) { // Check if x_cm was allocated (i.e., JOB != 'S')
             slicot_copy_symmetric_part(x_cm, x, n, uplo_upper, ldx, elem_size); // Copy full result to symmetric storage
         }
         // A, E, Q, Z (if FACT='N'), SCALE, SEP, FERR, ALPHAR, ALPHAI, BETA modified directly in original arrays.
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