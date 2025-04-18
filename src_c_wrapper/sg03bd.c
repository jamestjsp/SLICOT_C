/**
 * @file sg03bd.c
 * @brief C wrapper implementation for SLICOT routine SG03BD
 *
 * This file provides a C wrapper implementation for the SLICOT routine SG03BD,
 * which solves stable continuous- or discrete-time generalized
 * Lyapunov equations for the Cholesky factor U of the solution
 * X = op(U)'*op(U).
 */

 #include <stdlib.h>
 #include <ctype.h>
 #include <stddef.h> // For size_t
 #include <stdio.h>  // For toupper if not in ctype.h
 #include <string.h> // For memcpy

 // Include the header file for this wrapper
 #include "sg03bd.h" // Assuming sg03bd.h exists
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, MIN, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * A, E, Q, Z, B are input/output depending on flags.
  * SCALE, ALPHAR, ALPHAI, BETA are output.
  */
 extern void F77_FUNC(sg03bd, SG03BD)(
     const char* dico,       // CHARACTER*1 DICO
     const char* fact,       // CHARACTER*1 FACT
     const char* trans,      // CHARACTER*1 TRANS
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     double* a,              // DOUBLE PRECISION A(LDA,*) (in/out)
     const int* lda,         // INTEGER LDA
     double* e,              // DOUBLE PRECISION E(LDE,*) (in/out)
     const int* lde,         // INTEGER LDE
     double* q,              // DOUBLE PRECISION Q(LDQ,*) (in/out)
     const int* ldq,         // INTEGER LDQ
     double* z,              // DOUBLE PRECISION Z(LDZ,*) (in/out)
     const int* ldz,         // INTEGER LDZ
     double* b,              // DOUBLE PRECISION B(LDB,*) (in/out -> U)
     const int* ldb,         // INTEGER LDB
     double* scale,          // DOUBLE PRECISION SCALE (output)
     double* alphar,         // DOUBLE PRECISION ALPHAR(*) (output)
     double* alphai,         // DOUBLE PRECISION ALPHAI(*) (output)
     double* beta,           // DOUBLE PRECISION BETA(*) (output)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* info,              // INTEGER INFO (output)
     int dico_len,           // Hidden length
     int fact_len,           // Hidden length
     int trans_len           // Hidden length
 );


 /* C wrapper function definition */
 int slicot_sg03bd(char dico, char fact, char trans, int n, int m,
                   double* a, int lda, double* e, int lde,
                   double* q, int ldq, double* z, int ldz,
                   double* b, int ldb, double* scale,
                   double* alphar, double* alphai, double* beta,
                   int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = -1; /* Use -1 for workspace query */
     double dwork_query;
     double* dwork = NULL;

     const int dico_len = 1, fact_len = 1, trans_len = 1;

     char dico_upper = toupper(dico);
     char fact_upper = toupper(fact);
     char trans_upper = toupper(trans);

     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *e_cm = NULL, *q_cm = NULL, *z_cm = NULL, *b_cm = NULL;

     /* --- Input Parameter Validation --- */

     if (n < 0) { info = -4; goto cleanup; }
     if (m < 0) { info = -5; goto cleanup; }
     if (dico_upper != 'C' && dico_upper != 'D') { info = -1; goto cleanup; }
     if (fact_upper != 'F' && fact_upper != 'N') { info = -2; goto cleanup; }
     if (trans_upper != 'N' && trans_upper != 'T') { info = -3; goto cleanup; }

     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, n);
     int min_lde_f = MAX(1, n);
     int min_ldq_f = MAX(1, n);
     int min_ldz_f = MAX(1, n);
     int min_ldb_f = (trans_upper == 'N') ? MAX(1, MAX(m, n)) : MAX(1, n);

     if (row_major) {
         // For row-major C, LD is the number of columns
         int min_lda_rm_cols = n;
         int min_lde_rm_cols = n;
         int min_ldq_rm_cols = n;
         int min_ldz_rm_cols = n;
         int min_ldb_rm_cols = (trans_upper == 'N') ? n : m;
         if (lda < min_lda_rm_cols) { info = -7; goto cleanup; }
         if (lde < min_lde_rm_cols) { info = -9; goto cleanup; }
         if (ldq < min_ldq_rm_cols) { info = -11; goto cleanup; }
         if (ldz < min_ldz_rm_cols) { info = -13; goto cleanup; }
         if (ldb < min_ldb_rm_cols) { info = -15; goto cleanup; }
     } else {
         // For column-major C, LD is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -7; goto cleanup; }
         if (lde < min_lde_f) { info = -9; goto cleanup; }
         if (ldq < min_ldq_f) { info = -11; goto cleanup; }
         if (ldz < min_ldz_f) { info = -13; goto cleanup; }
         if (ldb < min_ldb_f) { info = -15; goto cleanup; }
     }

     /* --- Workspace Allocation --- */

     // Allocate DWORK based on query
     ldwork = -1; // Query mode
     F77_FUNC(sg03bd, SG03BD)(&dico_upper, &fact_upper, &trans_upper, &n, &m,
                              a, &lda, e, &lde, q, &ldq, z, &ldz, b, &ldb,
                              scale, alphar, alphai, beta,
                              &dwork_query, &ldwork, &info,
                              dico_len, fact_len, trans_len);

     if (info < 0 && info != -21) { info = info; goto cleanup; } // Query failed due to invalid argument
     info = 0; // Reset info after query

     // Get the required dwork size from query result
     ldwork = (int)dwork_query;
     // Check against minimum documented size
     // FIX: Nest MAX calls to handle three arguments correctly
     int min_ldwork = (fact_upper == 'N') ? MAX(1, MAX(4*n, 6*n-6)) : MAX(1, MAX(2*n, 6*n-6));
     ldwork = MAX(ldwork, min_ldwork); // This MAX call is correct (two arguments)

     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure

     /* --- Prepare Arrays and Call Fortran Routine --- */
     size_t elem_size = sizeof(double);

     // Determine sizes for potential copies
     size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
     size_t e_rows = n; size_t e_cols = n; size_t e_size = e_rows * e_cols;
     size_t q_rows = n; size_t q_cols = n; size_t q_size = q_rows * q_cols;
     size_t z_rows = n; size_t z_cols = n; size_t z_size = z_rows * z_cols;
     size_t b_rows = (trans_upper == 'N') ? m : n;
     size_t b_cols = (trans_upper == 'N') ? n : m;
     size_t b_size = b_rows * b_cols;

     if (row_major) {
         /* --- Row-Major Case --- */

         /* Allocate memory for column-major copies */
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (e_size > 0) { e_cm = (double*)malloc(e_size * elem_size); CHECK_ALLOC(e_cm); }
         if (q_size > 0) { q_cm = (double*)malloc(q_size * elem_size); CHECK_ALLOC(q_cm); }
         if (z_size > 0) { z_cm = (double*)malloc(z_size * elem_size); CHECK_ALLOC(z_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }

         /* Transpose C inputs to Fortran copies */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         if (e_size > 0) slicot_transpose_to_fortran(e, e_cm, e_rows, e_cols, elem_size);
         if (fact_upper == 'F' && q_size > 0) slicot_transpose_to_fortran(q, q_cm, q_rows, q_cols, elem_size); // Copy Q if factored
         // Z is only used if FACT='F', copy it then.
         if (fact_upper == 'F' && z_size > 0) slicot_transpose_to_fortran(z, z_cm, z_rows, z_cols, elem_size);
         // B is input only if FACT='F'
         if (fact_upper == 'F' && b_size > 0) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, elem_size);

         /* Fortran leading dimensions */
         int lda_f = (a_rows > 0) ? a_rows : 1;
         int lde_f = (e_rows > 0) ? e_rows : 1;
         int ldq_f = (q_rows > 0) ? q_rows : 1;
         int ldz_f = (z_rows > 0) ? z_rows : 1;
         int ldb_f = (b_rows > 0) ? b_rows : 1;

         /* Call the Fortran routine */
         F77_FUNC(sg03bd, SG03BD)(&dico_upper, &fact_upper, &trans_upper, &n, &m,
                                  a_cm, &lda_f, e_cm, &lde_f,
                                  (fact_upper == 'F' ? q_cm : q), &ldq_f, // Pass original q if FACT='N'
                                  (fact_upper == 'F' ? z_cm : z), &ldz_f, // Pass original z if FACT='N'
                                  (fact_upper == 'F' ? b_cm : b), &ldb_f, // Pass original b if FACT='N' (b is output U then)
                                  scale, alphar, alphai, beta,
                                  dwork, &ldwork, &info,
                                  dico_len, fact_len, trans_len);

         /* Copy back results from column-major temps to original row-major arrays */
         // Copy back eigenvalues regardless of info for potential diagnostics
         // Copy back A_s, E_s, Q_vec, Z_vec if FACT='N' and successful/certain errors
         if (fact_upper == 'N' && (info == 0 || info == 1 || info == 3 || info == 5 || info == 6 || info == 7)) {
             if (a_size > 0) slicot_transpose_to_c(a_cm, a, a_rows, a_cols, elem_size); // A_s
             if (e_size > 0) slicot_transpose_to_c(e_cm, e, e_rows, e_cols, elem_size); // E_s
             if (q_size > 0) slicot_transpose_to_c(q_cm, q, q_rows, q_cols, elem_size); // Q vectors
             if (z_size > 0) slicot_transpose_to_c(z_cm, z, z_rows, z_cols, elem_size); // Z vectors
         }
         // Copy back solution U (which is in b_cm) only if successful or nearly singular
         if (info == 0 || info == 1) {
              // Output U is N x N, regardless of TRANS value
              size_t u_rows_out = n; size_t u_cols_out = n;
              // Ensure b_cm was allocated (it should be if b_size > 0)
              if (b_cm && (size_t)u_rows_out * u_cols_out > 0) {
                 slicot_transpose_to_c(b_cm, b, u_rows_out, u_cols_out, elem_size);
              }
         }
         // SCALE, ALPHAR, ALPHAI, BETA modified directly (passed by address)
         /* Temps freed in cleanup */

     } else {
         /* --- Column-Major Case --- */

         /* Call the Fortran routine directly with user-provided arrays */
         F77_FUNC(sg03bd, SG03BD)(&dico_upper, &fact_upper, &trans_upper, &n, &m,
                                  a, &lda, e, &lde, q, &ldq, z, &ldz, b, &ldb,
                                  scale, alphar, alphai, beta,
                                  dwork, &ldwork, &info,
                                  dico_len, fact_len, trans_len);
         // A, E, Q, Z, B, SCALE, ALPHAR, ALPHAI, BETA modified in place.
     }

 cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(a_cm);
     free(e_cm);
     free(q_cm);
     free(z_cm);
     free(b_cm);

     return info;
 }
