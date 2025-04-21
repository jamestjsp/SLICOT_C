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
     double* a,              // DOUBLE PRECISION A(LDA,*) (in/out if FACT='N')
     const int* lda,         // INTEGER LDA
     double* e,              // DOUBLE PRECISION E(LDE,*) (in/out if FACT='N')
     const int* lde,         // INTEGER LDE
     double* q,              // DOUBLE PRECISION Q(LDQ,*) (in if FACT='F', out if FACT='N')
     const int* ldq,         // INTEGER LDQ
     double* z,              // DOUBLE PRECISION Z(LDZ,*) (in if FACT='F', out if FACT='N')
     const int* ldz,         // INTEGER LDZ
     double* b,              // DOUBLE PRECISION B(LDB,*) (in if FACT='F', out if FACT='N' -> U)
     const int* ldb,         // INTEGER LDB
     double* scale,          // DOUBLE PRECISION SCALE (output)
     double* alphar,         // DOUBLE PRECISION ALPHAR(*) (output if FACT='N')
     double* alphai,         // DOUBLE PRECISION ALPHAI(*) (output if FACT='N')
     double* beta,           // DOUBLE PRECISION BETA(*) (output if FACT='N')
     double* dwork,          // DOUBLE PRECISION DWORK(*) (workspace)
     const int* ldwork,      // INTEGER LDWORK
     int* info,              // INTEGER INFO (output)
     int dico_len,           // Hidden length
     int fact_len,           // Hidden length
     int trans_len           // Hidden length
 );


 /* C wrapper function definition */
 SLICOT_C_WRAPPER_API
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
     // No iwork needed

     const int dico_len = 1, fact_len = 1, trans_len = 1;

     char dico_upper = toupper(dico);
     char fact_upper = toupper(fact);
     char trans_upper = toupper(trans);

     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *e_cm = NULL, *q_cm = NULL, *z_cm = NULL, *b_cm = NULL;
     double *a_ptr, *e_ptr, *q_ptr, *z_ptr, *b_ptr;
     int lda_f, lde_f, ldq_f, ldz_f, ldb_f;


     /* --- Input Parameter Validation --- */

     if (dico_upper != 'C' && dico_upper != 'D') { info = -1; goto cleanup; }
     if (fact_upper != 'F' && fact_upper != 'N') { info = -2; goto cleanup; }
     if (trans_upper != 'N' && trans_upper != 'T') { info = -3; goto cleanup; }
     if (n < 0) { info = -4; goto cleanup; }
     if (m < 0) { info = -5; goto cleanup; }

     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, n);
     int min_lde_f = MAX(1, n);
     int min_ldq_f = MAX(1, n);
     int min_ldz_f = MAX(1, n);
     // B is N x M if TRANS='T', M x N if TRANS='N' (input for FACT='F')
     // B is N x N if TRANS='T', N x N if TRANS='N' (output U for FACT='N')
     int b_rows_in = (trans_upper == 'N') ? m : n;
     int b_cols_in = (trans_upper == 'N') ? n : m;
     int min_ldb_f_in = MAX(1, b_rows_in);
     int min_ldb_f_out = MAX(1, n); // Output U is always N x N
     int min_ldb_f = (fact_upper == 'F') ? min_ldb_f_in : min_ldb_f_out;

     if (row_major) {
         // For row-major C, LD is the number of columns
         int min_lda_rm_cols = n;
         int min_lde_rm_cols = n;
         int min_ldq_rm_cols = n;
         int min_ldz_rm_cols = n;
         int b_cols_rm_in = b_cols_in;
         int b_cols_rm_out = n;
         int min_ldb_rm_cols = (fact_upper == 'F') ? b_cols_rm_in : b_cols_rm_out;

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

     /* --- Prepare arrays for column-major format if using row-major --- */
     size_t elem_size = sizeof(double);

     // Determine sizes for potential copies
     size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
     size_t e_rows = n; size_t e_cols = n; size_t e_size = e_rows * e_cols;
     size_t q_rows = n; size_t q_cols = n; size_t q_size = q_rows * q_cols;
     size_t z_rows = n; size_t z_cols = n; size_t z_size = z_rows * z_cols;
     // Use input dimensions for B copy allocation if FACT='F'
     size_t b_copy_rows = (fact_upper == 'F') ? b_rows_in : n;
     size_t b_copy_cols = (fact_upper == 'F') ? b_cols_in : n;
     size_t b_copy_size = b_copy_rows * b_copy_cols;

     if (row_major) {
         /* Allocate memory for column-major copies */
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); } // A in/out if FACT='N'
         if (e_size > 0) { e_cm = (double*)malloc(e_size * elem_size); CHECK_ALLOC(e_cm); } // E in/out if FACT='N'
         if (fact_upper == 'F' && q_size > 0) { q_cm = (double*)malloc(q_size * elem_size); CHECK_ALLOC(q_cm); } // Q input if FACT='F', output if 'N'
         if (fact_upper == 'F' && z_size > 0) { z_cm = (double*)malloc(z_size * elem_size); CHECK_ALLOC(z_cm); } // Z input if FACT='F', output if 'N'
         if (b_copy_size > 0) { b_cm = (double*)malloc(b_copy_size * elem_size); CHECK_ALLOC(b_cm); } // B input if FACT='F', output if 'N'

         /* Transpose C inputs to Fortran copies */
         if (a_cm) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         if (e_cm) slicot_transpose_to_fortran(e, e_cm, e_rows, e_cols, elem_size);
         if (fact_upper == 'F' && q_cm) slicot_transpose_to_fortran(q, q_cm, q_rows, q_cols, elem_size);
         if (fact_upper == 'F' && z_cm) slicot_transpose_to_fortran(z, z_cm, z_rows, z_cols, elem_size);
         if (fact_upper == 'F' && b_cm) slicot_transpose_to_fortran(b, b_cm, b_rows_in, b_cols_in, elem_size);

         /* Fortran leading dimensions */
         lda_f = (a_rows > 0) ? a_rows : 1;
         lde_f = (e_rows > 0) ? e_rows : 1;
         ldq_f = (q_rows > 0) ? q_rows : 1;
         ldz_f = (z_rows > 0) ? z_rows : 1;
         ldb_f = (b_copy_rows > 0) ? b_copy_rows : 1; // Use copy rows for LD

         /* Set pointers */
         a_ptr = a_cm;
         e_ptr = e_cm;
         q_ptr = (fact_upper == 'F') ? q_cm : q; // Use original Q if output (FACT='N')
         z_ptr = (fact_upper == 'F') ? z_cm : z; // Use original Z if output (FACT='N')
         b_ptr = b_cm;                           // Use copy for B (input or output)

     } else {
         /* Column-major case - use original arrays */
         lda_f = lda;
         lde_f = lde;
         ldq_f = ldq;
         ldz_f = ldz;
         ldb_f = ldb;
         a_ptr = a;
         e_ptr = e;
         q_ptr = q;
         z_ptr = z;
         b_ptr = b;
     }


     /* --- Workspace Allocation --- */

     // Perform workspace query with small temporary array
     double dwork_temp[2]; // Small array to receive query result
     ldwork = -1; // Query mode
     
     F77_FUNC(sg03bd, SG03BD)(&dico_upper, &fact_upper, &trans_upper, &n, &m,
                              a_ptr, &lda_f, e_ptr, &lde_f, q_ptr, &ldq_f, z_ptr, &ldz_f, b_ptr, &ldb_f,
                              scale, alphar, alphai, beta,
                              dwork_temp, &ldwork, &info,
                              dico_len, fact_len, trans_len);

     if (info < 0 && info != -21) { goto cleanup; } // Query failed due to invalid argument (allow INFO=-21 from query)
     info = 0; // Reset info after query

     // Get the required dwork size from query result
     ldwork = (int)dwork_temp[0]; // First element contains optimal ldwork
     
     // Check against minimum documented size
     int min_ldwork = (fact_upper == 'N') ? MAX(1, MAX(4*n, 6*n-6)) : MAX(1, MAX(2*n, 6*n-6));
     ldwork = MAX(ldwork, min_ldwork);

     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure


     /* --- Call the computational routine --- */
     F77_FUNC(sg03bd, SG03BD)(&dico_upper, &fact_upper, &trans_upper, &n, &m,
                              a_ptr, &lda_f, e_ptr, &lde_f, q_ptr, &ldq_f, z_ptr, &ldz_f, b_ptr, &ldb_f,
                              scale, alphar, alphai, beta,
                              dwork, &ldwork, &info,
                              dico_len, fact_len, trans_len);

     /* --- Copy results back to row-major format if needed --- */
     if (row_major && (info == 0 || info == 1 || info == 3 || info == 5 || info == 6 || info == 7)) {
         size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
         size_t e_rows = n; size_t e_cols = n; size_t e_size = e_rows * e_cols;
         size_t q_rows = n; size_t q_cols = n; size_t q_size = q_rows * q_cols;
         size_t z_rows = n; size_t z_cols = n; size_t z_size = z_rows * z_cols;

         // Copy back A_s, E_s, Q_vec, Z_vec if FACT='N'
         if (fact_upper == 'N') {
             if (a_cm) slicot_transpose_to_c(a_cm, a, a_rows, a_cols, elem_size); // A_s
             if (e_cm) slicot_transpose_to_c(e_cm, e, e_rows, e_cols, elem_size); // E_s
             if (q_size > 0) slicot_transpose_to_c(q_ptr, q, q_rows, q_cols, elem_size); // Q vectors (q_ptr points to original q)
             if (z_size > 0) slicot_transpose_to_c(z_ptr, z, z_rows, z_cols, elem_size); // Z vectors (z_ptr points to original z)
         }

         // Copy back solution U (which is in b_cm) only if successful or nearly singular (info 0 or 1)
         if (info == 0 || info == 1) {
              size_t u_rows_out = n; size_t u_cols_out = n; // Output U is N x N
              size_t u_size_out = u_rows_out * u_cols_out;
              if (b_cm && u_size_out > 0) { // Check if b_cm was allocated
                 slicot_transpose_to_c(b_cm, b, u_rows_out, u_cols_out, elem_size);
              }
         }
         // SCALE, ALPHAR, ALPHAI, BETA modified directly
     }
     // In column-major case, A, E, Q, Z, B, SCALE, ALPHAR, ALPHAI, BETA modified in place.

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