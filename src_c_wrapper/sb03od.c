/**
 * @file sb03od.c
 * @brief C wrapper implementation for SLICOT routine SB03OD
 *
 * This file provides a C wrapper implementation for the SLICOT routine SB03OD,
 * which solves stable continuous- or discrete-time Lyapunov equations
 * for the Cholesky factor U of the solution X = op(U)'*op(U).
 */

 #include <stdlib.h>
 #include <ctype.h>
 #include <stddef.h> // For size_t

 // Include the header file for this wrapper
 #include "sb03od.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, MIN, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * A, Q, B are input/output depending on flags.
  * SCALE, WR, WI are output.
  */
 extern void F77_FUNC(sb03od, SB03OD)(
     const char* dico,       // CHARACTER*1 DICO
     const char* fact,       // CHARACTER*1 FACT
     const char* trans,      // CHARACTER*1 TRANS
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     double* a,              // DOUBLE PRECISION A(LDA,*) (in/out if FACT='N')
     const int* lda,         // INTEGER LDA
     double* q,              // DOUBLE PRECISION Q(LDQ,*) (in if FACT='F', out if FACT='N')
     const int* ldq,         // INTEGER LDQ
     double* b,              // DOUBLE PRECISION B(LDB,*) (in/out -> U)
     const int* ldb,         // INTEGER LDB
     double* scale,          // DOUBLE PRECISION SCALE (output)
     double* wr,             // DOUBLE PRECISION WR(*) (output if FACT='N')
     double* wi,             // DOUBLE PRECISION WI(*) (output if FACT='N')
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* info,              // INTEGER INFO (output)
     int dico_len,           // Hidden length
     int fact_len,           // Hidden length
     int trans_len           // Hidden length
 );


 /* C wrapper function definition */
 SLICOT_EXPORT
 int slicot_sb03od(char dico, char fact, char trans, int n, int m,
                   double* a, int lda, double* q, int ldq,
                   double* b, int ldb, double* scale,
                   double* wr, double* wi, int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = -1; /* Use -1 for workspace query */
     double dwork_query;
     double* dwork = NULL;
     // No iwork needed for this routine

     const int dico_len = 1, fact_len = 1, trans_len = 1;

     char dico_upper = toupper(dico);
     char fact_upper = toupper(fact);
     char trans_upper = toupper(trans);

     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *q_cm = NULL, *b_cm = NULL;
     double *a_ptr, *q_ptr, *b_ptr;
     int lda_f, ldq_f, ldb_f;

     /* --- Input Parameter Validation --- */

     if (dico_upper != 'C' && dico_upper != 'D') { info = -1; goto cleanup; }
     if (fact_upper != 'F' && fact_upper != 'N') { info = -2; goto cleanup; }
     if (trans_upper != 'N' && trans_upper != 'T') { info = -3; goto cleanup; }
     if (n < 0) { info = -4; goto cleanup; }
     if (m < 0) { info = -5; goto cleanup; }


     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, n);
     int min_ldq_f = MAX(1, n);
     int min_ldb_f = (trans_upper == 'N') ? MAX(1, MAX(m, n)) : MAX(1, n); // B stores U (NxN) or B(MxN)

     if (row_major) {
         // For row-major C, LD is the number of columns
         int min_lda_rm_cols = n;
         int min_ldq_rm_cols = n;
         int min_ldb_rm_cols = (trans_upper == 'N') ? n : m; // B stores U (N,N) or B(M,N)
         if (lda < min_lda_rm_cols) { info = -7; goto cleanup; }
         if (ldq < min_ldq_rm_cols) { info = -9; goto cleanup; }
         if (ldb < min_ldb_rm_cols) { info = -11; goto cleanup; }
     } else {
         // For column-major C, LD is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -7; goto cleanup; }
         if (ldq < min_ldq_f) { info = -9; goto cleanup; }
         if (ldb < min_ldb_f) { info = -11; goto cleanup; }
     }

     /* --- Prepare arrays for column-major format if using row-major --- */
     size_t elem_size = sizeof(double);
     if (row_major) {
         /* Allocate memory for column-major copies */
         size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
         size_t q_rows = n; size_t q_cols = n; size_t q_size = q_rows * q_cols;
         size_t b_rows_in = (trans_upper == 'N') ? m : n;
         size_t b_cols_in = (trans_upper == 'N') ? n : m;
         size_t b_size_in = b_rows_in * b_cols_in; // Input B size

         // Allocate based on usage (FACT flags)
         if (fact_upper == 'N' && a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); } // A in/out if FACT='N'
         else if (a_size > 0)                 { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); } // A input if FACT='F'
         if (q_size > 0)                      { q_cm = (double*)malloc(q_size * elem_size); CHECK_ALLOC(q_cm); } // Q in/out
         if (b_size_in > 0)                   { b_cm = (double*)malloc(b_size_in * elem_size); CHECK_ALLOC(b_cm); } // B in/out -> U


         /* Transpose C inputs to Fortran copies */
         if (a_cm) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         if (fact_upper == 'F' && q_cm) { // Q is input only if FACT='F'
             slicot_transpose_to_fortran(q, q_cm, q_rows, q_cols, elem_size);
         }
         if (b_cm) slicot_transpose_to_fortran(b, b_cm, b_rows_in, b_cols_in, elem_size);


         /* Fortran leading dimensions */
         lda_f = (a_rows > 0) ? a_rows : 1;
         ldq_f = (q_rows > 0) ? q_rows : 1;
         ldb_f = (b_rows_in > 0) ? b_rows_in : 1; // Use rows of input B for Fortran LD

         /* Set pointers */
         a_ptr = a_cm;
         q_ptr = q_cm;
         b_ptr = b_cm;

     } else {
         /* Column-major case - use original arrays */
         lda_f = lda;
         ldq_f = ldq;
         ldb_f = ldb;
         a_ptr = a;
         q_ptr = q;
         b_ptr = b;
     }


     /* --- Workspace Allocation --- */

     // Perform workspace query for DWORK
     double dwork_temp[1]; // Small array to receive query result
     ldwork = -1; // Query mode
     F77_FUNC(sb03od, SB03OD)(&dico_upper, &fact_upper, &trans_upper, &n, &m,
                              a_ptr, &lda_f, q_ptr, &ldq_f, b_ptr, &ldb_f,
                              scale, wr, wi,
                              dwork_temp, &ldwork, &info,
                              dico_len, fact_len, trans_len);

     if (info < 0 && info != -16) { goto cleanup; } // Query failed due to invalid argument (allow INFO=-16 from query)
     info = 0; // Reset info after query

     // Get the required dwork size from query result
     ldwork = (int)dwork_temp[0]; // First element contains optimal ldwork
     // Check against minimum documented size
     int min_ldwork = (m > 0) ? MAX(1, 4 * n) : 1; // If m=0, B is not referenced, min ldwork is 1
     ldwork = MAX(ldwork, min_ldwork);

     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure

     /* --- Call the computational routine --- */
     F77_FUNC(sb03od, SB03OD)(&dico_upper, &fact_upper, &trans_upper, &n, &m,
                              a_ptr, &lda_f,           // Pass A ptr
                              q_ptr, &ldq_f,           // Pass Q ptr
                              b_ptr, &ldb_f,           // Pass B ptr
                              scale, wr, wi,          // Pass output pointers
                              dwork, &ldwork, &info,
                              dico_len, fact_len, trans_len);

     /* --- Copy results back to row-major format if needed --- */
     if (row_major && (info == 0 || info == 1)) { // Copy back even if INFO=1 (nearly singular)
         size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
         size_t q_rows = n; size_t q_cols = n; size_t q_size = q_rows * q_cols;
         // Result U is N x N upper triangular, stored in B/b_cm
         size_t u_rows_out = n; size_t u_cols_out = n;
         size_t b_rows_in = (trans_upper == 'N') ? m : n; // Original input rows of B
         size_t b_cols_in = (trans_upper == 'N') ? n : m; // Original input cols of B
         size_t b_size_in = b_rows_in * b_cols_in; // Original input size


         if (fact_upper == 'N' && a_cm) slicot_transpose_to_c(a_cm, a, a_rows, a_cols, elem_size); // Schur factor S
         if (fact_upper == 'N' && q_cm) slicot_transpose_to_c(q_cm, q, q_rows, q_cols, elem_size); // Schur vectors Q
         if (b_cm) {
             slicot_transpose_to_c(b_cm, b, u_rows_out, u_cols_out, elem_size); // Copy NxN result U
         }
         // SCALE, WR, WI modified directly
     }
     // In column-major case, A, Q, B, SCALE, WR, WI are modified in place.

 cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(a_cm);
     free(q_cm);
     free(b_cm);

     return info;
 }