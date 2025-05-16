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
     // double dwork_query; // Not used if query result directly to dwork_temp
     double* dwork = NULL;
     // No iwork needed for this routine
 
     const int dico_len = 1, fact_len = 1, trans_len = 1;
 
     char dico_upper = toupper(dico);
     char fact_upper = toupper(fact);
     char trans_upper = toupper(trans);
 
     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *q_cm = NULL, *b_cm = NULL;
     double *a_ptr, *q_ptr, *b_ptr;
     int lda_f_fortran, ldq_f_fortran, ldb_f_fortran; // Fortran leading dimensions
 
     /* --- Input Parameter Validation --- */
     /* Moved to the beginning */
     if (dico_upper != 'C' && dico_upper != 'D') { info = -1; goto cleanup; }
     if (fact_upper != 'F' && fact_upper != 'N') { info = -2; goto cleanup; }
     if (trans_upper != 'N' && trans_upper != 'T') { info = -3; goto cleanup; }
     if (n < 0) { info = -4; goto cleanup; }
     if (m < 0) { info = -5; goto cleanup; }
 
     // Add NULL checks for essential pointers
     // These checks are now before any potential use of a, q, b, scale, wr, wi
     if (n > 0 && a == NULL) { info = -6; goto cleanup; }
     if (n > 0 && q == NULL) { info = -8; goto cleanup; }
     if (((m > 0 && n > 0) || (n > 0 && trans_upper == 'N') || (n > 0 && trans_upper == 'T')) && b == NULL) { info = -10; goto cleanup; }
     if (scale == NULL) { info = -12; goto cleanup; }
     if (n > 0 && wr == NULL) { info = -13; goto cleanup; }
     if (n > 0 && wi == NULL) { info = -14; goto cleanup; }
 
     // Check leading dimensions based on storage order
     int min_lda_f_val = MAX(1, n);
     int min_ldq_f_val = MAX(1, n);
     int min_ldb_f_val = (trans_upper == 'N') ? MAX(1, MAX(m, n)) : MAX(1, n);
 
     if (row_major) {
         // For row-major C, LD is the number of columns
         int min_lda_rm_cols = n;
         int min_ldq_rm_cols = n;
         // For input B (row-major): M x N if TRANS='N', N x M if TRANS='T'
         // For output U (row-major): N x N
         // The 'ldb' from C is num_cols of B_io.
         // If TRANS='N', input B is M rows, N cols. ldb (C) should be >= N.
         // If TRANS='T', input B is N rows, M cols. ldb (C) should be >= M.
         // The wrapper copies U (NxN) back to B_io, so B_io's ldb must accommodate N cols for U.
         int min_ldb_rm_cols_input = (trans_upper == 'N') ? n : m;
         int min_ldb_rm_cols_output_u = n;
         // The ldb passed by the user must be valid for their input B matrix.
         // The wrapper will ensure B_io can hold output U.
         if (lda < min_lda_rm_cols && n > 0) { info = -7; goto cleanup; }
         if (ldq < min_ldq_rm_cols && n > 0) { info = -9; goto cleanup; }
         if (ldb < min_ldb_rm_cols_input && ((trans_upper == 'N' && n > 0) || (trans_upper == 'T' && m > 0))) { info = -11; goto cleanup; }
         if (ldb < min_ldb_rm_cols_output_u && n > 0) { /* Also ensure ldb can hold output U */ info = -11; goto cleanup; }


     } else {
         // For column-major C, LD is the number of rows (Fortran style)
         if (lda < min_lda_f_val) { info = -7; goto cleanup; }
         if (ldq < min_ldq_f_val) { info = -9; goto cleanup; }
         if (ldb < min_ldb_f_val) { info = -11; goto cleanup; }
     }
 
     /* --- Prepare arrays for column-major format if using row-major --- */
     size_t elem_size = sizeof(double);
     if (row_major) {
         /* Calculate Fortran leading dimensions and column counts for _cm arrays */
         lda_f_fortran = MAX(1, n);
         ldq_f_fortran = MAX(1, n);
         ldb_f_fortran = (trans_upper == 'N') ? MAX(1, MAX(m, n)) : MAX(1, n);
 
         size_t a_cm_rows = lda_f_fortran; size_t a_cm_cols = MAX(1,n);
         size_t q_cm_rows = ldq_f_fortran; size_t q_cm_cols = MAX(1,n);
         size_t b_cm_rows = ldb_f_fortran;
         size_t b_cm_alloc_cols = (trans_upper == 'N') ? MAX(1,n) : MAX(1,MAX(n,m));
         if (n == 0 && trans_upper == 'T' && m == 0) b_cm_alloc_cols = 1; // Ensure min 1x1 if all zero
         else if (n == 0 && trans_upper == 'N') b_cm_alloc_cols = 1;
 
         size_t a_size = a_cm_rows * a_cm_cols;
         size_t q_size = q_cm_rows * q_cm_cols;
         size_t b_size = b_cm_rows * b_cm_alloc_cols;
 
         // Allocate based on usage (FACT flags)
         if (n > 0) { // A and Q are N x N
             if (fact_upper == 'N') { // A is in/out, Q is out
                 a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm);
                 q_cm = (double*)malloc(q_size * elem_size); CHECK_ALLOC(q_cm);
             } else { // FACT == 'F', A is in, Q is in
                 a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm);
                 q_cm = (double*)malloc(q_size * elem_size); CHECK_ALLOC(q_cm);
             }
         }
         // B is M x N or N x M input, N x N output (U)
         // Allocate b_cm if b is not NULL (it can be NULL if M=0 and N=0, though checks prevent N<0, M<0)
         // The Fortran routine might not touch B if M=0.
         // If n=0, b_ptr can be NULL. If m=0 and n>0, b_ptr is used for output U (zero).
         if (b_size > 0 && b != NULL) { 
              b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm);
         } else if (b_size == 0 && n == 0 && m == 0 && b != NULL) { 
            // Case where n=0, m=0, b_size might be 1*1 due to MAX(1,...) in alloc_cols.
            // If b is not NULL, and b_size is calculated > 0, allocate.
            // This specific condition might be too narrow, b_size > 0 && b != NULL should cover it.
         }


         /* Transpose C inputs to Fortran copies */
         if (a_cm && a) slicot_transpose_to_fortran_with_ld(a, a_cm, n, n, lda, lda_f_fortran, elem_size);
         if (q_cm && q) { // Q is input if FACT='F', Q is output if FACT='N' (but not transposed back if FACT='F')
             if (fact_upper == 'F') { // Q is input only if FACT='F'
                  slicot_transpose_to_fortran_with_ld(q, q_cm, n, n, ldq, ldq_f_fortran, elem_size);
             }
             // If FACT='N', Q is output, so q_cm is just an empty buffer for Fortran.
             // It will be transposed back later.
         }
         if (b_cm && b) {
             size_t b_math_rows_in = (trans_upper == 'N') ? m : n;
             size_t b_math_cols_in = (trans_upper == 'N') ? n : m;
             if (m > 0 || (trans_upper == 'T' && n > 0)) { // Only transpose if B has input data
                  slicot_transpose_to_fortran_with_ld(b, b_cm, b_math_rows_in, b_math_cols_in, ldb, ldb_f_fortran, elem_size);
             }
         }
 
          /* Set pointers */
          a_ptr = a_cm;
          q_ptr = q_cm;
          b_ptr = b_cm;
 
      } else {
          /* Column-major case - use original arrays */
          lda_f_fortran = lda;
          ldq_f_fortran = ldq;
          ldb_f_fortran = ldb;
          a_ptr = a;
          q_ptr = q;
          b_ptr = b;
      }
 
 
      /* --- Workspace Allocation --- */
 
      // Perform workspace query for DWORK
      double dwork_temp[1]; // Small array to receive query result
      ldwork = -1; // Query mode
      F77_FUNC(sb03od, SB03OD)(&dico_upper, &fact_upper, &trans_upper, &n, &m,
                              a_ptr, &lda_f_fortran, q_ptr, &ldq_f_fortran, b_ptr, &ldb_f_fortran, // Args 6-11
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
                              a_ptr, &lda_f_fortran,           
                              q_ptr, &ldq_f_fortran,           
                              b_ptr, &ldb_f_fortran,           
                               scale, wr, wi,          // Pass output pointers
                               dwork, &ldwork, &info,
                               dico_len, fact_len, trans_len);
 
      /* --- Copy results back to row-major format if needed --- */
      if (row_major && (info == 0 || info == 1)) { // Copy back even if INFO=1 (nearly singular)
         if (fact_upper == 'N' && a_cm && a) slicot_transpose_to_c_with_ld(a_cm, a, n, n, lda_f_fortran, lda, elem_size); // Schur factor S
         if (fact_upper == 'N' && q_cm && q) slicot_transpose_to_c_with_ld(q_cm, q, n, n, ldq_f_fortran, ldq, elem_size); // Schur vectors Q
          if (b_cm && b && n > 0) { // U is N x N
             slicot_transpose_to_c_with_ld(b_cm, b, n, n, ldb_f_fortran, ldb, elem_size); // Copy NxN result U
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