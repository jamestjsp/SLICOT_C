/**
 * @file sb02md.c
 * @brief C wrapper implementation for SLICOT routine SB02MD
 *
 * This file provides a C wrapper implementation for the SLICOT routine SB02MD,
 * which solves continuous- or discrete-time algebraic Riccati
 * equations using the Schur vectors method.
 */

 #include <stdlib.h>
 #include <ctype.h>
 #include <stddef.h> // For size_t
 #include <string.h> // For memcpy
 
 // Include the header file for this wrapper
 #include "sb02md.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines, slicot_copy_symmetric_part
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions
 
 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * A is input/output if DICO='D'. G is input. Q is input/output (solution X).
  * RCOND, WR, WI, S, U are output.
  */
 extern void F77_FUNC(sb02md, SB02MD)(
     const char* dico,       // CHARACTER*1 DICO
     const char* hinv,       // CHARACTER*1 HINV
     const char* uplo,       // CHARACTER*1 UPLO
     const char* scal,       // CHARACTER*1 SCAL
     const char* sort,       // CHARACTER*1 SORT
     const int* n,           // INTEGER N
     double* a,              // DOUBLE PRECISION A(LDA,*) (in/out if DICO='D')
     const int* lda,         // INTEGER LDA
     const double* g,        // DOUBLE PRECISION G(LDG,*)
     const int* ldg,         // INTEGER LDG
     double* q,              // DOUBLE PRECISION Q(LDQ,*) (in/out -> X)
     const int* ldq,         // INTEGER LDQ
     double* rcond,          // DOUBLE PRECISION RCOND (output)
     double* wr,             // DOUBLE PRECISION WR(*) (output)
     double* wi,             // DOUBLE PRECISION WI(*) (output)
     double* s,              // DOUBLE PRECISION S(LDS,*) (output)
     const int* lds,         // INTEGER LDS
     double* u,              // DOUBLE PRECISION U(LDU,*) (output)
     const int* ldu,         // INTEGER LDU
     int* iwork,             // INTEGER IWORK(*)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* bwork,             // LOGICAL BWORK(*) -> int*
     int* info,              // INTEGER INFO (output)
     int dico_len,           // Hidden length
     int hinv_len,           // Hidden length
     int uplo_len,           // Hidden length
     int scal_len,           // Hidden length
     int sort_len            // Hidden length
 );
 
 
 /* C wrapper function definition */
 int slicot_sb02md(char dico, char hinv, char uplo, char scal, char sort,
                   int n, double* a, int lda, const double* g, int ldg,
                   double* q, int ldq, double* rcond,
                   double* wr, double* wi,
                   double* s, int lds, double* u, int ldu,
                   int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = -1; /* Use -1 for workspace query */
     double dwork_query;
     double* dwork = NULL;
     int* iwork = NULL;
     int* bwork = NULL; // Map LOGICAL to int
     int work_size_2n = 0;
 
     const int dico_len = 1, hinv_len = 1, uplo_len = 1, scal_len = 1, sort_len = 1;
 
     char dico_upper = toupper(dico);
     char hinv_upper = toupper(hinv);
     char uplo_upper = toupper(uplo);
     char scal_upper = toupper(scal);
     char sort_upper = toupper(sort);
 
     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *g_cm = NULL, *q_cm = NULL;
     double *s_cm = NULL, *u_cm = NULL;
 
     /* --- Input Parameter Validation --- */
 
     if (n < 0) { info = -6; goto cleanup; }
     if (dico_upper != 'C' && dico_upper != 'D') { info = -1; goto cleanup; }
     if (dico_upper == 'D' && hinv_upper != 'D' && hinv_upper != 'I') { info = -2; goto cleanup; }
     if (uplo_upper != 'U' && uplo_upper != 'L') { info = -3; goto cleanup; }
     if (scal_upper != 'G' && scal_upper != 'N') { info = -4; goto cleanup; }
     if (sort_upper != 'S' && sort_upper != 'U') { info = -5; goto cleanup; }
 
     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, n);
     int min_ldg_f = MAX(1, n);
     int min_ldq_f = MAX(1, n);
     int min_lds_f = MAX(1, 2 * n);
     int min_ldu_f = MAX(1, 2 * n);
 
     if (row_major) {
         // For row-major C, LDA is the number of columns
         int min_lda_rm_cols = n;
         int min_ldg_rm_cols = n;
         int min_ldq_rm_cols = n;
         int min_lds_rm_cols = 2 * n;
         int min_ldu_rm_cols = 2 * n;
         if (lda < min_lda_rm_cols) { info = -8; goto cleanup; }
         if (ldg < min_ldg_rm_cols) { info = -10; goto cleanup; }
         if (ldq < min_ldq_rm_cols) { info = -12; goto cleanup; }
         if (lds < min_lds_rm_cols) { info = -17; goto cleanup; }
         if (ldu < min_ldu_rm_cols) { info = -19; goto cleanup; }
     } else {
         // For column-major C, LDA is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -8; goto cleanup; }
         if (ldg < min_ldg_f) { info = -10; goto cleanup; }
         if (ldq < min_ldq_f) { info = -12; goto cleanup; }
         if (lds < min_lds_f) { info = -17; goto cleanup; }
         if (ldu < min_ldu_f) { info = -19; goto cleanup; }
     }
 
     /* --- Workspace Allocation --- */
 
     // Allocate IWORK and BWORK (size 2*N)
     work_size_2n = MAX(1, 2 * n);
     iwork = (int*)malloc((size_t)work_size_2n * sizeof(int));
     CHECK_ALLOC(iwork);
     bwork = (int*)malloc((size_t)work_size_2n * sizeof(int)); // Use int for LOGICAL
     CHECK_ALLOC(bwork);
 
     // Allocate DWORK based on query
     ldwork = -1; // Query mode
     F77_FUNC(sb02md, SB02MD)(&dico_upper, &hinv_upper, &uplo_upper, &scal_upper, &sort_upper,
                              &n, a, &lda, g, &ldg, q, &ldq, rcond, wr, wi,
                              s, &lds, u, &ldu, iwork, &dwork_query, &ldwork,
                              bwork, &info,
                              dico_len, hinv_len, uplo_len, scal_len, sort_len);
 
     if (info < 0) { goto cleanup; } // Query failed due to invalid argument
     info = 0; // Reset info after query
 
     // Get the required dwork size from query result
     ldwork = (int)dwork_query;
     // Check against minimum documented size
     int min_ldwork = (dico_upper == 'C') ? MAX(2, 6 * n) : MAX(3, 6 * n);
     ldwork = MAX(ldwork, min_ldwork);
 
     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure
 
     /* --- Prepare Arrays and Call Fortran Routine --- */
     size_t elem_size = sizeof(double);
 
     // Determine sizes for potential copies (moved before if/else)
     size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
     size_t g_rows = n; size_t g_cols = n; size_t g_size = g_rows * g_cols;
     size_t q_rows = n; size_t q_cols = n; size_t q_size = q_rows * q_cols;
     size_t s_rows = 2*n; size_t s_cols = 2*n; size_t s_size = s_rows * s_cols;
     size_t u_rows = 2*n; size_t u_cols = 2*n; size_t u_size = u_rows * u_cols;
 
     if (row_major) {
         /* --- Row-Major Case --- */
 
         /* Allocate memory for column-major copies */
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (g_size > 0) { g_cm = (double*)malloc(g_size * elem_size); CHECK_ALLOC(g_cm); }
         if (q_size > 0) { q_cm = (double*)malloc(q_size * elem_size); CHECK_ALLOC(q_cm); }
         if (s_size > 0) { s_cm = (double*)malloc(s_size * elem_size); CHECK_ALLOC(s_cm); }
         if (u_size > 0) { u_cm = (double*)malloc(u_size * elem_size); CHECK_ALLOC(u_cm); }
 
         /* Transpose C (row-major) inputs to Fortran (column-major) copies */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         if (g_size > 0) slicot_transpose_symmetric_to_fortran(g, g_cm, n, uplo_upper, elem_size); // Use symmetric transpose
         if (q_size > 0) slicot_transpose_symmetric_to_fortran(q, q_cm, n, uplo_upper, elem_size); // Use symmetric transpose
 
         /* Fortran leading dimensions */
         int lda_f = (a_rows > 0) ? a_rows : 1;
         int ldg_f = (g_rows > 0) ? g_rows : 1;
         int ldq_f = (q_rows > 0) ? q_rows : 1;
         int lds_f = (s_rows > 0) ? s_rows : 1;
         int ldu_f = (u_rows > 0) ? u_rows : 1;
 
         /* Call the Fortran routine */
         F77_FUNC(sb02md, SB02MD)(&dico_upper, &hinv_upper, &uplo_upper, &scal_upper, &sort_upper,
                                  &n,
                                  a_cm, &lda_f,           // Pass CM A (in/out if D='D')
                                  g_cm, &ldg_f,           // Pass CM G (in)
                                  q_cm, &ldq_f,           // Pass CM Q (in/out -> X)
                                  rcond, wr, wi,          // Pass output pointers
                                  s_cm, &lds_f,           // Pass CM S (out)
                                  u_cm, &ldu_f,           // Pass CM U (out)
                                  iwork, dwork, &ldwork, bwork, &info,
                                  dico_len, hinv_len, uplo_len, scal_len, sort_len);
 
         /* Copy back results from column-major temps to original row-major arrays */
         if (info == 0 || info == 5) { // Copy back even if INFO=5
             // Copy back A only if discrete time case
             if (dico_upper == 'D' && a_size > 0) {
                 slicot_transpose_to_c(a_cm, a, a_rows, a_cols, elem_size);
             }
             // Copy back Q (contains solution X)
             if (q_size > 0) slicot_transpose_symmetric_to_c(q_cm, q, n, uplo_upper, elem_size); // Symmetric copy back
             // Copy back S and U
             if (s_size > 0) slicot_transpose_to_c(s_cm, s, s_rows, s_cols, elem_size);
             if (u_size > 0) slicot_transpose_to_c(u_cm, u, u_rows, u_cols, elem_size);
             // RCOND, WR, WI are filled directly.
         }
         /* Column-major copies will be freed in cleanup */
 
     } else {
         /* --- Column-Major Case --- */
         // Need to copy symmetric inputs G and Q to ensure full matrix is passed if needed
         if (g_size > 0) { g_cm = (double*)malloc(g_size * elem_size); CHECK_ALLOC(g_cm); slicot_copy_symmetric_part(g, g_cm, n, uplo_upper, ldg, elem_size); }
         if (q_size > 0) { q_cm = (double*)malloc(q_size * elem_size); CHECK_ALLOC(q_cm); slicot_copy_symmetric_part(q, q_cm, n, uplo_upper, ldq, elem_size); }
 
 
         /* Call the Fortran routine directly with user-provided arrays (or copies) */
         F77_FUNC(sb02md, SB02MD)(&dico_upper, &hinv_upper, &uplo_upper, &scal_upper, &sort_upper,
                                  &n,
                                  a, &lda,                // Pass original A
                                  (g_cm ? g_cm : g), &ldg, // Pass copy of G if created
                                  (q_cm ? q_cm : q), &ldq, // Pass copy of Q (in/out -> X) if created
                                  rcond, wr, wi,          // Pass output pointers
                                  s, &lds,                // Pass original S
                                  u, &ldu,                // Pass original U
                                  iwork, dwork, &ldwork, bwork, &info,
                                  dico_len, hinv_len, uplo_len, scal_len, sort_len);
 
         // Copy back solution X from q_cm to q if column major copy was made
         if (info == 0 || info == 5) {
              if (q_size > 0 && q_cm) slicot_copy_symmetric_part(q_cm, q, n, uplo_upper, ldq, elem_size);
         }
         // A (if DICO='D'), RCOND, WR, WI, S, U are modified in place.
     }
 
 cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(iwork);
     free(bwork);
     free(a_cm);
     free(g_cm);
     free(q_cm);
     free(s_cm);
     free(u_cm);
 
     return info;
 }
 