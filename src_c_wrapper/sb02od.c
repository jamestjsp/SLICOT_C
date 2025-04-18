/**
 * @file sb02od.c
 * @brief C wrapper implementation for SLICOT routine SB02OD
 *
 * This file provides a C wrapper implementation for the SLICOT routine SB02OD,
 * which solves continuous- or discrete-time algebraic Riccati equations
 * using the generalized Schur vectors method.
 */

 #include <stdlib.h>
 #include <ctype.h>
 #include <stddef.h> // For size_t
 #include <string.h> // For memcpy
 
 // Include the header file for this wrapper
 #include "sb02od.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, MIN, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines, slicot_copy_symmetric_part
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions
 
 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * Input matrices are const. Output matrices are X, S, T, U.
  * Output scalars are RCOND, ALFAR, ALFAI, BETA.
  */
 extern void F77_FUNC(sb02od, SB02OD)(
     const char* dico,       // CHARACTER*1 DICO
     const char* jobb,       // CHARACTER*1 JOBB
     const char* fact,       // CHARACTER*1 FACT
     const char* uplo,       // CHARACTER*1 UPLO
     const char* jobl,       // CHARACTER*1 JOBL
     const char* sort,       // CHARACTER*1 SORT
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     const int* p,           // INTEGER P
     const double* a,        // DOUBLE PRECISION A(LDA,*)
     const int* lda,         // INTEGER LDA
     const double* b,        // DOUBLE PRECISION B(LDB,*) (or G)
     const int* ldb,         // INTEGER LDB
     const double* q,        // DOUBLE PRECISION Q(LDQ,*) (or C)
     const int* ldq,         // INTEGER LDQ
     const double* r,        // DOUBLE PRECISION R(LDR,*) (or D)
     const int* ldr,         // INTEGER LDR
     const double* l,        // DOUBLE PRECISION L(LDL,*)
     const int* ldl,         // INTEGER LDL
     double* rcond,          // DOUBLE PRECISION RCOND (output)
     double* x,              // DOUBLE PRECISION X(LDX,*) (output)
     const int* ldx,         // INTEGER LDX
     double* alfar,          // DOUBLE PRECISION ALFAR(*) (output)
     double* alfai,          // DOUBLE PRECISION ALFAI(*) (output)
     double* beta,           // DOUBLE PRECISION BETA(*) (output)
     double* s,              // DOUBLE PRECISION S(LDS,*) (output)
     const int* lds,         // INTEGER LDS
     double* t,              // DOUBLE PRECISION T(LDT,*) (output)
     const int* ldt,         // INTEGER LDT
     double* u,              // DOUBLE PRECISION U(LDU,*) (output)
     const int* ldu,         // INTEGER LDU
     const double* tol,      // DOUBLE PRECISION TOL
     int* iwork,             // INTEGER IWORK(*)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* bwork,             // LOGICAL BWORK(*) -> int*
     int* info,              // INTEGER INFO (output)
     int dico_len,           // Hidden length
     int jobb_len,           // Hidden length
     int fact_len,           // Hidden length
     int uplo_len,           // Hidden length
     int jobl_len,           // Hidden length
     int sort_len            // Hidden length
 );
 
 
 /* C wrapper function definition */
 SLICOT_C_WRAPPER_API
 int slicot_sb02od(char dico, char jobb, char fact, char uplo, char jobl, char sort,
                   int n, int m, int p,
                   const double* a, int lda, const double* b, int ldb,
                   const double* q, int ldq, const double* r, int ldr,
                   const double* l, int ldl, double* rcond,
                   double* x, int ldx, double* alfar, double* alfai, double* beta,
                   double* s, int lds, double* t, int ldt,
                   double* u, int ldu, double tol, int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = -1; /* Use -1 for workspace query */
     double dwork_query;
     double* dwork = NULL;
     int* iwork = NULL;
     int* bwork = NULL; // Map LOGICAL to int
     int iwork_size = 0;
     int bwork_size = 0;
 
     const int dico_len = 1, jobb_len = 1, fact_len = 1, uplo_len = 1;
     const int jobl_len = 1, sort_len = 1;
 
     char dico_upper = toupper(dico);
     char jobb_upper = toupper(jobb);
     char fact_upper = toupper(fact);
     char uplo_upper = toupper(uplo);
     char jobl_upper = toupper(jobl);
     char sort_upper = toupper(sort);
 
     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *b_cm = NULL, *q_cm = NULL, *r_cm = NULL, *l_cm = NULL;
     double *x_cm = NULL, *s_cm = NULL, *t_cm = NULL, *u_cm = NULL;
 
     /* --- Input Parameter Validation --- */
 
     if (n < 0) { info = -7; goto cleanup; }
     if (m < 0) { info = -8; goto cleanup; }
     if (p < 0) { info = -9; goto cleanup; }
     if (dico_upper != 'C' && dico_upper != 'D') { info = -1; goto cleanup; }
     if (jobb_upper != 'B' && jobb_upper != 'G') { info = -2; goto cleanup; }
     if (fact_upper != 'N' && fact_upper != 'C' && fact_upper != 'D' && fact_upper != 'B') { info = -3; goto cleanup; }
     if (uplo_upper != 'U' && uplo_upper != 'L') { info = -4; goto cleanup; }
     if (jobl_upper != 'Z' && jobl_upper != 'N') { info = -5; goto cleanup; }
     if (sort_upper != 'S' && sort_upper != 'U') { info = -6; goto cleanup; }
     // TOL check done by Fortran
 
     // Determine effective dimensions based on flags for validation
     int eff_ldq_f = 1, eff_ldr_f = 1, eff_ldl_f = 1;
     int eff_ldq_rm_cols = 1, eff_ldr_rm_cols = 1, eff_ldl_rm_cols = 1;
 
     int min_lda_f = MAX(1, n);
     int min_ldb_f = MAX(1, n); // For B or G
     int min_ldx_f = MAX(1, n);
     int min_lds_f = (jobb_upper == 'B') ? MAX(1, 2 * n + m) : MAX(1, 2 * n);
     int min_ldt_f = (jobb_upper == 'B') ? MAX(1, 2 * n + m) : ((dico_upper == 'D') ? MAX(1, 2 * n) : 1);
     int min_ldu_f = MAX(1, 2 * n);
 
     if (fact_upper == 'N' || fact_upper == 'D') { // Q is N x N
         eff_ldq_f = MAX(1, n); eff_ldq_rm_cols = n;
     } else { // Q contains C (P x N)
         eff_ldq_f = MAX(1, p); eff_ldq_rm_cols = n;
     }
 
     if (jobb_upper == 'B') {
         if (fact_upper == 'N' || fact_upper == 'C') { // R is M x M
             eff_ldr_f = MAX(1, m); eff_ldr_rm_cols = m;
         } else { // R contains D (P x M)
             eff_ldr_f = MAX(1, p); eff_ldr_rm_cols = m;
         }
         if (jobl_upper == 'N') { // L is N x M
             eff_ldl_f = MAX(1, n); eff_ldl_rm_cols = m;
         }
     } // else R, L not used
 
     // Check leading dimensions based on storage order
     if (row_major) {
         int min_lda_rm_cols = n;
         int min_ldb_rm_cols = (jobb_upper == 'B') ? m : n; // B or G
         int min_ldx_rm_cols = n;
         int min_lds_rm_cols = (jobb_upper == 'B') ? (2 * n + m) : (2 * n);
         int min_ldt_rm_cols = (dico_upper == 'C' && jobb_upper == 'G') ? 1 : ((jobb_upper == 'B') ? (2 * n + m) : (2 * n));
         int min_ldu_rm_cols = 2 * n;
 
         if (lda < min_lda_rm_cols) { info = -11; goto cleanup; }
         if (ldb < min_ldb_rm_cols) { info = -13; goto cleanup; }
         if (ldq < eff_ldq_rm_cols) { info = -15; goto cleanup; }
         if (jobb_upper == 'B' && ldr < eff_ldr_rm_cols) { info = -17; goto cleanup; }
         if (jobb_upper == 'B' && jobl_upper == 'N' && ldl < eff_ldl_rm_cols) { info = -19; goto cleanup; }
         if (ldx < min_ldx_rm_cols) { info = -21; goto cleanup; }
         if (lds < min_lds_rm_cols) { info = -26; goto cleanup; }
         if (ldt < min_ldt_rm_cols) { info = -28; goto cleanup; }
         if (ldu < min_ldu_rm_cols) { info = -30; goto cleanup; }
     } else {
         if (lda < min_lda_f) { info = -11; goto cleanup; }
         if (ldb < min_ldb_f) { info = -13; goto cleanup; }
         if (ldq < eff_ldq_f) { info = -15; goto cleanup; }
         if (jobb_upper == 'B' && ldr < eff_ldr_f) { info = -17; goto cleanup; }
         if (jobb_upper == 'B' && jobl_upper == 'N' && ldl < eff_ldl_f) { info = -19; goto cleanup; }
         if (ldx < min_ldx_f) { info = -21; goto cleanup; }
         if (lds < min_lds_f) { info = -26; goto cleanup; }
         if (ldt < min_ldt_f) { info = -28; goto cleanup; }
         if (ldu < min_ldu_f) { info = -30; goto cleanup; }
     }
 
     /* --- Workspace Allocation --- */
 
     // Allocate IWORK and BWORK
     iwork_size = (jobb_upper == 'B') ? MAX(1, MAX(m, 2 * n)) : MAX(1, 2 * n);
     bwork_size = MAX(1, 2 * n);
     iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
     CHECK_ALLOC(iwork);
     bwork = (int*)malloc((size_t)bwork_size * sizeof(int)); // Use int for LOGICAL
     CHECK_ALLOC(bwork);
 
     // Allocate DWORK based on query
     ldwork = -1; // Query mode
     F77_FUNC(sb02od, SB02OD)(&dico_upper, &jobb_upper, &fact_upper, &uplo_upper, &jobl_upper, &sort_upper,
                              &n, &m, &p, a, &lda, b, &ldb, q, &ldq, r, &ldr, l, &ldl,
                              rcond, x, &ldx, alfar, alfai, beta, s, &lds, t, &ldt, u, &ldu,
                              &tol, iwork, &dwork_query, &ldwork, bwork, &info,
                              dico_len, jobb_len, fact_len, uplo_len, jobl_len, sort_len);
 
     if (info < 0) { goto cleanup; } // Query failed due to invalid argument
     info = 0; // Reset info after query
 
     // Get the required dwork size from query result
     ldwork = (int)dwork_query;
     // Check against minimum documented size
     int min_ldwork = 3;
     if (jobb_upper == 'G') {
         min_ldwork = (dico_upper == 'C') ? MAX(3, 6 * n) : MAX(7 * (2 * n + 1) + 16, 16 * n);
     } else { // JOBB = 'B'
         min_ldwork = MAX(7 * (2 * n + 1) + 16, MAX(16 * n, MAX(2 * n + m, 3 * m)));
     }
     ldwork = MAX(ldwork, min_ldwork);
 
 
     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure
 
     /* --- Prepare Arrays and Call Fortran Routine --- */
     size_t elem_size = sizeof(double);
 
     // Determine sizes for potential copies (moved before if/else)
     size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
     size_t b_rows = n; size_t b_cols = (jobb_upper == 'B') ? m : n; size_t b_size = b_rows * b_cols; // B or G
     size_t q_rows = (fact_upper == 'N' || fact_upper == 'D') ? n : p; size_t q_cols = n; size_t q_size = q_rows * q_cols; // Q or C
     size_t r_rows = (fact_upper == 'N' || fact_upper == 'C') ? m : p; size_t r_cols = m; size_t r_size = r_rows * r_cols; // R or D
     size_t l_rows = n; size_t l_cols = m; size_t l_size = l_rows * l_cols;
     size_t x_rows = n; size_t x_cols = n; size_t x_size = x_rows * x_cols;
     int s_t_rows = (jobb_upper == 'B') ? (2*n+m) : (2*n); // Rows for S and T
     int s_cols = s_t_rows; size_t s_size = (size_t)s_t_rows * s_cols;
     int t_cols = 2*n; size_t t_size = (size_t)s_t_rows * t_cols; // T is S_rows x 2N
     size_t u_rows = 2*n; size_t u_cols = 2*n; size_t u_size = u_rows * u_cols;
 
     if (row_major) {
         /* --- Row-Major Case --- */
 
         /* Allocate memory for column-major copies */
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (q_size > 0) { q_cm = (double*)malloc(q_size * elem_size); CHECK_ALLOC(q_cm); }
         if (jobb_upper == 'B' && r_size > 0) { r_cm = (double*)malloc(r_size * elem_size); CHECK_ALLOC(r_cm); }
         if (jobb_upper == 'B' && jobl_upper == 'N' && l_size > 0) { l_cm = (double*)malloc(l_size * elem_size); CHECK_ALLOC(l_cm); }
         if (x_size > 0) { x_cm = (double*)malloc(x_size * elem_size); CHECK_ALLOC(x_cm); }
         if (s_size > 0) { s_cm = (double*)malloc(s_size * elem_size); CHECK_ALLOC(s_cm); }
         if (t_size > 0 && !(dico_upper == 'C' && jobb_upper == 'G')) { t_cm = (double*)malloc(t_size * elem_size); CHECK_ALLOC(t_cm); }
         if (u_size > 0) { u_cm = (double*)malloc(u_size * elem_size); CHECK_ALLOC(u_cm); }
 
         /* Transpose C inputs to Fortran copies */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         if (jobb_upper == 'G' && b_size > 0) slicot_transpose_symmetric_to_fortran(b, b_cm, n, uplo_upper, elem_size);
         else if (jobb_upper == 'B' && b_size > 0) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, elem_size);
         if ((fact_upper == 'N' || fact_upper == 'D') && q_size > 0) slicot_transpose_symmetric_to_fortran(q, q_cm, n, uplo_upper, elem_size);
         else if ((fact_upper == 'C' || fact_upper == 'B') && q_size > 0) slicot_transpose_to_fortran(q, q_cm, q_rows, q_cols, elem_size);
         if (jobb_upper == 'B') {
             if ((fact_upper == 'N' || fact_upper == 'C') && r_size > 0) slicot_transpose_symmetric_to_fortran(r, r_cm, m, uplo_upper, elem_size);
             else if ((fact_upper == 'D' || fact_upper == 'B') && r_size > 0) slicot_transpose_to_fortran(r, r_cm, r_rows, r_cols, elem_size);
             if (jobl_upper == 'N' && l_size > 0) slicot_transpose_to_fortran(l, l_cm, l_rows, l_cols, elem_size);
         }
 
         /* Fortran leading dimensions */
         int lda_f = (a_rows > 0) ? a_rows : 1;
         int ldb_f = (b_rows > 0) ? b_rows : 1;
         int ldq_f = (q_rows > 0) ? q_rows : 1;
         int ldr_f = (jobb_upper == 'B' && r_rows > 0) ? r_rows : 1;
         int ldl_f = (jobb_upper == 'B' && jobl_upper == 'N' && l_rows > 0) ? l_rows : 1;
         int ldx_f = (x_rows > 0) ? x_rows : 1;
         int lds_f = (s_t_rows > 0) ? s_t_rows : 1; // Use calculated rows for S/T
         int ldt_f = (s_t_rows > 0) ? s_t_rows : 1;
         int ldu_f = (u_rows > 0) ? u_rows : 1;
 
 
         /* Call the Fortran routine */
         F77_FUNC(sb02od, SB02OD)(&dico_upper, &jobb_upper, &fact_upper, &uplo_upper, &jobl_upper, &sort_upper,
                                  &n, &m, &p,
                                  a_cm, &lda_f, b_cm, &ldb_f, q_cm, &ldq_f,
                                  (jobb_upper == 'B' ? r_cm : NULL), &ldr_f,
                                  (jobb_upper == 'B' && jobl_upper == 'N' ? l_cm : NULL), &ldl_f,
                                  rcond,
                                  x_cm, &ldx_f, alfar, alfai, beta,
                                  s_cm, &lds_f,
                                  (!(dico_upper == 'C' && jobb_upper == 'G') ? t_cm : NULL), &ldt_f,
                                  u_cm, &ldu_f,
                                  &tol, iwork, dwork, &ldwork, bwork, &info,
                                  dico_len, jobb_len, fact_len, uplo_len, jobl_len, sort_len);
 
         /* Copy back results from column-major temps to original row-major arrays */
         if (info == 0 || info == 6) { // Copy back even if INFO=6
             if (x_size > 0) slicot_transpose_symmetric_to_c(x_cm, x, n, uplo_upper, elem_size); // X is symmetric
             if (s_size > 0) slicot_transpose_to_c(s_cm, s, s_t_rows, s_cols, elem_size); // Use s_t_rows
             if (t_size > 0 && t_cm) slicot_transpose_to_c(t_cm, t, s_t_rows, t_cols, elem_size); // Use s_t_rows
             if (u_size > 0) slicot_transpose_to_c(u_cm, u, u_rows, u_cols, elem_size);
             // RCOND, ALFAR, ALFAI, BETA modified directly
         }
         /* Temps freed in cleanup */
 
     } else {
         /* --- Column-Major Case --- */
         // Need to copy symmetric inputs G, Q, R to ensure full matrix is passed if needed
         if (jobb_upper == 'G' && b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); slicot_copy_symmetric_part(b, b_cm, n, uplo_upper, ldb, elem_size); }
         if ((fact_upper == 'N' || fact_upper == 'D') && q_size > 0) { q_cm = (double*)malloc(q_size * elem_size); CHECK_ALLOC(q_cm); slicot_copy_symmetric_part(q, q_cm, n, uplo_upper, ldq, elem_size); }
         if (jobb_upper == 'B' && (fact_upper == 'N' || fact_upper == 'C') && r_size > 0) { r_cm = (double*)malloc(r_size * elem_size); CHECK_ALLOC(r_cm); slicot_copy_symmetric_part(r, r_cm, m, uplo_upper, ldr, elem_size); }
 
 
         /* Call the Fortran routine directly with user-provided arrays (or copies) */
         F77_FUNC(sb02od, SB02OD)(&dico_upper, &jobb_upper, &fact_upper, &uplo_upper, &jobl_upper, &sort_upper,
                                  &n, &m, &p,
                                  a, &lda,
                                  (jobb_upper == 'G' ? (b_cm ? b_cm : b) : b), &ldb, // Use copy of G if made
                                  ((fact_upper == 'N' || fact_upper == 'D') ? (q_cm ? q_cm : q) : q), &ldq, // Use copy of Q if made
                                  (jobb_upper == 'B' ? ((fact_upper == 'N' || fact_upper == 'C') ? (r_cm ? r_cm : r) : r) : NULL), &ldr, // Use copy of R if made
                                  (jobb_upper == 'B' && jobl_upper == 'N' ? l : NULL), &ldl,
                                  rcond, x, &ldx, alfar, alfai, beta, s, &lds, t, &ldt, u, &ldu,
                                  &tol, iwork, dwork, &ldwork, bwork, &info,
                                  dico_len, jobb_len, fact_len, uplo_len, jobl_len, sort_len);
         // Output arrays X, S, T, U and scalars modified directly
         // Copy back symmetric solution X from temp if needed
         if (info == 0 || info == 6) {
              // Check if q_cm was allocated before using it
              if (q_cm && (fact_upper == 'N' || fact_upper == 'D') && q_size > 0) {
                  slicot_copy_symmetric_part(q_cm, x, n, uplo_upper, ldx, elem_size); // Copy X back to original X
              }
         }
     }
 
 cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(iwork);
     free(bwork);
     free(a_cm);
     free(b_cm); // Might be G copy
     free(q_cm); // Might be Q copy or C
     free(r_cm); // Might be R copy
     free(l_cm);
     free(x_cm);
     free(s_cm);
     free(t_cm);
     free(u_cm);
 
     return info;
 }
 