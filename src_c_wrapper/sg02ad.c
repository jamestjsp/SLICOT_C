/**
 * @file sg02ad.c
 * @brief C wrapper implementation for SLICOT routine SG02AD
 *
 * This file provides a C wrapper implementation for the SLICOT routine SG02AD,
 * which solves continuous- or discrete-time algebraic Riccati equations
 * for descriptor systems using the generalized Schur vectors method.
 */

 #include <stdlib.h>
 #include <ctype.h>
 #include <stddef.h> // For size_t
 #include <string.h> // For memcpy
 #include <stdio.h> // For toupper if not in ctype.h

 // Include the header file for this wrapper
 #include "sg02ad.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, MIN, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines, slicot_copy_symmetric_part
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * Input matrices are const. Output matrices X, S, T, U are non-const.
  * Output scalars RCONDU, ALFAR, ALFAI, BETA, IWARN are pointers.
  */
 extern void F77_FUNC(sg02ad, SG02AD)(
     const char* dico,       // CHARACTER*1 DICO
     const char* jobb,       // CHARACTER*1 JOBB
     const char* fact,       // CHARACTER*1 FACT
     const char* uplo,       // CHARACTER*1 UPLO
     const char* jobl,       // CHARACTER*1 JOBL
     const char* scal,       // CHARACTER*1 SCAL
     const char* sort,       // CHARACTER*1 SORT
     const char* acc,        // CHARACTER*1 ACC
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     const int* p,           // INTEGER P
     const double* a,        // DOUBLE PRECISION A(LDA,*)
     const int* lda,         // INTEGER LDA
     const double* e,        // DOUBLE PRECISION E(LDE,*)
     const int* lde,         // INTEGER LDE
     const double* b,        // DOUBLE PRECISION B(LDB,*) (or G)
     const int* ldb,         // INTEGER LDB
     const double* q,        // DOUBLE PRECISION Q(LDQ,*) (or C)
     const int* ldq,         // INTEGER LDQ
     const double* r,        // DOUBLE PRECISION R(LDR,*) (or D)
     const int* ldr,         // INTEGER LDR
     const double* l,        // DOUBLE PRECISION L(LDL,*)
     const int* ldl,         // INTEGER LDL
     double* rcondu,         // DOUBLE PRECISION RCONDU (output)
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
     int* iwarn,             // INTEGER IWARN (output)
     int* info,              // INTEGER INFO (output)
     int dico_len,           // Hidden length
     int jobb_len,           // Hidden length
     int fact_len,           // Hidden length
     int uplo_len,           // Hidden length
     int jobl_len,           // Hidden length
     int scal_len,           // Hidden length
     int sort_len,           // Hidden length
     int acc_len             // Hidden length
 );


 /* C wrapper function definition */
 int slicot_sg02ad(char dico, char jobb, char fact, char uplo, char jobl, char scal, char sort, char acc,
                   int n, int m, int p,
                   const double* a, int lda, const double* e, int lde,
                   const double* b, int ldb, const double* q, int ldq,
                   const double* r, int ldr, const double* l, int ldl,
                   double* rcondu, double* x, int ldx,
                   double* alfar, double* alfai, double* beta,
                   double* s, int lds, double* t, int ldt,
                   double* u, int ldu, double tol, int* iwarn, int row_major)
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
     const int jobl_len = 1, scal_len = 1, sort_len = 1, acc_len = 1;

     char dico_upper = toupper(dico);
     char jobb_upper = toupper(jobb);
     char fact_upper = toupper(fact);
     char uplo_upper = toupper(uplo);
     char jobl_upper = toupper(jobl);
     char scal_upper = toupper(scal);
     char sort_upper = toupper(sort);
     char acc_upper = toupper(acc);

     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *e_cm = NULL, *b_cm = NULL, *q_cm = NULL, *r_cm = NULL, *l_cm = NULL;
     double *x_cm = NULL, *s_cm = NULL, *t_cm = NULL, *u_cm = NULL;

     /* --- Input Parameter Validation --- */

     if (n < 0) { info = -9; goto cleanup; }
     if (m < 0) { info = -10; goto cleanup; }
     if (p < 0) { info = -11; goto cleanup; }
     if (dico_upper != 'C' && dico_upper != 'D') { info = -1; goto cleanup; }
     if (jobb_upper != 'B' && jobb_upper != 'G') { info = -2; goto cleanup; }
     if (fact_upper != 'N' && fact_upper != 'C' && fact_upper != 'D' && fact_upper != 'B') { info = -3; goto cleanup; }
     if (uplo_upper != 'U' && uplo_upper != 'L') { info = -4; goto cleanup; }
     if (jobl_upper != 'Z' && jobl_upper != 'N') { info = -5; goto cleanup; }
     if (scal_upper != 'G' && scal_upper != 'N') { info = -6; goto cleanup; }
     if (sort_upper != 'S' && sort_upper != 'U') { info = -7; goto cleanup; }
     if (acc_upper != 'R' && acc_upper != 'N') { info = -8; goto cleanup; }
     // TOL check done by Fortran

     // Determine effective dimensions based on flags for validation
     int eff_ldq_f = 1, eff_ldr_f = 1, eff_ldl_f = 1;
     int eff_ldq_rm_cols = 1, eff_ldr_rm_cols = 1, eff_ldl_rm_cols = 1;

     int min_lda_f = MAX(1, n);
     int min_lde_f = MAX(1, n);
     int min_ldb_f = MAX(1, n); // For B or G
     int min_ldx_f = MAX(1, n);
     int min_lds_f = (jobb_upper == 'B') ? MAX(1, 2 * n + m) : MAX(1, 2 * n);
     int min_ldt_f = (jobb_upper == 'B') ? MAX(1, 2 * n + m) : MAX(1, 2 * n); // T same size as S
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
         int min_lde_rm_cols = n;
         int min_ldb_rm_cols = (jobb_upper == 'B') ? m : n; // B or G
         int min_ldx_rm_cols = n;
         int min_lds_rm_cols = (jobb_upper == 'B') ? (2 * n + m) : (2 * n);
         int min_ldt_rm_cols = (jobb_upper == 'B') ? (2 * n + m) : (2 * n); // T same size as S
         int min_ldu_rm_cols = 2 * n;

         if (lda < min_lda_rm_cols) { info = -13; goto cleanup; }
         if (lde < min_lde_rm_cols) { info = -15; goto cleanup; }
         if (ldb < min_ldb_rm_cols) { info = -17; goto cleanup; }
         if (ldq < eff_ldq_rm_cols) { info = -19; goto cleanup; }
         if (jobb_upper == 'B' && ldr < eff_ldr_rm_cols) { info = -21; goto cleanup; }
         if (jobb_upper == 'B' && jobl_upper == 'N' && ldl < eff_ldl_rm_cols) { info = -23; goto cleanup; }
         if (ldx < min_ldx_rm_cols) { info = -25; goto cleanup; }
         if (lds < min_lds_rm_cols) { info = -30; goto cleanup; }
         if (ldt < min_ldt_rm_cols) { info = -32; goto cleanup; }
         if (ldu < min_ldu_rm_cols) { info = -34; goto cleanup; }
     } else {
         if (lda < min_lda_f) { info = -13; goto cleanup; }
         if (lde < min_lde_f) { info = -15; goto cleanup; }
         if (ldb < min_ldb_f) { info = -17; goto cleanup; }
         if (ldq < eff_ldq_f) { info = -19; goto cleanup; }
         if (jobb_upper == 'B' && ldr < eff_ldr_f) { info = -21; goto cleanup; }
         if (jobb_upper == 'B' && jobl_upper == 'N' && ldl < eff_ldl_f) { info = -23; goto cleanup; }
         if (ldx < min_ldx_f) { info = -25; goto cleanup; }
         if (lds < min_lds_f) { info = -30; goto cleanup; }
         if (ldt < min_ldt_f) { info = -32; goto cleanup; }
         if (ldu < min_ldu_f) { info = -34; goto cleanup; }
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
     F77_FUNC(sg02ad, SG02AD)(&dico_upper, &jobb_upper, &fact_upper, &uplo_upper, &jobl_upper, &scal_upper, &sort_upper, &acc_upper,
                              &n, &m, &p, a, &lda, e, &lde, b, &ldb, q, &ldq, r, &ldr, l, &ldl,
                              rcondu, x, &ldx, alfar, alfai, beta, s, &lds, t, &ldt, u, &ldu,
                              &tol, iwork, &dwork_query, &ldwork, bwork, iwarn, &info,
                              dico_len, jobb_len, fact_len, uplo_len, jobl_len, scal_len, sort_len, acc_len);

     if (info < 0 && info != -39) { info = info; goto cleanup; } // Query failed due to invalid argument
     info = 0; // Reset info after query

     // Get the required dwork size from query result
     ldwork = (int)dwork_query;
     // Check against minimum documented size
     int min_ldwork = 1;
     if (jobb_upper == 'G') {
         min_ldwork = MAX(7 * (2 * n + 1) + 16, 16 * n);
     } else { // JOBB = 'B'
         min_ldwork = MAX(7 * (2 * n + 1) + 16, MAX(16 * n, MAX(2 * n + m, 3 * m)));
     }
     ldwork = MAX(ldwork, min_ldwork);


     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure

     /* --- Prepare Arrays and Call Fortran Routine --- */
     size_t elem_size = sizeof(double);

     // Determine sizes for potential copies
     size_t a_size = (size_t)n * n;
     size_t e_size = (size_t)n * n;
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
         if (e_size > 0) { e_cm = (double*)malloc(e_size * elem_size); CHECK_ALLOC(e_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (q_size > 0) { q_cm = (double*)malloc(q_size * elem_size); CHECK_ALLOC(q_cm); }
         if (jobb_upper == 'B' && r_size > 0) { r_cm = (double*)malloc(r_size * elem_size); CHECK_ALLOC(r_cm); }
         if (jobb_upper == 'B' && jobl_upper == 'N' && l_size > 0) { l_cm = (double*)malloc(l_size * elem_size); CHECK_ALLOC(l_cm); }
         if (x_size > 0) { x_cm = (double*)malloc(x_size * elem_size); CHECK_ALLOC(x_cm); }
         if (s_size > 0) { s_cm = (double*)malloc(s_size * elem_size); CHECK_ALLOC(s_cm); }
         if (t_size > 0) { t_cm = (double*)malloc(t_size * elem_size); CHECK_ALLOC(t_cm); }
         if (u_size > 0) { u_cm = (double*)malloc(u_size * elem_size); CHECK_ALLOC(u_cm); }

         /* --- FIX: Define row/col variables needed for transposition --- */
         int a_rows = n; int a_cols = n;
         int e_rows = n; int e_cols = n;
         /* Note: b_rows, b_cols, q_rows, q_cols etc. are already defined above */
         /* -------------------------------------------------------------- */

         /* Transpose C inputs to Fortran copies */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size); // Line 262 (now uses defined vars)
         if (e_size > 0) slicot_transpose_to_fortran(e, e_cm, e_rows, e_cols, elem_size); // Line 263 (now uses defined vars)
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
         int lda_f = (a_rows > 0) ? a_rows : 1; // Line 275 (now uses defined vars)
         int lde_f = (e_rows > 0) ? e_rows : 1; // Line 276 (now uses defined vars)
         int ldb_f = (b_rows > 0) ? b_rows : 1;
         int ldq_f = (q_rows > 0) ? q_rows : 1;
         int ldr_f = (jobb_upper == 'B' && r_rows > 0) ? r_rows : 1;
         int ldl_f = (jobb_upper == 'B' && jobl_upper == 'N' && l_rows > 0) ? l_rows : 1;
         int ldx_f = (x_rows > 0) ? x_rows : 1;
         int lds_f = (s_t_rows > 0) ? s_t_rows : 1;
         int ldt_f = (s_t_rows > 0) ? s_t_rows : 1;
         int ldu_f = (u_rows > 0) ? u_rows : 1;


         /* Call the Fortran routine */
         F77_FUNC(sg02ad, SG02AD)(&dico_upper, &jobb_upper, &fact_upper, &uplo_upper, &jobl_upper, &scal_upper, &sort_upper, &acc_upper,
                                  &n, &m, &p,
                                  a_cm, &lda_f, e_cm, &lde_f, b_cm, &ldb_f, q_cm, &ldq_f,
                                  (jobb_upper == 'B' ? r_cm : NULL), &ldr_f,
                                  (jobb_upper == 'B' && jobl_upper == 'N' ? l_cm : NULL), &ldl_f,
                                  rcondu, x_cm, &ldx_f, alfar, alfai, beta,
                                  s_cm, &lds_f, t_cm, &ldt_f, u_cm, &ldu_f,
                                  &tol, iwork, dwork, &ldwork, bwork, iwarn, &info,
                                  dico_len, jobb_len, fact_len, uplo_len, jobl_len, scal_len, sort_len, acc_len);

         /* Copy back results from column-major temps to original row-major arrays */
         if (info == 0 || info == 7) { // Copy back even if INFO=7
             if (x_size > 0) slicot_transpose_symmetric_to_c(x_cm, x, n, uplo_upper, elem_size); // X is symmetric
             if (s_size > 0) slicot_transpose_to_c(s_cm, s, s_t_rows, s_cols, elem_size);
             if (t_size > 0) slicot_transpose_to_c(t_cm, t, s_t_rows, t_cols, elem_size);
             if (u_size > 0) slicot_transpose_to_c(u_cm, u, u_rows, u_cols, elem_size);
             // RCONDU, ALFAR, ALFAI, BETA, IWARN modified directly
         }
         /* Temps freed in cleanup */

     } else {
         /* --- Column-Major Case --- */
         // Need to copy symmetric inputs G, Q, R to ensure full matrix is passed if needed
         if (jobb_upper == 'G' && b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); slicot_copy_symmetric_part(b, b_cm, n, uplo_upper, ldb, elem_size); }
         if ((fact_upper == 'N' || fact_upper == 'D') && q_size > 0) { q_cm = (double*)malloc(q_size * elem_size); CHECK_ALLOC(q_cm); slicot_copy_symmetric_part(q, q_cm, n, uplo_upper, ldq, elem_size); }
         if (jobb_upper == 'B' && (fact_upper == 'N' || fact_upper == 'C') && r_size > 0) { r_cm = (double*)malloc(r_size * elem_size); CHECK_ALLOC(r_cm); slicot_copy_symmetric_part(r, r_cm, m, uplo_upper, ldr, elem_size); }


         /* Call the Fortran routine directly with user-provided arrays (or copies) */
         F77_FUNC(sg02ad, SG02AD)(&dico_upper, &jobb_upper, &fact_upper, &uplo_upper, &jobl_upper, &scal_upper, &sort_upper, &acc_upper,
                                  &n, &m, &p,
                                  a, &lda, e, &lde,
                                  (jobb_upper == 'G' ? (b_cm ? b_cm : b) : b), &ldb, // Use copy of G if made
                                  ((fact_upper == 'N' || fact_upper == 'D') ? (q_cm ? q_cm : q) : q), &ldq, // Use copy of Q if made
                                  (jobb_upper == 'B' ? ((fact_upper == 'N' || fact_upper == 'C') ? (r_cm ? r_cm : r) : r) : NULL), &ldr, // Use copy of R if made
                                  (jobb_upper == 'B' && jobl_upper == 'N' ? l : NULL), &ldl,
                                  rcondu, x, &ldx, alfar, alfai, beta, s, &lds, t, &ldt, u, &ldu,
                                  &tol, iwork, dwork, &ldwork, bwork, iwarn, &info,
                                  dico_len, jobb_len, fact_len, uplo_len, jobl_len, scal_len, sort_len, acc_len);
         // Output arrays X, S, T, U and scalars modified directly
         // Copy back symmetric solution X from temp if needed
         if (info == 0 || info == 7) {
              // Check if q_cm was allocated before using it
              // FIX: The logic here was attempting to copy q_cm to x, which is incorrect.
              // The solution X is already computed directly into the 'x' array in the column-major case.
              // If X needs to be made fully symmetric (if only one triangle was computed),
              // that logic would need to be added here, but the Fortran routine likely returns
              // the full symmetric X already. We remove the incorrect copy.
              // if (q_cm && (fact_upper == 'N' || fact_upper == 'D') && q_size > 0) {
              //     slicot_copy_symmetric_part(q_cm, x, n, uplo_upper, ldx, elem_size); // Incorrect: This copies Q to X
              // }
              // If X needs explicit symmetrization (unlikely for this routine, but for illustration):
              // slicot_copy_symmetric_part(x, x, n, uplo_upper, ldx, elem_size); // Copies relevant triangle of X to fill the other triangle
         }
     }

 cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(iwork);
     free(bwork);
     free(a_cm); free(e_cm); free(b_cm); free(q_cm); free(r_cm); free(l_cm);
     free(x_cm); free(s_cm); free(t_cm); free(u_cm);

     return info;
 }
