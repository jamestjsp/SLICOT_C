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
 #include "slicot_utils.h" // Assumed to contain MAX, MIN, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines, slicot_copy_symmetric_part, slicot_transpose_symmetric_to_fortran, slicot_transpose_symmetric_to_c
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
     const double* a,        // DOUBLE PRECISION A(LDA,*) (in)
     const int* lda,         // INTEGER LDA
     const double* b,        // DOUBLE PRECISION B(LDB,*) (in, B or G)
     const int* ldb,         // INTEGER LDB
     const double* q,        // DOUBLE PRECISION Q(LDQ,*) (in, Q or C)
     const int* ldq,         // INTEGER LDQ
     const double* r,        // DOUBLE PRECISION R(LDR,*) (in, R or D)
     const int* ldr,         // INTEGER LDR
     const double* l,        // DOUBLE PRECISION L(LDL,*) (in)
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
     const double *a_ptr, *b_ptr, *q_ptr, *r_ptr, *l_ptr; // Input pointers
     double *x_ptr, *s_ptr, *t_ptr, *u_ptr; // Output pointers
     int lda_f, ldb_f, ldq_f, ldr_f, ldl_f, ldx_f, lds_f, ldt_f, ldu_f;


     /* --- Input Parameter Validation --- */

     if (dico_upper != 'C' && dico_upper != 'D') { info = -1; goto cleanup; }
     if (jobb_upper != 'B' && jobb_upper != 'G') { info = -2; goto cleanup; }
     if (fact_upper != 'N' && fact_upper != 'C' && fact_upper != 'D' && fact_upper != 'B') { info = -3; goto cleanup; }
     if (uplo_upper != 'U' && uplo_upper != 'L') { info = -4; goto cleanup; }
     if (jobl_upper != 'Z' && jobl_upper != 'N') { info = -5; goto cleanup; }
     if (sort_upper != 'S' && sort_upper != 'U') { info = -6; goto cleanup; }
     if (n < 0) { info = -7; goto cleanup; }
     if (m < 0) { info = -8; goto cleanup; }
     if (p < 0) { info = -9; goto cleanup; }
     // TOL check done by Fortran

     // Determine effective dimensions based on flags for validation
     int q_rows_f = (fact_upper == 'N' || fact_upper == 'D') ? n : p;
     int r_rows_f = (jobb_upper == 'B' && (fact_upper == 'N' || fact_upper == 'C')) ? m : (jobb_upper == 'B' ? p : 1);
     int l_rows_f = (jobb_upper == 'B' && jobl_upper == 'N') ? n : 1;

     int min_lda_f = MAX(1, n);
     int min_ldb_f = MAX(1, n); // For B or G
     int min_ldq_f = MAX(1, q_rows_f); // For Q or C
     int min_ldr_f = MAX(1, r_rows_f); // For R or D (only if JOBB='B')
     int min_ldl_f = MAX(1, l_rows_f); // For L (only if JOBB='B', JOBL='N')
     int min_ldx_f = MAX(1, n);
     int s_t_rows_f = (jobb_upper == 'B') ? (2 * n + m) : (2 * n);
     int min_lds_f = MAX(1, s_t_rows_f);
     int min_ldt_f = (dico_upper == 'C' && jobb_upper == 'G') ? 1 : MAX(1, s_t_rows_f); // T not ref'd if DICO='C',JOBB='G'
     int min_ldu_f = MAX(1, 2 * n);

     if (row_major) {
         // Determine effective column counts for row-major LDA checks
         int b_cols_rm = (jobb_upper == 'B') ? m : n; // B(n,m) or G(n,n)
         int q_cols_rm = n; // Q(n,n) or C(p,n)
         int r_cols_rm = (jobb_upper == 'B') ? m : 1; // R(m,m) or D(p,m)
         int l_cols_rm = (jobb_upper == 'B' && jobl_upper == 'N') ? m : 1; // L(n,m)
         int s_cols_rm = s_t_rows_f; // S is square
         int t_cols_rm = (dico_upper == 'C' && jobb_upper == 'G') ? 1 : s_t_rows_f; // T is s_t_rows x s_t_rows (or not ref)
         int u_cols_rm = 2 * n; // U is 2n x 2n

         if (lda < n) { info = -11; goto cleanup; } // A(n,n)
         if (ldb < b_cols_rm) { info = -13; goto cleanup; }
         if (ldq < q_cols_rm) { info = -15; goto cleanup; }
         if (jobb_upper == 'B' && ldr < r_cols_rm) { info = -17; goto cleanup; }
         if (jobb_upper == 'B' && jobl_upper == 'N' && ldl < l_cols_rm) { info = -19; goto cleanup; }
         if (ldx < n) { info = -21; goto cleanup; } // X(n,n)
         if (lds < s_cols_rm) { info = -26; goto cleanup; }
         if (ldt < t_cols_rm) { info = -28; goto cleanup; }
         if (ldu < u_cols_rm) { info = -30; goto cleanup; }

     } else { // Column-major checks
         if (lda < min_lda_f) { info = -11; goto cleanup; }
         if (ldb < min_ldb_f) { info = -13; goto cleanup; }
         if (ldq < min_ldq_f) { info = -15; goto cleanup; }
         if (jobb_upper == 'B' && ldr < min_ldr_f) { info = -17; goto cleanup; }
         if (jobb_upper == 'B' && jobl_upper == 'N' && ldl < min_ldl_f) { info = -19; goto cleanup; }
         if (ldx < min_ldx_f) { info = -21; goto cleanup; }
         if (lds < min_lds_f) { info = -26; goto cleanup; }
         if (ldt < min_ldt_f) { info = -28; goto cleanup; }
         if (ldu < min_ldu_f) { info = -30; goto cleanup; }
     }

     /* --- Prepare arrays for column-major format if using row-major --- */
     size_t elem_size = sizeof(double);

     // Determine sizes for potential copies (moved before if/else)
     size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
     size_t b_rows = n; size_t b_cols = (jobb_upper == 'B') ? m : n; size_t b_size = b_rows * b_cols; // B or G
     size_t q_rows = q_rows_f; size_t q_cols = n; size_t q_size = q_rows * q_cols; // Q or C
     size_t r_rows = r_rows_f; size_t r_cols = m; size_t r_size = r_rows * r_cols; // R or D
     size_t l_rows = n; size_t l_cols = m; size_t l_size = (jobb_upper == 'B' && jobl_upper == 'N') ? l_rows * l_cols : 0; // L
     size_t x_rows = n; size_t x_cols = n; size_t x_size = x_rows * x_cols; // X
     int s_t_rows = s_t_rows_f; // Rows for S and T
     int s_cols = s_t_rows; size_t s_size = (size_t)s_t_rows * s_cols; // S
     int t_cols = s_t_rows; size_t t_size = (size_t)s_t_rows * t_cols; // T square (or not ref'd)
     size_t u_rows = 2*n; size_t u_cols = 2*n; size_t u_size = u_rows * u_cols; // U

     if (row_major) {
         /* Allocate memory for column-major copies */
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (q_size > 0) { q_cm = (double*)malloc(q_size * elem_size); CHECK_ALLOC(q_cm); }
         if (jobb_upper == 'B' && r_size > 0) { r_cm = (double*)malloc(r_size * elem_size); CHECK_ALLOC(r_cm); }
         if (l_size > 0) { l_cm = (double*)malloc(l_size * elem_size); CHECK_ALLOC(l_cm); }
         // Output arrays
         if (x_size > 0) { x_cm = (double*)malloc(x_size * elem_size); CHECK_ALLOC(x_cm); }
         if (s_size > 0) { s_cm = (double*)malloc(s_size * elem_size); CHECK_ALLOC(s_cm); }
         if (t_size > 0 && !(dico_upper == 'C' && jobb_upper == 'G')) { t_cm = (double*)malloc(t_size * elem_size); CHECK_ALLOC(t_cm); } // T not ref'd if DICO='C',JOBB='G'
         if (u_size > 0) { u_cm = (double*)malloc(u_size * elem_size); CHECK_ALLOC(u_cm); }

         /* Transpose C inputs to Fortran copies */
         if (a_cm) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         if (jobb_upper == 'G' && b_cm) slicot_transpose_symmetric_to_fortran(b, b_cm, n, uplo_upper, elem_size); // G is symmetric
         else if (jobb_upper == 'B' && b_cm) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, elem_size); // B is general
         if ((fact_upper == 'N' || fact_upper == 'D') && q_cm) slicot_transpose_symmetric_to_fortran(q, q_cm, n, uplo_upper, elem_size); // Q is symmetric
         else if ((fact_upper == 'C' || fact_upper == 'B') && q_cm) slicot_transpose_to_fortran(q, q_cm, q_rows, q_cols, elem_size); // C is general
         if (jobb_upper == 'B') {
             if ((fact_upper == 'N' || fact_upper == 'C') && r_cm) slicot_transpose_symmetric_to_fortran(r, r_cm, m, uplo_upper, elem_size); // R is symmetric
             else if ((fact_upper == 'D' || fact_upper == 'B') && r_cm) slicot_transpose_to_fortran(r, r_cm, r_rows, r_cols, elem_size); // D is general
             if (l_cm) slicot_transpose_to_fortran(l, l_cm, l_rows, l_cols, elem_size); // L is general
         }

         /* Fortran leading dimensions */
         lda_f = (a_rows > 0) ? a_rows : 1;
         ldb_f = (b_rows > 0) ? b_rows : 1;
         ldq_f = (q_rows > 0) ? q_rows : 1;
         ldr_f = (jobb_upper == 'B' && r_rows > 0) ? r_rows : 1;
         ldl_f = (l_rows > 0) ? l_rows : 1;
         ldx_f = (x_rows > 0) ? x_rows : 1;
         lds_f = (s_t_rows > 0) ? s_t_rows : 1;
         ldt_f = (s_t_rows > 0) ? s_t_rows : 1;
         ldu_f = (u_rows > 0) ? u_rows : 1;

         /* Set pointers */
         a_ptr = a_cm;
         b_ptr = b_cm;
         q_ptr = q_cm;
         r_ptr = (jobb_upper == 'B') ? r_cm : NULL; // Pass NULL if not used
         l_ptr = l_cm;                            // Pass NULL if not used (l_size=0)
         x_ptr = x_cm;
         s_ptr = s_cm;
         t_ptr = t_cm;                            // Pass NULL if not used
         u_ptr = u_cm;

     } else {
         /* Column-major case - use original arrays, potentially copy symmetric inputs */
         lda_f = lda;
         ldb_f = ldb;
         ldq_f = ldq;
         ldr_f = ldr;
         ldl_f = ldl;
         ldx_f = ldx;
         lds_f = lds;
         ldt_f = ldt;
         ldu_f = ldu;
         a_ptr = a;
         b_ptr = b;
         q_ptr = q;
         r_ptr = r;
         l_ptr = l;
         x_ptr = x;
         s_ptr = s;
         t_ptr = t;
         u_ptr = u;

         // Need to copy symmetric inputs G, Q, R to ensure full matrix is passed
         if (jobb_upper == 'G' && b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); slicot_copy_symmetric_part(b, b_cm, n, uplo_upper, ldb, elem_size); b_ptr = b_cm; ldb_f = n;}
         if ((fact_upper == 'N' || fact_upper == 'D') && q_size > 0) { q_cm = (double*)malloc(q_size * elem_size); CHECK_ALLOC(q_cm); slicot_copy_symmetric_part(q, q_cm, n, uplo_upper, ldq, elem_size); q_ptr = q_cm; ldq_f = n;}
         if (jobb_upper == 'B' && (fact_upper == 'N' || fact_upper == 'C') && r_size > 0) { r_cm = (double*)malloc(r_size * elem_size); CHECK_ALLOC(r_cm); slicot_copy_symmetric_part(r, r_cm, m, uplo_upper, ldr, elem_size); r_ptr = r_cm; ldr_f = m;}

         // Adjust pointers for NULL cases
         if (jobb_upper != 'B') r_ptr = NULL;
         if (!(jobb_upper == 'B' && jobl_upper == 'N')) l_ptr = NULL;
         if (dico_upper == 'C' && jobb_upper == 'G') t_ptr = NULL;
     }


     /* --- Workspace Allocation --- */

     // Allocate IWORK and BWORK
     iwork_size = (jobb_upper == 'B') ? MAX(1, MAX(m, 2 * n)) : MAX(1, 2 * n);
     bwork_size = MAX(1, 2 * n);
     iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
     CHECK_ALLOC(iwork);
     bwork = (int*)malloc((size_t)bwork_size * sizeof(int)); // Use int for LOGICAL
     CHECK_ALLOC(bwork);

     // Calculate the required minimum workspace size based on documentation
     if (jobb_upper == 'G') {
         if (dico_upper == 'C') {
             ldwork = MAX(3, 6 * n);
         } else { // DICO = 'D'
             ldwork = MAX(7 * (2 * n + 1) + 16, 16 * n);
         }
     } else { // JOBB = 'B'
         ldwork = MAX(7 * (2 * n + 1) + 16, MAX(16 * n, MAX(2 * n + m, 3 * m)));
     }

     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure


     /* --- Call the computational routine --- */
     F77_FUNC(sb02od, SB02OD)(&dico_upper, &jobb_upper, &fact_upper, &uplo_upper, &jobl_upper, &sort_upper,
                              &n, &m, &p,
                              a_ptr, &lda_f, b_ptr, &ldb_f, q_ptr, &ldq_f,
                              r_ptr, &ldr_f, l_ptr, &ldl_f,
                              rcond, x_ptr, &ldx_f, alfar, alfai, beta,
                              s_ptr, &lds_f, t_ptr, &ldt_f, u_ptr, &ldu_f,
                              &tol, iwork, dwork, &ldwork, bwork, &info,
                              dico_len, jobb_len, fact_len, uplo_len, jobl_len, sort_len);

     /* --- Copy results back to row-major format if needed --- */
     if (row_major && (info == 0 || info == 6)) { // Copy back even if INFO=6
         if (x_cm) slicot_transpose_symmetric_to_c(x_cm, x, n, uplo_upper, elem_size); // X is symmetric
         if (s_cm) slicot_transpose_to_c(s_cm, s, s_t_rows, s_cols, elem_size);
         if (t_cm) slicot_transpose_to_c(t_cm, t, s_t_rows, t_cols, elem_size);
         if (u_cm) slicot_transpose_to_c(u_cm, u, u_rows, u_cols, elem_size);
         // RCOND, ALFAR, ALFAI, BETA modified directly
     } else if (!row_major && (info == 0 || info == 6)) {
         // Output arrays X, S, T, U and scalars modified directly in original arrays.
         // No copy back needed from temps as temps were only for *input* symmetric matrices.
     }

 cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(iwork);
     free(bwork);
     // Input copies
     free(a_cm);
     free(b_cm); // Might be G copy
     free(q_cm); // Might be Q or C copy
     free(r_cm); // Might be R copy
     free(l_cm);
     // Output copies
     free(x_cm);
     free(s_cm);
     free(t_cm);
     free(u_cm);

     return info;
 }