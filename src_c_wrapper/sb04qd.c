/**
 * @file sb04qd.c
 * @brief C wrapper implementation for SLICOT routine SB04QD
 *
 * This file provides a C wrapper implementation for the SLICOT routine SB04QD,
 * which solves the discrete-time Sylvester equation X + AXB = C
 * using the Hessenberg-Schur method.
 */

 #include <stdlib.h>
 #include <stddef.h> // For size_t

 // Include the header file for this wrapper
 #include "sb04qd.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * A, B, C are input/output. Z is output.
  */
 extern void F77_FUNC(sb04qd, SB04QD)(
     const int* n,     // INTEGER N
     const int* m,     // INTEGER M
     double* a,        // DOUBLE PRECISION A(LDA,N)
     const int* lda,   // INTEGER LDA
     double* b,        // DOUBLE PRECISION B(LDB,M)
     const int* ldb,   // INTEGER LDB
     double* c,        // DOUBLE PRECISION C(LDC,M) -> X(N,M)
     const int* ldc,   // INTEGER LDC
     double* z,        // DOUBLE PRECISION Z(LDZ,M)
     const int* ldz,   // INTEGER LDZ
     int* iwork,      // INTEGER IWORK(*)
     double* dwork,    // DOUBLE PRECISION DWORK(*)
     const int* ldwork,// INTEGER LDWORK
     int* info         // INTEGER INFO
 );

 /* C wrapper function definition */
 SLICOT_EXPORT
 int slicot_sb04qd(int n, int m,
                   double* a, int lda, double* b, int ldb,
                   double* c, int ldc, double* z, int ldz,
                   int row_major)
 {
     int info = 0;
     double* dwork = NULL;
     int ldwork = 0; // Initialize
     int ldwork_query_flag = -1; // For workspace query

     int* iwork = NULL;
     int iwork_size = 0;

     double* a_cm = NULL;
     double* b_cm = NULL;
     double* c_cm = NULL;
     double* z_cm = NULL;
     double *a_ptr, *b_ptr, *c_ptr, *z_ptr;
     int lda_f, ldb_f, ldc_f, ldz_f;


     /* --- Input Parameter Validation --- */
     if (n < 0) { info = -1; goto cleanup; }
     if (m < 0) { info = -2; goto cleanup; }
     if (a == NULL && n > 0) { info = -3; goto cleanup; }
     // For LDA (4th arg):
     int min_lda_f_val = MAX(1, n);
     if (row_major) { if (n > 0 && lda < n) { info = -4; goto cleanup; } }
     else { if (lda < min_lda_f_val && n > 0) { info = -4; goto cleanup; } }

     if (b == NULL && m > 0) { info = -5; goto cleanup; }
     // For LDB (6th arg):
     int min_ldb_f_val = MAX(1, m);
     if (row_major) { if (m > 0 && ldb < m) { info = -6; goto cleanup; } }
     else { if (ldb < min_ldb_f_val && m > 0) { info = -6; goto cleanup; } }

     if (c == NULL && n > 0 && m > 0) { info = -7; goto cleanup; }
     // For LDC (8th arg):
     int min_ldc_f_val = MAX(1, n);
     if (row_major) { if (m > 0 && ldc < m) { info = -8; goto cleanup; } } // C (NxM) needs M cols in RM
     else { if (ldc < min_ldc_f_val && n > 0 && m > 0) { info = -8; goto cleanup; } }

     if (z == NULL && m > 0) { info = -9; goto cleanup; }
     // For LDZ (10th arg):
     int min_ldz_f_val = MAX(1, m);
     if (row_major) { if (m > 0 && ldz < m) { info = -10; goto cleanup; } } // Z (MxM) needs M cols in RM
     else { if (ldz < min_ldz_f_val && m > 0) { info = -10; goto cleanup; } }

     if (info != 0) { goto cleanup; }

     /* --- Allocate Integer Workspace IWORK --- */
     // From SB04QD.html example, IWORK dimension is 4*N
     if (n > 0) {
         iwork_size = 4 * n;
     } else { // N = 0
         iwork_size = 1; // Minimum size 1, SLICOT often robust to IWORK(1) if N=0
     }
     iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
     CHECK_ALLOC(iwork);

     /* --- Workspace Allocation for DWORK --- */
     double dwork_query_result[1];
     int info_query_ws; // Use a separate info variable for workspace query

     // Use appropriate dummy LDs for the query call, matching Fortran expectations
     int lda_q = MAX(1, n); int ldb_q = MAX(1, m);
     int ldc_q = MAX(1, n); int ldz_q = MAX(1, m);

     F77_FUNC(sb04qd, SB04QD)(&n, &m, NULL, &lda_q, NULL, &ldb_q, NULL, &ldc_q, NULL, &ldz_q,
                              iwork, /* IWORK is needed by SB04QD even for query based on some routines */
                              dwork_query_result, &ldwork_query_flag, &info_query_ws);

     if (info_query_ws == 0) { // Query successful
         ldwork = (int)dwork_query_result[0];
     }
     // Ensure ldwork meets minimum documented formula regardless of query result or if query failed
     // LDWORK = MAX(1, 2*N*N + 9*N, 5*M, N + M)
     int term_2n2_9n = (n > 0) ? (2 * n * n + 9 * n) : 0;
     int term_5m    = (m > 0) ? (5 * m) : 0;
     int term_n_m   = n + m;
     int calc_ldwork = 1; // Start with MAX(1, ...)
     calc_ldwork = MAX(calc_ldwork, term_2n2_9n);
     calc_ldwork = MAX(calc_ldwork, term_5m);
     calc_ldwork = MAX(calc_ldwork, term_n_m);

     if (info_query_ws == 0) { // If query was successful, take the larger of queried and calculated
        ldwork = MAX(ldwork, calc_ldwork);
     } else { // Query failed (e.g. info_query_ws = -13 for LDWORK) or was not supported, use calculated
        ldwork = calc_ldwork;
     }
     ldwork = MAX(1, ldwork); // Final safety ensure at least 1

     if (ldwork > 0) {
         dwork = (double*)malloc((size_t)ldwork * sizeof(double));
         CHECK_ALLOC(dwork);
     } else { // Should not happen if MAX(1,...) is used.
         dwork = NULL; // Or handle as an error.
     }


     /* --- Prepare arrays for column-major format if using row-major --- */
     size_t elem_size = sizeof(double);
     if (row_major) {
         lda_f = MAX(1, n); ldb_f = MAX(1, m);
         ldc_f = MAX(1, n); ldz_f = MAX(1, m);

         size_t a_elems = (size_t)n * n; if (n == 0) a_elems = 0;
         size_t b_elems = (size_t)m * m; if (m == 0) b_elems = 0;
         size_t c_elems = (size_t)n * m; if (n == 0 || m == 0) c_elems = 0;
         size_t z_elems = (size_t)m * m; if (m == 0) z_elems = 0;

         if (a_elems > 0) { a_cm = (double*)malloc(a_elems * elem_size); CHECK_ALLOC(a_cm); slicot_transpose_to_fortran_with_ld(a, a_cm, n, n, lda, lda_f, elem_size); }
         if (b_elems > 0) { b_cm = (double*)malloc(b_elems * elem_size); CHECK_ALLOC(b_cm); slicot_transpose_to_fortran_with_ld(b, b_cm, m, m, ldb, ldb_f, elem_size); }
         if (c_elems > 0) { c_cm = (double*)malloc(c_elems * elem_size); CHECK_ALLOC(c_cm); slicot_transpose_to_fortran_with_ld(c, c_cm, n, m, ldc, ldc_f, elem_size); }
         if (z_elems > 0) { z_cm = (double*)malloc(z_elems * elem_size); CHECK_ALLOC(z_cm); /* Z is output only */ }

         a_ptr = (a_elems > 0) ? a_cm : NULL;
         b_ptr = (b_elems > 0) ? b_cm : NULL;
         c_ptr = (c_elems > 0) ? c_cm : NULL;
         z_ptr = (z_elems > 0) ? z_cm : NULL;
     } else { // Column-major
         a_ptr = a; b_ptr = b; c_ptr = c; z_ptr = z;
         lda_f = lda; ldb_f = ldb; ldc_f = ldc; ldz_f = ldz;

         if (n == 0) a_ptr = NULL;
         if (m == 0) b_ptr = NULL;
         if (n == 0 || m == 0) c_ptr = NULL;
         if (m == 0) z_ptr = NULL;
     }

     // Ensure LDs are at least 1 for dummy/zero-size cases if pointers are not NULL
     if (a_ptr == NULL) lda_f = MAX(1, lda_f); // Or just 1 if ptr is NULL
     if (b_ptr == NULL) ldb_f = MAX(1, ldb_f);
     if (c_ptr == NULL) ldc_f = MAX(1, ldc_f);
     if (z_ptr == NULL) ldz_f = MAX(1, ldz_f);

     lda_f = MAX(1, lda_f); ldb_f = MAX(1, ldb_f);
     ldc_f = MAX(1, ldc_f); ldz_f = MAX(1, ldz_f);


     /* --- Call the computational routine --- */
     F77_FUNC(sb04qd, SB04QD)(&n, &m, a_ptr, &lda_f, b_ptr, &ldb_f, c_ptr, &ldc_f, z_ptr, &ldz_f,
                              iwork, dwork, &ldwork, &info);

     /* --- Copy results back to row-major format if needed --- */
     if (row_major && info == 0) {
         // A and B are input/output per SB04QD.html
         size_t a_elems = (size_t)n*n; if (n==0) a_elems = 0;
         size_t b_elems = (size_t)m*m; if (m==0) b_elems = 0;
         size_t c_elems = (size_t)n * m; if (n == 0 || m == 0) c_elems = 0;
         size_t z_elems = (size_t)m * m; if (m == 0) z_elems = 0;

         if (a_ptr && a_elems > 0) { slicot_transpose_to_c_with_ld(a_ptr, a, n, n, lda_f, lda, elem_size); }
         if (b_ptr && b_elems > 0) { slicot_transpose_to_c_with_ld(b_ptr, b, m, m, ldb_f, ldb, elem_size); }
         if (c_ptr && c_elems > 0) { slicot_transpose_to_c_with_ld(c_ptr, c, n, m, ldc_f, ldc, elem_size); }
         if (z_ptr && z_elems > 0) { slicot_transpose_to_c_with_ld(z_ptr, z, m, m, ldz_f, ldz, elem_size); }
     }

 cleanup:
     free(dwork);
     free(iwork);
     if (row_major) {
        free(a_cm); free(b_cm); free(c_cm); free(z_cm);
     }
     return info;
 }