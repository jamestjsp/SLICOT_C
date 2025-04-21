/**
 * @file ab08nd.c
 * @brief C wrapper implementation for SLICOT routine AB08ND
 *
 * This file provides a C implementation of the SLICOT routine AB08ND wrapper,
 * which constructs a regular pencil for a given system and computes its
 * invariant zeros and Kronecker indices.
 * Refactored to align with ab01nd.c structure.
 */

 #include <stdlib.h>
 #include <math.h>   // For isnan, isinf (optional checks)
 #include <ctype.h>  // For toupper
 #include <stddef.h> // For size_t

 // Include the header file for this wrapper
 #include "ab08nd.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, MIN, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

 /*
  * Declare the external FORTRAN routine using the F77_FUNC macro.
  * Explicitly includes the hidden 'int' argument for the length of 'equil'.
  */
 extern void F77_FUNC(ab08nd, AB08ND)(
     const char* equil, const int* n, const int* m, const int* p,
     double* a, const int* lda, double* b, const int* ldb,
     double* c, const int* ldc, double* d, const int* ldd,
     int* nu, int* rank, int* dinfz, int* nkror, int* nkrol,
     int* infz, int* kronr, int* kronl,
     double* af, const int* ldaf, double* bf, const int* ldbf,
     const double* tol, int* iwork, double* dwork, const int* ldwork,
     int* info,
     int equil_len /* Hidden length argument for equil */);

 /* C wrapper for AB08ND */
 SLICOT_C_WRAPPER_API
 int slicot_ab08nd(char equil, int n, int m, int p,
                   double* a, int lda,
                   double* b, int ldb,
                   double* c, int ldc,
                   double* d, int ldd,
                   int* nu, int* rank,
                   int* dinfz, int* nkror, int* nkrol,
                   int* infz, int* kronr, int* kronl,
                   double* af, int ldaf,
                   double* bf, int ldbf,
                   double tol, int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = -1; /* Use -1 for workspace query */
     double dwork_query;
     double *dwork = NULL;
     int *iwork = NULL;
     int iwork_size = 0;
     const int equil_len = 1; // Fortran expects 1-based length for strings

     char equil_upper = toupper(equil);

     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL;
     double *af_cm = NULL, *bf_cm = NULL;

     /* Pointers to pass to Fortran */
     double *a_ptr, *b_ptr, *c_ptr, *d_ptr;
     double *af_ptr, *bf_ptr;
     int lda_f, ldb_f, ldc_f, ldd_f;
     int ldaf_f, ldbf_f;

     /* --- Input Parameter Validation --- */

     if (n < 0) { info = -2; goto cleanup; }
     if (m < 0) { info = -3; goto cleanup; }
     if (p < 0) { info = -4; goto cleanup; }
     if (equil_upper != 'S' && equil_upper != 'N') { info = -1; goto cleanup; }
     // Optional: Check TOL range if necessary
     // if (tol < 0.0) { info = -23; goto cleanup; }

     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, n);
     int min_ldb_f = MAX(1, n);
     int min_ldc_f = MAX(1, p);
     int min_ldd_f = MAX(1, p);
     int min_ldaf_f = MAX(1, n + m); // Fortran routine requires LDAF >= N+M
     int min_ldbf_f = MAX(1, n + p); // Fortran routine requires LDBF >= N+P

     if (row_major) {
         // For row-major C, LDA/LDB/LDC/LDD are the number of columns
         int min_lda_rm_cols = n;
         int min_ldb_rm_cols = m;
         int min_ldc_rm_cols = n;
         int min_ldd_rm_cols = m;
         // LDAF/LDBF (number of columns) check depends on NU, done after Fortran call.
         // However, the *values* passed for ldaf/ldbf must be >= Fortran's row requirements
         // because the wrapper uses them as rows for internal arrays.
         if (n > 0 && lda < min_lda_rm_cols) { info = -6; goto cleanup; }
         if (n > 0 && ldb < min_ldb_rm_cols) { info = -8; goto cleanup; }
         if (p > 0 && ldc < min_ldc_rm_cols) { info = -10; goto cleanup; }
         if (p > 0 && ldd < min_ldd_rm_cols) { info = -12; goto cleanup; }
     } else {
         // For column-major C, LDA/LDB/LDC/LDD are the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -6; goto cleanup; }
         if (ldb < min_ldb_f) { info = -8; goto cleanup; }
         if (ldc < min_ldc_f) { info = -10; goto cleanup; }
         if (ldd < min_ldd_f) { info = -12; goto cleanup; }
         // Check LDAF/LDBF against minimum Fortran requirements (rows)
         if (ldaf < min_ldaf_f) { info = -20; goto cleanup; }
         if (ldbf < min_ldbf_f) { info = -22; goto cleanup; }
     }

     /* --- Workspace Allocation --- */

     // Determine iwork size from AB08ND documentation: MAX(M,P)
     iwork_size = MAX(1, MAX(m, p)); // Ensure minimum size 1
     iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
     CHECK_ALLOC(iwork);

     // Call Fortran routine for DWORK workspace query
     ldwork = -1; // Query mode
     // Use dummy LDs for query if dimensions are 0
     int lda_q = row_major ? MAX(1, n) : lda;
     int ldb_q = row_major ? MAX(1, n) : ldb;
     int ldc_q = row_major ? MAX(1, p) : ldc;
     int ldd_q = row_major ? MAX(1, p) : ldd;
     int ldaf_q = ldaf; // Pass user-provided LDs
     int ldbf_q = ldbf;

     F77_FUNC(ab08nd, AB08ND)(&equil_upper, &n, &m, &p,
                              NULL, &lda_q, NULL, &ldb_q, // NULL arrays
                              NULL, &ldc_q, NULL, &ldd_q,
                              nu, rank, dinfz, nkror, nkrol, // Pass original pointers
                              NULL, NULL, NULL, // NULL integer arrays
                              NULL, &ldaf_q, NULL, &ldbf_q, // NULL output arrays
                              &tol, iwork, &dwork_query, &ldwork, &info, // ldwork = -1
                              equil_len /* Explicit length for equil */);

     if (info != 0) {
         // Query failed, likely due to invalid N, M, P, LDA etc. passed to query
         goto cleanup;
     }

     // Get the required dwork size from query result
     ldwork = (int)dwork_query;
     // Check against minimum documented size: MAX(1, N + MAX(N, M, P))
     int min_ldwork = MAX(1, n + MAX(n, MAX(m, p)));
     ldwork = MAX(ldwork, min_ldwork);

     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Use CHECK_ALLOC macro

     /* --- Prepare Arrays and Call Fortran Routine --- */
     size_t elem_size = sizeof(double);

     if (row_major) {
         /* --- Row-Major Case --- */

         /* Allocate memory for column-major copies */
         size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
         size_t b_rows = n; size_t b_cols = m; size_t b_size = b_rows * b_cols;
         size_t c_rows = p; size_t c_cols = n; size_t c_size = c_rows * c_cols;
         size_t d_rows = p; size_t d_cols = m; size_t d_size = d_rows * d_cols;
         // Fortran AF is LDAF x (N+MIN(P,M)), BF is LDBF x (N+M)
         // Allocate based on Fortran LDs and max possible columns
         size_t af_fort_rows = ldaf; size_t af_fort_cols = n + MIN(p, m); size_t af_size = af_fort_rows * af_fort_cols;
         size_t bf_fort_rows = ldbf; size_t bf_fort_cols = n + m; size_t bf_size = bf_fort_rows * bf_fort_cols;

         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
         if (d_size > 0) { d_cm = (double*)malloc(d_size * elem_size); CHECK_ALLOC(d_cm); }
         if (af_size > 0) { af_cm = (double*)malloc(af_size * elem_size); CHECK_ALLOC(af_cm); }
         if (bf_size > 0) { bf_cm = (double*)malloc(bf_size * elem_size); CHECK_ALLOC(bf_cm); }

         /* Convert row-major inputs to column-major */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         if (b_size > 0) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, elem_size);
         if (c_size > 0) slicot_transpose_to_fortran(c, c_cm, c_rows, c_cols, elem_size);
         if (d_size > 0) slicot_transpose_to_fortran(d, d_cm, d_rows, d_cols, elem_size);

         /* Fortran leading dimensions for the call */
         lda_f = MAX(1, a_rows);
         ldb_f = MAX(1, b_rows);
         ldc_f = MAX(1, c_rows);
         ldd_f = MAX(1, d_rows);
         ldaf_f = ldaf; // Use passed C ldaf as Fortran LDAF
         ldbf_f = ldbf; // Use passed C ldbf as Fortran LDBF

         /* Set pointers for Fortran call */
         a_ptr = a_cm; b_ptr = b_cm; c_ptr = c_cm; d_ptr = d_cm;
         af_ptr = af_cm; bf_ptr = bf_cm;

     } else {
         /* --- Column-Major Case --- */
         lda_f = lda; ldb_f = ldb; ldc_f = ldc; ldd_f = ldd;
         ldaf_f = ldaf; ldbf_f = ldbf;
         a_ptr = a; b_ptr = b; c_ptr = c; d_ptr = d;
         af_ptr = af; bf_ptr = bf;
     }

     /* Call the computational routine */
     F77_FUNC(ab08nd, AB08ND)(&equil_upper, &n, &m, &p,
                              a_ptr, &lda_f, b_ptr, &ldb_f,
                              c_ptr, &ldc_f, d_ptr, &ldd_f,
                              nu, rank, dinfz, nkror, nkrol,
                              infz, kronr, kronl,
                              af_ptr, &ldaf_f, bf_ptr, &ldbf_f,
                              &tol, iwork, dwork, &ldwork, &info,
                              equil_len);

     /* Copy back results if row_major */
     if (row_major && info == 0) {
          int nu_val = *nu;
          // Check if C caller allocated enough columns for the nu x nu result
          if (ldaf < nu_val) { info = -20; goto cleanup; } // Invalid LDAF (cols) for computed NU
          if (ldbf < nu_val) { info = -22; goto cleanup; } // Invalid LDBF (cols) for computed NU

          // Copy back AF and BF
          if (nu_val > 0) {
              // Fortran AF is ldaf_f x nu_val, BF is ldbf_f x nu_val
              // ldaf_f = C ldaf (rows allocated for af_cm)
              // ldbf_f = C ldbf (rows allocated for bf_cm)
              size_t af_src_rows = ldaf_f; size_t af_src_cols = nu_val;
              size_t bf_src_rows = ldbf_f; size_t bf_src_cols = nu_val;
              if (af_src_rows * af_src_cols > 0) {
                  slicot_transpose_to_c(af_cm, af, af_src_rows, af_src_cols, elem_size);
              }
              if (bf_src_rows * bf_src_cols > 0) {
                  slicot_transpose_to_c(bf_cm, bf, bf_src_rows, bf_src_cols, elem_size);
              }
          }
          /* Update original input matrices if changed by the routine (EQUIL='S') */
          if (equil_upper == 'S') {
              size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
              size_t b_rows = n; size_t b_cols = m; size_t b_size = b_rows * b_cols;
              size_t c_rows = p; size_t c_cols = n; size_t c_size = c_rows * c_cols;
              size_t d_rows = p; size_t d_cols = m; size_t d_size = d_rows * d_cols;

              if (a_size > 0) slicot_transpose_to_c(a_cm, a, a_rows, a_cols, elem_size);
              if (b_size > 0) slicot_transpose_to_c(b_cm, b, b_rows, b_cols, elem_size);
              if (c_size > 0) slicot_transpose_to_c(c_cm, c, c_rows, c_cols, elem_size);
              if (d_size > 0) slicot_transpose_to_c(d_cm, d, d_rows, d_cols, elem_size);
          }
     }

 cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(iwork);
     free(a_cm);
     free(b_cm);
     free(c_cm);
     free(d_cm);
     free(af_cm);
     free(bf_cm);

     return info;
 }
