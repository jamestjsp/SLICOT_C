/**
 * @file ab08nd.c
 * @brief C wrapper implementation for SLICOT routine AB08ND
 *
 * This file provides a C implementation of the SLICOT routine AB08ND wrapper,
 * which constructs a regular pencil for a given system and computes its
 * invariant zeros and Kronecker indices.
 *
 * Note: Explicitly passes Fortran hidden string length argument for 'equil'.
 */

 #include <stdlib.h>
 #include <math.h>
 #include <ctype.h> // For toupper
 
 #include "ab08nd.h"
 #include "slicot_utils.h" // Assumed to contain MAX, MIN, SLICOT_MEMORY_ERROR, CHECK_ALLOC, transpose routines
 #include "slicot_f77.h" /* For F77 call conventions */
 
 
 /*
  * Declare the FORTRAN routine
  * - Uses F77_FUNC for name mangling.
  * - Explicitly includes the hidden 'int' argument at the end
  * for the length of the 'equil' CHARACTER*1 argument.
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
                   double tol, int row_major) {
 
     int info = 0;
     int ldwork = -1; /* For workspace query */
     double dwork_query;
     double *dwork = NULL;
     int *iwork = NULL;
     int iwork_size = 0;
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL;
     double *af_cm = NULL, *bf_cm = NULL;
     int equil_len = 1; // Fortran expects 1-based length for strings
     char equil_upper = toupper(equil);
 
     /* --- Input Parameter Validation --- */
     // Get minimum Fortran leading dimensions from documentation
     int min_lda_f = MAX(1, n);
     int min_ldb_f = MAX(1, n);
     int min_ldc_f = MAX(1, p);
     int min_ldd_f = MAX(1, p);
     int min_ldaf_f = MAX(1, n + m); // From AB08ND docs
     int min_ldbf_f = MAX(1, n + p); // From AB08ND docs
 
     // For row-major C, the user provides LDA as the number of columns.
     int min_lda_rm_cols = n;
     int min_ldb_rm_cols = m;
     int min_ldc_rm_cols = n;
     int min_ldd_rm_cols = m;
 
     if (n < 0) { info = -2; goto cleanup; }
     if (m < 0) { info = -3; goto cleanup; }
     if (p < 0) { info = -4; goto cleanup; }
     if (equil_upper != 'S' && equil_upper != 'N') { info = -1; goto cleanup; }
 
     if (row_major) {
         if (lda < min_lda_rm_cols) { info = -6; goto cleanup; }
         if (ldb < min_ldb_rm_cols) { info = -8; goto cleanup; }
         if (ldc < min_ldc_rm_cols) { info = -10; goto cleanup; }
         if (ldd < min_ldd_rm_cols) { info = -12; goto cleanup; }
         // For row-major, ldaf/ldbf are number of columns in C array.
         // Check is done later after NU is known if they are sufficient.
     } else {
         if (lda < min_lda_f) { info = -6; goto cleanup; }
         if (ldb < min_ldb_f) { info = -8; goto cleanup; }
         if (ldc < min_ldc_f) { info = -10; goto cleanup; }
         if (ldd < min_ldd_f) { info = -12; goto cleanup; }
     }
     // Check LDAF/LDBF against minimum Fortran requirements
     // These are passed directly as Fortran LDs.
     if (ldaf < min_ldaf_f) { info = -20; goto cleanup; }
     if (ldbf < min_ldbf_f) { info = -22; goto cleanup; }
 
 
     /* --- Workspace Query --- */
     // Determine iwork size from AB08ND documentation: MAX(M,P)
     iwork_size = MAX(m, p);
     if (iwork_size < 1) iwork_size = 1; // Ensure at least size 1
     iwork = (int*)malloc(iwork_size * sizeof(int));
     // Use CHECK_ALLOC macro (assumed defined in slicot_utils.h)
     // It sets info = SLICOT_MEMORY_ERROR and jumps to cleanup on failure.
     CHECK_ALLOC(iwork);
 
     // Call Fortran routine for workspace query
     ldwork = -1; // Query mode
     F77_FUNC(ab08nd, AB08ND)(&equil_upper, &n, &m, &p,
                              NULL, &lda, NULL, &ldb, // NULL arrays
                              NULL, &ldc, NULL, &ldd,
                              nu, rank, dinfz, nkror, nkrol, // Pass original pointers
                              NULL, NULL, NULL, // NULL arrays
                              NULL, &ldaf, NULL, &ldbf, // NULL arrays
                              &tol, iwork, &dwork_query, &ldwork, &info, // ldwork = -1
                              equil_len /* Explicit length for equil */);
 
     if (info != 0) {
         // Query failed, likely due to invalid N, M, P, LDA etc. passed to query
         goto cleanup;
     }
 
     /* --- Allocate Workspace --- */
     ldwork = (int)dwork_query;
     if (ldwork < 1) ldwork = 1; // Ensure at least size 1
     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Use CHECK_ALLOC macro
 
     /* --- Prepare Arrays and Call Fortran Routine --- */
     if (row_major) {
         /* Allocate memory for column-major copies */
         size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
         size_t b_rows = n; size_t b_cols = m; size_t b_size = b_rows * b_cols;
         size_t c_rows = p; size_t c_cols = n; size_t c_size = c_rows * c_cols;
         size_t d_rows = p; size_t d_cols = m; size_t d_size = d_rows * d_cols;
         size_t af_fort_rows = ldaf; size_t af_fort_cols = n + MIN(p, m); size_t af_size = af_fort_rows * af_fort_cols;
         size_t bf_fort_rows = ldbf; size_t bf_fort_cols = n + m; size_t bf_size = bf_fort_rows * bf_fort_cols;
 
         // Allocate memory using CHECK_ALLOC
         if (a_size > 0) { a_cm = (double*)malloc(a_size * sizeof(double)); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * sizeof(double)); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * sizeof(double)); CHECK_ALLOC(c_cm); }
         if (d_size > 0) { d_cm = (double*)malloc(d_size * sizeof(double)); CHECK_ALLOC(d_cm); }
         if (af_size > 0) { af_cm = (double*)malloc(af_size * sizeof(double)); CHECK_ALLOC(af_cm); }
         if (bf_size > 0) { bf_cm = (double*)malloc(bf_size * sizeof(double)); CHECK_ALLOC(bf_cm); }
 
         /* Convert row-major inputs to column-major */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, sizeof(double));
         if (b_size > 0) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, sizeof(double));
         if (c_size > 0) slicot_transpose_to_fortran(c, c_cm, c_rows, c_cols, sizeof(double));
         if (d_size > 0) slicot_transpose_to_fortran(d, d_cm, d_rows, d_cols, sizeof(double));
 
         /* Fortran leading dimensions for the call */
         int lda_f = (a_rows > 0) ? a_rows : 1;
         int ldb_f = (b_rows > 0) ? b_rows : 1;
         int ldc_f = (c_rows > 0) ? c_rows : 1;
         int ldd_f = (d_rows > 0) ? d_rows : 1;
         int ldaf_f = ldaf; // Use passed C ldaf as Fortran LDAF
         int ldbf_f = ldbf; // Use passed C ldbf as Fortran LDBF
 
         /* Call FORTRAN routine */
         F77_FUNC(ab08nd, AB08ND)(&equil_upper, &n, &m, &p,
                                  a_cm, &lda_f, b_cm, &ldb_f,
                                  c_cm, &ldc_f, d_cm, &ldd_f,
                                  nu, rank, dinfz, nkror, nkrol,
                                  infz, kronr, kronl,
                                  af_cm, &ldaf_f, bf_cm, &ldbf_f,
                                  &tol, iwork, dwork, &ldwork, &info,
                                  equil_len);
 
         /* Convert column-major results back to row-major if needed */
         if (info == 0) {
              if (*nu > 0) {
                  size_t nu_val = (size_t)(*nu);
                  // Check target C array dimensions (ldaf/ldbf = #cols >= nu)
                  // and source Fortran LDs (ldaf_f/ldbf_f = #rows >= nu)
                  if (ldaf >= *nu && ldaf_f >= *nu) {
                      slicot_transpose_to_c(af_cm, af, nu_val, nu_val, sizeof(double));
                  } else {
                      info = -20; // Map custom error to Fortran-style: Invalid LDAF for computed NU
                      goto cleanup; // Go to cleanup since we cannot proceed
                  }
                  if (ldbf >= *nu && ldbf_f >= *nu) {
                      slicot_transpose_to_c(bf_cm, bf, nu_val, nu_val, sizeof(double));
                  } else {
                      info = -22; // Map custom error: Invalid LDBF for computed NU
                      goto cleanup; // Go to cleanup
                  }
              }
              /* Update original input matrices if changed by the routine (EQUIL='S') */
              if (equil_upper == 'S') { // Check info already done
                  if (a_size > 0) slicot_transpose_to_c(a_cm, a, a_rows, a_cols, sizeof(double));
                  if (b_size > 0) slicot_transpose_to_c(b_cm, b, b_rows, b_cols, sizeof(double));
                  if (c_size > 0) slicot_transpose_to_c(c_cm, c, c_rows, c_cols, sizeof(double));
                  if (d_size > 0) slicot_transpose_to_c(d_cm, d, d_rows, d_cols, sizeof(double));
              }
         }
         /* Temporary arrays a_cm etc. will be freed in cleanup */
 
     } else {
         /* Call FORTRAN routine directly with column-major arrays */
         F77_FUNC(ab08nd, AB08ND)(&equil_upper, &n, &m, &p,
                                   a, &lda, b, &ldb,
                                   c, &ldc, d, &ldd,
                                   nu, rank, dinfz, nkror, nkrol,
                                   infz, kronr, kronl,
                                   af, &ldaf, bf, &ldbf,
                                   &tol, iwork, dwork, &ldwork, &info,
                                   equil_len);
     }
 
 cleanup:
     /* --- Cleanup --- */
     // Free allocated memory (free(NULL) is safe)
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
 