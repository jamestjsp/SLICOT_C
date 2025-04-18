/**
 * @file tb05ad.c
 * @brief C wrapper implementation for SLICOT routine TB05AD
 *
 * This file provides a C wrapper implementation for the SLICOT routine TB05AD,
 * which computes the frequency response matrix G(freq) = C*(freq*I - A)^-1 * B
 * for a given state-space representation (A,B,C) at a specified complex frequency.
 */

 #include <stdlib.h>
 #include <ctype.h>
 #include <string.h> // For memcpy
 #include <stddef.h> // For size_t
 #include <complex.h> // For C99 complex types
 
 // Include the header file for this wrapper
 #include "tb05ad.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, MIN, CHECK_ALLOC, SLICOT_MEMORY_ERROR, slicot_complex_double, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions
 
 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * Note the COMPLEX*16 arguments are mapped via slicot_complex_double.
  * Hidden lengths for CHARACTER arguments are added at the end.
  */
 extern void F77_FUNC(tb05ad, TB05AD)(
     const char* baleig,     // CHARACTER*1 BALEIG
     const char* inita,      // CHARACTER*1 INITA
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     const int* p,           // INTEGER P
     const slicot_complex_double* freq, // COMPLEX*16 FREQ (in)
     double* a,              // DOUBLE PRECISION A(LDA,*) (in/out)
     const int* lda,         // INTEGER LDA
     double* b,              // DOUBLE PRECISION B(LDB,*) (in/out)
     const int* ldb,         // INTEGER LDB
     double* c,              // DOUBLE PRECISION C(LDC,*) (in/out)
     const int* ldc,         // INTEGER LDC
     double* rcond,          // DOUBLE PRECISION RCOND (output)
     slicot_complex_double* g,    // COMPLEX*16 G(LDG,*) (output)
     const int* ldg,         // INTEGER LDG
     double* evre,           // DOUBLE PRECISION EVRE(*) (output)
     double* evim,           // DOUBLE PRECISION EVIM(*) (output)
     slicot_complex_double* hinvb,// COMPLEX*16 HINVB(LDHINV,*) (output)
     const int* ldhinv,      // INTEGER LDHINV
     int* iwork,             // INTEGER IWORK(*)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     slicot_complex_double* zwork,// COMPLEX*16 ZWORK(*)
     const int* lzwork,      // INTEGER LZWORK
     int* info,              // INTEGER INFO (output)
     int baleig_len,         // Hidden length
     int inita_len           // Hidden length
 );
 
 
 /* C wrapper function definition */
 int slicot_tb05ad(char baleig, char inita, int n, int m, int p,
                   slicot_complex_double freq, /* C complex type */
                   double* a, int lda, double* b, int ldb,
                   double* c, int ldc, double* rcond,
                   slicot_complex_double* g, int ldg, /* C complex type */
                   double* evre, double* evim,
                   slicot_complex_double* hinvb, int ldhinv, /* C complex type */
                   int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = -1; /* Use -1 for workspace query */
     int lzwork = -1; /* Use -1 for workspace query */
     double dwork_query;
     slicot_complex_double zwork_query; /* Use standard complex type */
     double* dwork = NULL;
     int* iwork = NULL;
     slicot_complex_double* zwork = NULL; /* Use standard complex type */
     int iwork_size = 0;
 
     const int baleig_len = 1, inita_len = 1;
 
     char baleig_upper = toupper(baleig);
     char inita_upper = toupper(inita);
 
     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL;
     slicot_complex_double *g_cm = NULL, *hinvb_cm = NULL; /* Use standard complex type */
 
     /* --- Input Parameter Validation --- */
     if (n < 0) { info = -3; goto cleanup; }
     if (m < 0) { info = -4; goto cleanup; }
     if (p < 0) { info = -5; goto cleanup; }
     if (baleig_upper != 'N' && baleig_upper != 'C' && baleig_upper != 'B' &&
         baleig_upper != 'E' && baleig_upper != 'A') { info = -1; goto cleanup; }
     if (inita_upper != 'G' && inita_upper != 'H') { info = -2; goto cleanup; }
     // Further consistency checks (e.g., baleig='A' requires inita='G') could be added.
 
     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, n);
     int min_ldb_f = MAX(1, n);
     int min_ldc_f = MAX(1, p);
     int min_ldg_f = MAX(1, p);
     int min_ldhinv_f = MAX(1, n);
 
     if (row_major) {
         // For row-major C, LD is the number of columns
         int min_lda_rm_cols = n;
         int min_ldb_rm_cols = m;
         int min_ldc_rm_cols = n;
         int min_ldg_rm_cols = m;
         int min_ldhinv_rm_cols = m;
         if (lda < min_lda_rm_cols) { info = -8; goto cleanup; }
         if (ldb < min_ldb_rm_cols) { info = -10; goto cleanup; }
         if (ldc < min_ldc_rm_cols) { info = -12; goto cleanup; }
         if (ldg < min_ldg_rm_cols) { info = -15; goto cleanup; }
         if (ldhinv < min_ldhinv_rm_cols) { info = -19; goto cleanup; }
     } else {
         // For column-major C, LD is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -8; goto cleanup; }
         if (ldb < min_ldb_f) { info = -10; goto cleanup; }
         if (ldc < min_ldc_f) { info = -12; goto cleanup; }
         if (ldg < min_ldg_f) { info = -15; goto cleanup; }
         if (ldhinv < min_ldhinv_f) { info = -19; goto cleanup; }
     }
 
     /* --- Workspace Allocation --- */
 
     // Allocate IWORK
     iwork_size = n; // Size is N
     if (iwork_size < 1) iwork_size = 1;
     iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
     CHECK_ALLOC(iwork);
 
     // Perform workspace queries for DWORK and ZWORK
     ldwork = -1;
     lzwork = -1;
     slicot_complex_double freq_f = freq; // Assign C complex to standard complex type
 
     F77_FUNC(tb05ad, TB05AD)(&baleig_upper, &inita_upper, &n, &m, &p, &freq_f,
                              NULL, &lda, NULL, &ldb, NULL, &ldc, rcond,
                              NULL, &ldg, evre, evim, NULL, &ldhinv,
                              iwork, &dwork_query, &ldwork, &zwork_query, &lzwork, &info,
                              baleig_len, inita_len);
 
     if (info < 0 && info != -22 && info != -24) { info = info; goto cleanup; } // Query failed
     info = 0; // Reset info after query
 
     // Get required sizes from query results
     ldwork = (int)dwork_query;
     // Use cabs for complex magnitude if needed, but ZWORK query returns size in zwork_query.real
     // The Fortran interface likely returns the size in the real part of the complex query variable.
     // However, the standard way is to return it in DWORK(1) or similar.
     // Let's rely on the documented minimum sizes if the query seems ambiguous.
     // Fortran documentation for LZWORK query is often unclear in wrappers.
     // Let's calculate minimums based on documentation.
     if (baleig_upper == 'C' || baleig_upper == 'A') {
         lzwork = MAX(1, n*n + 2*n);
     } else {
         lzwork = MAX(1, n*n);
     }
 
     // Check DWORK minimum size based on documentation
     int min_ldwork = 1;
     if (inita_upper == 'G') {
         if (baleig_upper == 'N' || baleig_upper == 'B' || baleig_upper == 'E') {
             min_ldwork = MAX(1, n - 1 + MAX(n, MAX(m, p)));
         } else { // 'C' or 'A'
             min_ldwork = MAX(1, n + MAX(n, MAX(m - 1, p - 1)));
         }
     } else { // 'H'
         if (baleig_upper == 'C' || baleig_upper == 'A') {
             min_ldwork = MAX(1, 2 * n);
         } else { // 'N', 'B', 'E' not applicable or min is 1
             min_ldwork = 1;
         }
     }
     ldwork = MAX(ldwork, min_ldwork);
 
 
     // Allocate DWORK and ZWORK
     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork);
     zwork = (slicot_complex_double*)malloc((size_t)lzwork * sizeof(slicot_complex_double));
     CHECK_ALLOC(zwork);
 
     /* --- Prepare Arrays and Call Fortran Routine --- */
     size_t elem_size_double = sizeof(double);
     size_t elem_size_complex = sizeof(slicot_complex_double);
 
     // Determine sizes for potential copies
     size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
     size_t b_rows = n; size_t b_cols = m; size_t b_size = b_rows * b_cols;
     size_t c_rows = p; size_t c_cols = n; size_t c_size = c_rows * c_cols;
     size_t g_rows = p; size_t g_cols = m; size_t g_size = g_rows * g_cols;
     size_t hinvb_rows = n; size_t hinvb_cols = m; size_t hinvb_size = hinvb_rows * hinvb_cols;
 
     if (row_major) {
         /* --- Row-Major Case --- */
 
         /* Allocate memory for column-major copies */
         if (inita_upper == 'G') { // Only need copies if A, B, C are modified
             if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size_double); CHECK_ALLOC(a_cm); }
             if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size_double); CHECK_ALLOC(b_cm); }
             if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size_double); CHECK_ALLOC(c_cm); }
         }
         // Outputs G and HINVB always need temp storage for copy-back
         if (g_size > 0) { g_cm = (slicot_complex_double*)malloc(g_size * elem_size_complex); CHECK_ALLOC(g_cm); }
         if (hinvb_size > 0) { hinvb_cm = (slicot_complex_double*)malloc(hinvb_size * elem_size_complex); CHECK_ALLOC(hinvb_cm); }
 
         /* Transpose C inputs to Fortran copies (only if INITA='G') */
         double* a_ptr = a;
         double* b_ptr = b;
         double* c_ptr = c;
         if (inita_upper == 'G') {
             if (a_size > 0) { slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size_double); a_ptr = a_cm; }
             if (b_size > 0) { slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, elem_size_double); b_ptr = b_cm; }
             if (c_size > 0) { slicot_transpose_to_fortran(c, c_cm, c_rows, c_cols, elem_size_double); c_ptr = c_cm; }
         }
 
         /* Fortran leading dimensions */
         int lda_f = (a_rows > 0) ? a_rows : 1;
         int ldb_f = (b_rows > 0) ? b_rows : 1;
         int ldc_f = (c_rows > 0) ? c_rows : 1;
         int ldg_f = (g_rows > 0) ? g_rows : 1;
         int ldhinv_f = (hinvb_rows > 0) ? hinvb_rows : 1;
 
         /* Call the Fortran routine */
         F77_FUNC(tb05ad, TB05AD)(&baleig_upper, &inita_upper, &n, &m, &p, &freq_f,
                                  a_ptr, &lda_f, b_ptr, &ldb_f, c_ptr, &ldc_f, rcond,
                                  g_cm, &ldg_f, evre, evim, hinvb_cm, &ldhinv_f,
                                  iwork, dwork, &ldwork, zwork, &lzwork, &info,
                                  baleig_len, inita_len);
 
         /* Copy back results from column-major temps to original row-major arrays */
         if (info == 0 || info == 1 || info == 2) { // Copy back even on warnings/errors
             if (inita_upper == 'G') { // Copy back modified A, B, C
                  if (a_size > 0) slicot_transpose_to_c(a_cm, a, a_rows, a_cols, elem_size_double);
                  if (b_size > 0) slicot_transpose_to_c(b_cm, b, b_rows, b_cols, elem_size_double);
                  if (c_size > 0) slicot_transpose_to_c(c_cm, c, c_rows, c_cols, elem_size_double);
             }
             // Copy back complex outputs G and HINVB
             if (g_size > 0) slicot_transpose_to_c(g_cm, g, g_rows, g_cols, elem_size_complex);
             if (hinvb_size > 0) slicot_transpose_to_c(hinvb_cm, hinvb, hinvb_rows, hinvb_cols, elem_size_complex);
             // RCOND, EVRE, EVIM are modified directly
         }
         /* Temps freed in cleanup */
 
     } else {
         /* --- Column-Major Case --- */
         slicot_complex_double freq_f = freq; // Assign C complex to standard complex type
         /* Call the Fortran routine directly with user-provided arrays */
         F77_FUNC(tb05ad, TB05AD)(&baleig_upper, &inita_upper, &n, &m, &p, &freq_f,
                                  a, &lda, b, &ldb, c, &ldc, rcond,
                                  g, &ldg, /* Pass C complex* directly */
                                  evre, evim,
                                  hinvb, &ldhinv, /* Pass C complex* directly */
                                  iwork, dwork, &ldwork, zwork, &lzwork, &info,
                                  baleig_len, inita_len);
         // A, B, C, RCOND, G, EVRE, EVIM, HINVB modified in place.
     }
 
  cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(iwork);
     free(zwork);
     free(a_cm);
     free(b_cm);
     free(c_cm);
     free(g_cm);
     free(hinvb_cm);
 
     return info;
 }
 