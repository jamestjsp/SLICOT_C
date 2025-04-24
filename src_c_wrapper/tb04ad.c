/**
 * @file tb04ad.c
 * @brief C wrapper implementation for SLICOT routine TB04AD
 *
 * This file provides a C wrapper implementation for the SLICOT routine TB04AD,
 * which finds the transfer matrix T(s) of a given state-space
 * representation (A,B,C,D), expressed as row or column polynomial
 * vectors over monic least common denominators.
 */

 #include <stdlib.h>
 #include <ctype.h>
 #include <string.h> // For memcpy
 #include <stddef.h> // For size_t

 // Include the header file for this wrapper
 // #include "tb04ad.h" // Assuming a header file exists
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, MIN, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * Note the 3D array UCOEFF is passed as a flat pointer.
  * Hidden length for CHARACTER argument is added at the end.
  */
 extern void F77_FUNC(tb04ad, TB04AD)(
     const char* rowcol,     // CHARACTER*1 ROWCOL
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     const int* p,           // INTEGER P
     double* a,              // DOUBLE PRECISION A(LDA,*) (in/out)
     const int* lda,         // INTEGER LDA
     double* b,              // DOUBLE PRECISION B(LDB,*) (in/out) - Needs workspace
     const int* ldb,         // INTEGER LDB
     double* c,              // DOUBLE PRECISION C(LDC,*) (in/out) - Needs workspace
     const int* ldc,         // INTEGER LDC
     double* d,              // DOUBLE PRECISION D(LDD,*) (in) - Fortran modifies workspace part if ROWCOL='C'
     const int* ldd,         // INTEGER LDD
     int* nr,                // INTEGER NR (output)
     int* index,             // INTEGER INDEX(*) (output)
     double* dcoeff,         // DOUBLE PRECISION DCOEFF(LDDCOE,*) (output)
     const int* lddcoe,      // INTEGER LDDCOE
     double* ucoeff,         // DOUBLE PRECISION UCOEFF(LDUCO1,LDUCO2,*) (output) - Needs workspace
     const int* lduco1,      // INTEGER LDUCO1
     const int* lduco2,      // INTEGER LDUCO2
     const double* tol1,     // DOUBLE PRECISION TOL1
     const double* tol2,     // DOUBLE PRECISION TOL2
     int* iwork,             // INTEGER IWORK(*) (output - block orders)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* info,              // INTEGER INFO (output)
     int rowcol_len          // Hidden length
 );


 /* C wrapper function definition */
 SLICOT_EXPORT
 int slicot_tb04ad(char rowcol, int n, int m, int p,
                   double* a, int lda, double* b, int ldb,
                   double* c, int ldc, double* d, int ldd,
                   int* nr, int* index,
                   double* dcoeff, int lddcoe,
                   double* ucoeff, int lduco1, int lduco2,
                   double tol1, double tol2, int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = -1; /* Use -1 for workspace query */
     double dwork_query;
     double* dwork = NULL;
     int* iwork = NULL;
     int iwork_size = 0;
     int kpcoef = 0; // Max degree + 1 for output arrays

     const int rowcol_len = 1;
     char rowcol_upper = toupper(rowcol);

     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL;
     double *dcoeff_cm = NULL, *ucoeff_cm = NULL;
     double *a_ptr, *b_ptr, *c_ptr, *d_ptr;
     double *dcoeff_ptr, *ucoeff_ptr;
     int lda_f, ldb_f, ldc_f, ldd_f;
     int lddcoe_f, lduco1_f, lduco2_f;

     /* Determine dimensions based on ROWCOL */
     int porm = (rowcol_upper == 'R') ? p : m; // Size for INDEX, DCOEFF dim
     int porp = (rowcol_upper == 'R') ? m : p; // Size for UCOEFF dim
     int maxmp = MAX(m, p);                  // Max dimension for workspace arrays B, C, D

     /* --- Input Parameter Validation --- */
     if (rowcol_upper != 'R' && rowcol_upper != 'C') { info = -1; goto cleanup; }
     if (n < 0) { info = -2; goto cleanup; }
     if (m < 0) { info = -3; goto cleanup; }
     if (p < 0) { info = -4; goto cleanup; }
     // TOL1, TOL2 checks done by Fortran

     // Check leading dimensions based on storage order
     // Note: Fortran workspace requirements may mandate larger LDs than input/output sizes.
     int min_lda_f = MAX(1, n);
     int min_ldb_f = (rowcol_upper == 'R') ? MAX(1, n) : MAX(1, n); // Needs workspace if 'C' and p > m
     int min_ldc_f = (rowcol_upper == 'R') ? MAX(1, p) : MAX(1, maxmp); // Needs workspace if 'C' and m > p
     int min_ldd_f = (rowcol_upper == 'R') ? MAX(1, p) : MAX(1, maxmp); // Needs workspace if 'C'
     // Output arrays
     int min_lddcoe_f = MAX(1, porm);
     int min_lduco1_f = MAX(1, porm);
     int min_lduco2_f = MAX(1, porp);

     if (row_major) {
         // For row-major C, LD is the number of columns
         int min_lda_rm_cols = n;
         int min_ldb_rm_cols = (rowcol_upper == 'R') ? m : maxmp; // Needs workspace
         int min_ldc_rm_cols = n;
         int min_ldd_rm_cols = (rowcol_upper == 'R') ? m : maxmp; // Needs workspace
         // Output arrays (slices)
         int min_lddcoe_rm_rows = porm;
         int min_lduco1_rm_rows = porm;
         int min_lduco2_rm_cols = porp;

         if (lda < min_lda_rm_cols) { info = -6; goto cleanup; }
         if (ldb < min_ldb_rm_cols) { info = -8; goto cleanup; }
         if (ldc < min_ldc_rm_cols) { info = -10; goto cleanup; }
         if (ldd < min_ldd_rm_cols) { info = -12; goto cleanup; }
         if (lddcoe < min_lddcoe_rm_rows) { info = -16; goto cleanup; }
         if (lduco1 < min_lduco1_rm_rows) { info = -18; goto cleanup; }
         if (lduco2 < min_lduco2_rm_cols) { info = -19; goto cleanup; }

     } else {
         // For column-major C, LD is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -6; goto cleanup; }
         if (ldb < min_ldb_f) { info = -8; goto cleanup; }
         if (ldc < min_ldc_f) { info = -10; goto cleanup; }
         if (ldd < min_ldd_f) { info = -12; goto cleanup; }
         if (lddcoe < min_lddcoe_f) { info = -16; goto cleanup; }
         if (lduco1 < min_lduco1_f) { info = -18; goto cleanup; }
         if (lduco2 < min_lduco2_f) { info = -19; goto cleanup; }
     }

      /* --- Prepare arrays for column-major format if using row-major --- */
     size_t elem_size = sizeof(double);
     if (row_major) {
         // Allocate CM copies for inputs A, B, C, D (considering workspace needs)
         size_t a_rows_f = n; size_t a_cols_f = n; size_t a_size = a_rows_f * a_cols_f;
         size_t b_rows_f = n; size_t b_cols_f = (rowcol_upper == 'R') ? m : maxmp; size_t b_size = b_rows_f * b_cols_f;
         size_t c_rows_f = (rowcol_upper == 'R') ? p : maxmp; size_t c_cols_f = n; size_t c_size = c_rows_f * c_cols_f;
         size_t d_rows_f = (rowcol_upper == 'R') ? p : maxmp; size_t d_cols_f = (rowcol_upper == 'R') ? m : maxmp; size_t d_size = d_rows_f * d_cols_f;

         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
         if (d_size > 0) { d_cm = (double*)malloc(d_size * elem_size); CHECK_ALLOC(d_cm); }

         // Transpose relevant input parts to CM copies
         if (n > 0 && n > 0 && a_cm) slicot_transpose_to_fortran(a, a_cm, n, n, elem_size);
         if (n > 0 && m > 0 && b_cm) slicot_transpose_to_fortran(b, b_cm, n, m, elem_size);
         if (p > 0 && n > 0 && c_cm) slicot_transpose_to_fortran(c, c_cm, p, n, elem_size);
         if (p > 0 && m > 0 && d_cm) slicot_transpose_to_fortran(d, d_cm, p, m, elem_size);

         // Allocate CM copies for outputs DCOEFF, UCOEFF
         // Sizes depend on kpcoef which is unknown until after call. Cannot pre-allocate perfectly.
         // Allocate based on input LDs * estimated max degree (N). This is an upper bound.
         kpcoef = n + 1; // Upper bound for max degree + 1
         size_t dcoeff_sz_est = (size_t)lddcoe * kpcoef;
         size_t ucoeff_sz_est = (size_t)lduco1 * lduco2 * kpcoef;

         if (lddcoe > 0 && kpcoef > 0) { dcoeff_cm = (double*)malloc((size_t)lddcoe * kpcoef * elem_size); CHECK_ALLOC(dcoeff_cm); }
         if (ucoeff_sz_est > 0) { ucoeff_cm = (double*)malloc(ucoeff_sz_est * elem_size); CHECK_ALLOC(ucoeff_cm); }

         // Set Fortran leading dimensions
         lda_f = (n > 0) ? n : 1;
         ldb_f = (n > 0) ? n : 1;
         ldc_f = (rowcol_upper == 'R') ? ((p > 0) ? p : 1) : ((maxmp > 0) ? maxmp : 1);
         ldd_f = (rowcol_upper == 'R') ? ((p > 0) ? p : 1) : ((maxmp > 0) ? maxmp : 1);
         lddcoe_f = (porm > 0) ? porm : 1;
         lduco1_f = (porm > 0) ? porm : 1;
         lduco2_f = (porp > 0) ? porp : 1;

         // Set pointers
         a_ptr = a_cm; b_ptr = b_cm; c_ptr = c_cm; d_ptr = d_cm;
         dcoeff_ptr = dcoeff_cm; ucoeff_ptr = ucoeff_cm;

     } else {
         /* Column-major case - use original arrays */
         lda_f = lda; ldb_f = ldb; ldc_f = ldc; ldd_f = ldd;
         lddcoe_f = lddcoe; lduco1_f = lduco1; lduco2_f = lduco2;
         a_ptr = a; b_ptr = b; c_ptr = c; d_ptr = d;
         dcoeff_ptr = dcoeff; ucoeff_ptr = ucoeff;
     }

     /* --- Workspace allocation --- */

     // Allocate IWORK
     iwork_size = n + maxmp;
     if (iwork_size < 1) iwork_size = 1;
     iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
     CHECK_ALLOC(iwork);

     // Calculate DWORK size directly using the formula
     int mp = (rowcol_upper == 'R') ? m : p;
     int pm = (rowcol_upper == 'R') ? p : m;
     
     int inner_term1 = n * mp + 2 * n + MAX(n, mp);
     int inner_term2 = 3 * mp;
     int inner_term3 = pm;
     
     // LDWORK >= MAX(1, N*(N + 1) + MAX(N*MP + 2*N + MAX(N,MP), 3*MP, PM))
     ldwork = MAX(1, n * (n + 1) + MAX(MAX(inner_term1, inner_term2), inner_term3));

     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure


     /* --- Call the computational routine --- */
     F77_FUNC(tb04ad, TB04AD)(&rowcol_upper, &n, &m, &p,
                              a_ptr, &lda_f, b_ptr, &ldb_f, c_ptr, &ldc_f, d_ptr, &ldd_f,
                              nr, index, dcoeff_ptr, &lddcoe_f, ucoeff_ptr, &lduco1_f, &lduco2_f,
                              &tol1, &tol2, iwork, dwork, &ldwork, &info,
                              rowcol_len);

     /* --- Copy results back to row-major format if needed --- */
     if (row_major && info == 0) {
          int nr_val = *nr; // Get the computed minimal order
          kpcoef = 0;       // Recalculate max degree for output copy
          if (porm > 0 && nr_val > 0) {
               for (int i = 0; i < porm; ++i) {
                   kpcoef = MAX(kpcoef, index[i]); // index is output
               }
               kpcoef += 1;
          } else {
               kpcoef = 1;
          }
          if (kpcoef <= 0) kpcoef = 1;

          // Copy back A, B, C (potentially modified)
          if (nr_val > 0) {
              size_t a_rows_f = n; size_t a_cols_f = n; size_t a_size = a_rows_f * a_cols_f;
              size_t b_rows_f = n; size_t b_cols_f = (rowcol_upper == 'R') ? m : maxmp; size_t b_size = b_rows_f * b_cols_f;
              size_t c_rows_f = (rowcol_upper == 'R') ? p : maxmp; size_t c_cols_f = n; size_t c_size = c_rows_f * c_cols_f;

              if (a_size > 0 && a_cm) slicot_transpose_to_c_with_ld(a_cm, a, nr_val, nr_val, n, lda, elem_size);
              if (b_size > 0 && m > 0 && b_cm) slicot_transpose_to_c_with_ld(b_cm, b, nr_val, m, n, ldb, elem_size);
              if (c_size > 0 && p > 0 && c_cm) slicot_transpose_to_c_with_ld(c_cm, c, p, nr_val, p, ldc, elem_size);
          }

          // Copy back 2D output array DCOEFF
          if (porm > 0 && kpcoef > 0 && dcoeff_cm) {
              slicot_transpose_to_c_with_ld(dcoeff_cm, dcoeff, porm, kpcoef, porm, lddcoe, elem_size);
          }

          // Copy back 3D output array UCOEFF slice by slice
          size_t ucoeff_slice_size = (size_t)porm * porp; // Size of one PORM x PORP slice
          for (int k = 0; k < kpcoef; ++k) {
              if (ucoeff_slice_size > 0 && ucoeff_cm) {
                  slicot_transpose_to_c_with_ld(ucoeff_cm + k * ucoeff_slice_size, 
                                               ucoeff + k * ucoeff_slice_size, 
                                               porm, porp, porm, lduco1, elem_size);
              }
          }
     }
     // In column-major case, A, B, C, NR, INDEX, DCOEFF, UCOEFF, IWORK modified in place.

  cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(iwork);
     free(a_cm); free(b_cm); free(c_cm); free(d_cm);
     free(dcoeff_cm); free(ucoeff_cm);

     return info;
 }