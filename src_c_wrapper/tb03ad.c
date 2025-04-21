/**
 * @file tb03ad.c
 * @brief C wrapper implementation for SLICOT routine TB03AD
 *
 * This file provides a C wrapper implementation for the SLICOT routine TB03AD,
 * which finds a relatively prime left or right polynomial matrix
 * representation equivalent to a given state-space representation.
 */

 #include <stdlib.h>
 #include <ctype.h>
 #include <string.h> // For memcpy
 #include <stddef.h> // For size_t

 // Include the header file for this wrapper
 // #include "tb03ad.h" // Assuming a header file exists
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, MIN, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * Note the 3D arrays PCOEFF, QCOEFF, VCOEFF are passed as flat pointers.
  * Hidden lengths for CHARACTER arguments are added at the end.
  */
 extern void F77_FUNC(tb03ad, TB03AD)(
     const char* leri,       // CHARACTER*1 LERI
     const char* equil,      // CHARACTER*1 EQUIL
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     const int* p,           // INTEGER P
     double* a,              // DOUBLE PRECISION A(LDA,*) (in/out)
     const int* lda,         // INTEGER LDA
     double* b,              // DOUBLE PRECISION B(LDB,*) (in/out) - Needs workspace
     const int* ldb,         // INTEGER LDB
     double* c,              // DOUBLE PRECISION C(LDC,*) (in/out) - Needs workspace
     const int* ldc,         // INTEGER LDC
     double* d,              // DOUBLE PRECISION D(LDD,*) (in) - Fortran modifies workspace part
     const int* ldd,         // INTEGER LDD
     int* nr,                // INTEGER NR (output)
     int* index,             // INTEGER INDEX(*) (output)
     double* pcoeff,         // DOUBLE PRECISION PCOEFF(LDPCO1,LDPCO2,*) (output)
     const int* ldpco1,      // INTEGER LDPCO1
     const int* ldpco2,      // INTEGER LDPCO2
     double* qcoeff,         // DOUBLE PRECISION QCOEFF(LDQCO1,LDQCO2,*) (output) - Needs workspace
     const int* ldqco1,      // INTEGER LDQCO1
     const int* ldqco2,      // INTEGER LDQCO2
     double* vcoeff,         // DOUBLE PRECISION VCOEFF(LDVCO1,LDVCO2,*) (output)
     const int* ldvco1,      // INTEGER LDVCO1
     const int* ldvco2,      // INTEGER LDVCO2
     const double* tol,      // DOUBLE PRECISION TOL
     int* iwork,             // INTEGER IWORK(*) (output - block orders)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* info,              // INTEGER INFO (output)
     int leri_len,           // Hidden length
     int equil_len           // Hidden length
 );


 /* C wrapper function definition */
 SLICOT_C_WRAPPER_API
 int slicot_tb03ad(char leri, char equil, int n, int m, int p,
                   double* a, int lda, double* b, int ldb,
                   double* c, int ldc, double* d, int ldd,
                   int* nr, int* index,
                   double* pcoeff, int ldpco1, int ldpco2,
                   double* qcoeff, int ldqco1, int ldqco2,
                   double* vcoeff, int ldvco1, int ldvco2,
                   double tol, int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = -1; /* Use -1 for workspace query */
     double dwork_query;
     double* dwork = NULL;
     int* iwork = NULL;
     int iwork_size = 0;
     int kpcoef = 0; // Max degree + 1 for output arrays

     const int leri_len = 1, equil_len = 1;
     char leri_upper = toupper(leri);
     char equil_upper = toupper(equil);

     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL;
     double *pcoeff_cm = NULL, *qcoeff_cm = NULL, *vcoeff_cm = NULL;
     double *a_ptr, *b_ptr, *c_ptr, *d_ptr;
     double *pcoeff_ptr, *qcoeff_ptr, *vcoeff_ptr;
     int lda_f, ldb_f, ldc_f, ldd_f;
     int ldpco1_f, ldpco2_f, ldqco1_f, ldqco2_f, ldvco1_f, ldvco2_f;


     /* Determine dimensions based on LERI */
     int porm = (leri_upper == 'L') ? p : m; // Size for INDEX, PCOEFF dims
     int porp = (leri_upper == 'L') ? m : p; // Size for QCOEFF dims
     int maxmp = MAX(m, p);                  // Max dimension for workspace arrays B, C, D

     /* --- Input Parameter Validation --- */
     if (leri_upper != 'L' && leri_upper != 'R') { info = -1; goto cleanup; }
     if (equil_upper != 'S' && equil_upper != 'N') { info = -2; goto cleanup; }
     if (n < 0) { info = -3; goto cleanup; }
     if (m < 0) { info = -4; goto cleanup; }
     if (p < 0) { info = -5; goto cleanup; }
     // TOL check done by Fortran

     // Check leading dimensions based on storage order
     // Note: Fortran workspace requirements may mandate larger LDs than input/output sizes.
     int min_lda_f = MAX(1, n);
     int min_ldb_f = MAX(1, n); // Needs workspace if p > m
     int min_ldc_f = MAX(1, maxmp); // Needs workspace if m > p
     int min_ldd_f = MAX(1, maxmp); // Needs workspace

     // Output array dimensions validation (minimal checks)
     int min_ldpco1_f = MAX(1, porm);
     int min_ldpco2_f = MAX(1, porm);
     int min_ldqco1_f = (leri_upper == 'L') ? MAX(1, p) : MAX(1, maxmp); // Needs workspace
     int min_ldqco2_f = (leri_upper == 'L') ? MAX(1, m) : MAX(1, maxmp); // Needs workspace
     int min_ldvco1_f = MAX(1, porm);
     int min_ldvco2_f = MAX(1, n); // Max dimension needed is N

     if (row_major) {
         // For row-major C, LD is the number of columns
         int min_lda_rm_cols = n;
         int min_ldb_rm_cols = maxmp; // Needs workspace
         int min_ldc_rm_cols = n;
         int min_ldd_rm_cols = maxmp; // Needs workspace

         // Check 3D array dimensions (interpreted as rows/cols of each slice)
         int min_ldpco1_rm_rows = porm;
         int min_ldpco2_rm_cols = porm;
         int min_ldqco1_rm_rows = (leri_upper == 'L') ? p : maxmp;
         int min_ldqco2_rm_cols = (leri_upper == 'L') ? m : maxmp;
         int min_ldvco1_rm_rows = porm;
         int min_ldvco2_rm_cols = n; // Max dimension needed is N

         if (lda < min_lda_rm_cols) { info = -7; goto cleanup; }
         if (ldb < min_ldb_rm_cols) { info = -9; goto cleanup; }
         if (ldc < min_ldc_rm_cols) { info = -11; goto cleanup; }
         if (ldd < min_ldd_rm_cols) { info = -13; goto cleanup; }
         if (ldpco1 < min_ldpco1_rm_rows) { info = -17; goto cleanup; }
         if (ldpco2 < min_ldpco2_rm_cols) { info = -18; goto cleanup; }
         if (ldqco1 < min_ldqco1_rm_rows) { info = -20; goto cleanup; }
         if (ldqco2 < min_ldqco2_rm_cols) { info = -21; goto cleanup; }
         if (ldvco1 < min_ldvco1_rm_rows) { info = -23; goto cleanup; }
         if (ldvco2 < min_ldvco2_rm_cols) { info = -24; goto cleanup; }

     } else {
         // For column-major C, LD is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -7; goto cleanup; }
         if (ldb < min_ldb_f) { info = -9; goto cleanup; }
         if (ldc < min_ldc_f) { info = -11; goto cleanup; }
         if (ldd < min_ldd_f) { info = -13; goto cleanup; }
         if (ldpco1 < min_ldpco1_f) { info = -17; goto cleanup; }
         if (ldpco2 < min_ldpco2_f) { info = -18; goto cleanup; }
         if (ldqco1 < min_ldqco1_f) { info = -20; goto cleanup; }
         if (ldqco2 < min_ldqco2_f) { info = -21; goto cleanup; }
         if (ldvco1 < min_ldvco1_f) { info = -23; goto cleanup; }
         if (ldvco2 < min_ldvco2_f) { info = -24; goto cleanup; }
     }

     /* --- Prepare arrays for column-major format if using row-major --- */
     size_t elem_size = sizeof(double);
     if (row_major) {
         // Allocate CM copies for inputs A, B, C, D (considering workspace needs)
         size_t a_rows_f = n; size_t a_cols_f = n; size_t a_size = a_rows_f * a_cols_f;
         size_t b_rows_f = n; size_t b_cols_f = maxmp; size_t b_size = b_rows_f * b_cols_f;
         size_t c_rows_f = maxmp; size_t c_cols_f = n; size_t c_size = c_rows_f * c_cols_f;
         size_t d_rows_f = maxmp; size_t d_cols_f = maxmp; size_t d_size = d_rows_f * d_cols_f;

         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
         if (d_size > 0) { d_cm = (double*)malloc(d_size * elem_size); CHECK_ALLOC(d_cm); }

         // Transpose relevant input parts to CM copies
         if (n > 0 && n > 0 && a_cm) slicot_transpose_to_fortran(a, a_cm, n, n, elem_size);
         if (n > 0 && m > 0 && b_cm) slicot_transpose_to_fortran(b, b_cm, n, m, elem_size);
         if (p > 0 && n > 0 && c_cm) slicot_transpose_to_fortran(c, c_cm, p, n, elem_size);
         if (p > 0 && m > 0 && d_cm) slicot_transpose_to_fortran(d, d_cm, p, m, elem_size);

         // Allocate CM copies for outputs PCOEFF, QCOEFF, VCOEFF
         // Sizes depend on kpcoef which is unknown until after call. Cannot pre-allocate perfectly.
         // Allocate based on input LDs * estimated max degree (N). This is an upper bound.
         kpcoef = n + 1; // Upper bound for max degree + 1
         size_t pcoeff_sz_est = (size_t)ldpco1 * ldpco2 * kpcoef;
         size_t qcoeff_sz_est = (size_t)ldqco1 * ldqco2 * kpcoef;
         size_t vcoeff_sz_est = (size_t)ldvco1 * ldvco2 * kpcoef;

         if (pcoeff_sz_est > 0) { pcoeff_cm = (double*)malloc(pcoeff_sz_est * elem_size); CHECK_ALLOC(pcoeff_cm); }
         if (qcoeff_sz_est > 0) { qcoeff_cm = (double*)malloc(qcoeff_sz_est * elem_size); CHECK_ALLOC(qcoeff_cm); }
         if (vcoeff_sz_est > 0) { vcoeff_cm = (double*)malloc(vcoeff_sz_est * elem_size); CHECK_ALLOC(vcoeff_cm); }

         // Set Fortran leading dimensions
         lda_f = (n > 0) ? n : 1;
         ldb_f = (n > 0) ? n : 1;
         ldc_f = (maxmp > 0) ? maxmp : 1;
         ldd_f = (maxmp > 0) ? maxmp : 1;
         ldpco1_f = (porm > 0) ? porm : 1;
         ldpco2_f = (porm > 0) ? porm : 1;
         ldqco1_f = (leri_upper == 'L') ? ((p > 0) ? p : 1) : ((maxmp > 0) ? maxmp : 1);
         ldqco2_f = (leri_upper == 'L') ? ((m > 0) ? m : 1) : ((maxmp > 0) ? maxmp : 1);
         ldvco1_f = (porm > 0) ? porm : 1;
         ldvco2_f = (n > 0) ? n : 1;

         // Set pointers
         a_ptr = a_cm; b_ptr = b_cm; c_ptr = c_cm; d_ptr = d_cm;
         pcoeff_ptr = pcoeff_cm; qcoeff_ptr = qcoeff_cm; vcoeff_ptr = vcoeff_cm;

     } else {
         /* Column-major case - use original arrays */
         lda_f = lda; ldb_f = ldb; ldc_f = ldc; ldd_f = ldd;
         ldpco1_f = ldpco1; ldpco2_f = ldpco2;
         ldqco1_f = ldqco1; ldqco2_f = ldqco2;
         ldvco1_f = ldvco1; ldvco2_f = ldvco2;
         a_ptr = a; b_ptr = b; c_ptr = c; d_ptr = d;
         pcoeff_ptr = pcoeff; qcoeff_ptr = qcoeff; vcoeff_ptr = vcoeff;
     }


     /* --- Workspace allocation --- */

     // Allocate IWORK
     iwork_size = n + maxmp;
     if (iwork_size < 1) iwork_size = 1;
     iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
     CHECK_ALLOC(iwork);

     // Perform workspace query for DWORK
     ldwork = -1; // Query mode
     int nr_dummy; // Dummy output for query
     F77_FUNC(tb03ad, TB03AD)(&leri_upper, &equil_upper, &n, &m, &p,
                              a_ptr, &lda_f, b_ptr, &ldb_f, c_ptr, &ldc_f, d_ptr, &ldd_f,
                              &nr_dummy, index, pcoeff_ptr, &ldpco1_f, &ldpco2_f,
                              qcoeff_ptr, &ldqco1_f, &ldqco2_f, vcoeff_ptr, &ldvco1_f, &ldvco2_f,
                              &tol, iwork, &dwork_query, &ldwork, &info,
                              leri_len, equil_len);

     if (info < 0 && info != -28) { info = info; goto cleanup; } // Query failed due to invalid argument
     info = 0; // Reset info after query

     // Get the required dwork size from query result
     ldwork = (int)dwork_query;
     // Check against minimum documented size
     int min_ldwork_pm = (leri_upper == 'L') ? p : m;
     int min_ldwork = MAX(1, MAX(n + MAX(n, MAX(3*m, 3*p)), min_ldwork_pm * (min_ldwork_pm + 2)));
     ldwork = MAX(ldwork, min_ldwork);

     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure

     /* --- Call the computational routine --- */
     F77_FUNC(tb03ad, TB03AD)(&leri_upper, &equil_upper, &n, &m, &p,
                              a_ptr, &lda_f, b_ptr, &ldb_f, c_ptr, &ldc_f, d_ptr, &ldd_f,
                              nr, index, pcoeff_ptr, &ldpco1_f, &ldpco2_f,
                              qcoeff_ptr, &ldqco1_f, &ldqco2_f, vcoeff_ptr, &ldvco1_f, &ldvco2_f,
                              &tol, iwork, dwork, &ldwork, &info,
                              leri_len, equil_len);

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
              size_t b_rows_f = n; size_t b_cols_f = maxmp; size_t b_size = b_rows_f * b_cols_f;
              size_t c_rows_f = maxmp; size_t c_cols_f = n; size_t c_size = c_rows_f * c_cols_f;
              if (a_size > 0 && a_cm) slicot_transpose_to_c(a_cm, a, nr_val, nr_val, elem_size);
              if (b_size > 0 && m > 0 && b_cm) slicot_transpose_to_c(b_cm, b, nr_val, m, elem_size);
              if (c_size > 0 && p > 0 && c_cm) slicot_transpose_to_c(c_cm, c, p, nr_val, elem_size);
          }

          // Copy back 3D output arrays slice by slice
          size_t pcoeff_slice_size = (size_t)porm * porm;
          size_t qcoeff_slice_size = (leri_upper == 'L') ? (size_t)p * m : (size_t)maxmp * maxmp; // Use output Q dims (PxM) or workspace size? Assume PxM.
          size_t vcoeff_slice_size = (size_t)porm * n; // Use max N
          size_t q_rows_out = p; size_t q_cols_out = m; size_t q_out_slice = q_rows_out * q_cols_out;
          size_t v_rows_out = porm; size_t v_cols_out = nr_val; size_t v_out_slice = v_rows_out * v_cols_out;

          for (int k = 0; k < kpcoef; ++k) {
              if (pcoeff_slice_size > 0 && pcoeff_cm) {
                  slicot_transpose_to_c(pcoeff_cm + k * pcoeff_slice_size, pcoeff + k * pcoeff_slice_size, porm, porm, elem_size);
              }
              if (q_out_slice > 0 && qcoeff_cm) {
                  slicot_transpose_to_c(qcoeff_cm + k * qcoeff_slice_size, qcoeff + k * q_out_slice, q_rows_out, q_cols_out, elem_size);
              }
              if (v_out_slice > 0 && vcoeff_cm) {
                  slicot_transpose_to_c(vcoeff_cm + k * vcoeff_slice_size, vcoeff + k * v_out_slice, v_rows_out, v_cols_out, elem_size);
              }
          }
     }
     // In column-major case, A, B, C, NR, INDEX, PCOEFF, QCOEFF, VCOEFF, IWORK modified in place.

  cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(iwork);
     free(a_cm); free(b_cm); free(c_cm); free(d_cm);
     free(pcoeff_cm); free(qcoeff_cm); free(vcoeff_cm);

     return info;
 }