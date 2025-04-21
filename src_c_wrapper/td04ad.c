/**
 * @file td04ad.c
 * @brief C wrapper implementation for SLICOT routine TD04AD
 *
 * This file provides a C wrapper implementation for the SLICOT routine TD04AD,
 * which finds a minimal state-space representation (A,B,C,D) for a
 * proper transfer matrix T(s) given as row or column polynomial
 * vectors over denominator polynomials.
 */

 #include <stdlib.h>
 #include <ctype.h>
 #include <string.h> // For memcpy
 #include <stddef.h> // For size_t

 // Include the header file for this wrapper
 // #include "td04ad.h" // Assuming a header file exists
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, MIN, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * Note the 3D array UCOEFF is passed as a flat pointer.
  * Hidden length for CHARACTER argument is added at the end.
  * UCOEFF is const in C wrapper but non-const in Fortran if ROWCOL='C'.
  */
 extern void F77_FUNC(td04ad, TD04AD)(
     const char* rowcol,     // CHARACTER*1 ROWCOL
     const int* m,           // INTEGER M
     const int* p,           // INTEGER P
     const int* index,       // INTEGER INDEX(*) (in)
     const double* dcoeff,   // DOUBLE PRECISION DCOEFF(LDDCOE,*) (in)
     const int* lddcoe,      // INTEGER LDDCOE
     double* ucoeff,         // DOUBLE PRECISION UCOEFF(LDUCO1,LDUCO2,*) (in) - Fortran modifies if ROWCOL='C'
     const int* lduco1,      // INTEGER LDUCO1
     const int* lduco2,      // INTEGER LDUCO2
     int* nr,                // INTEGER NR (output)
     double* a,              // DOUBLE PRECISION A(LDA,*) (output)
     const int* lda,         // INTEGER LDA
     double* b,              // DOUBLE PRECISION B(LDB,*) (output) - Needs workspace
     const int* ldb,         // INTEGER LDB
     double* c,              // DOUBLE PRECISION C(LDC,*) (output) - Needs workspace
     const int* ldc,         // INTEGER LDC
     double* d,              // DOUBLE PRECISION D(LDD,*) (output) - Needs workspace
     const int* ldd,         // INTEGER LDD
     const double* tol,      // DOUBLE PRECISION TOL
     int* iwork,             // INTEGER IWORK(*) (output - block orders)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* info,              // INTEGER INFO (output)
     int rowcol_len          // Hidden length
 );


 /* C wrapper function definition */
 SLICOT_C_WRAPPER_API
 int slicot_td04ad(char rowcol, int m, int p, const int* index,
                   const double* dcoeff, int lddcoe,
                   const double* ucoeff, int lduco1, int lduco2,
                   int* nr,
                   double* a, int lda, double* b, int ldb,
                   double* c, int ldc, double* d, int ldd,
                   double tol, int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = -1; /* Use -1 for workspace query */
     double dwork_query;
     double* dwork = NULL;
     int* iwork = NULL;
     int iwork_size = 0;
     int n_calc = 0; // Max possible order N = sum(index)
     int kdcoef = 0; // Max degree + 1

     const int rowcol_len = 1;
     char rowcol_upper = toupper(rowcol);

     /* Pointers for column-major copies if needed */
     double *dcoeff_cm = NULL, *ucoeff_cm = NULL; // For inputs (need mutable copy for UCOEFF)
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL; // For outputs
     double *a_ptr, *b_ptr, *c_ptr, *d_ptr;
     const double *dcoeff_ptr; // DCOEFF is const input
     double *ucoeff_ptr;      // UCOEFF needs mutable copy
     int lda_f, ldb_f, ldc_f, ldd_f;
     int lddcoe_f, lduco1_f, lduco2_f;

     /* Determine dimensions based on ROWCOL */
     int porm = (rowcol_upper == 'R') ? p : m; // Dimension for INDEX, DCOEFF
     int porp = (rowcol_upper == 'R') ? m : p; // Other dimension for UCOEFF
     int maxmp = MAX(m, p);

     /* --- Input Parameter Validation --- */
     if (rowcol_upper != 'R' && rowcol_upper != 'C') { info = -1; goto cleanup; }
     if (m < 0) { info = -2; goto cleanup; }
     if (p < 0) { info = -3; goto cleanup; }
     if (!index) { info = -4; goto cleanup; } // Check index pointer
     if (!dcoeff) { info = -5; goto cleanup; } // Check dcoeff pointer
     if (!ucoeff) { info = -7; goto cleanup; } // Check ucoeff pointer
     // Check other output pointers
     if (!nr || !a || !b || !c || !d) { info = -99; goto cleanup; } // Custom code for NULL output pointers
     // TOL check done by Fortran

     // Calculate kdcoef (max degree + 1) and N (sum of degrees)
     if (porm > 0) {
         for (int i = 0; i < porm; ++i) {
             if (index[i] < 0) { info = -4; goto cleanup; } // Degrees must be non-negative
             kdcoef = MAX(kdcoef, index[i]);
             n_calc += index[i];
         }
         kdcoef += 1;
     } else {
         kdcoef = 1; // If porm=0, n=0, need kdcoef=1 for array bounds
         n_calc = 0;
     }
     if (kdcoef <= 0) kdcoef = 1; // Ensure kdcoef is at least 1

     // Check leading dimensions based on storage order (use calculated N)
     // Note: Fortran LDA, LDB are based on N=sum(index), but output A, B are NR x NR, NR x M
     int min_lda_f = MAX(1, n_calc);
     int min_ldb_f = MAX(1, n_calc); // Needs workspace if p > m
     int min_ldc_f = MAX(1, maxmp); // Needs workspace
     int min_ldd_f = (rowcol_upper == 'R') ? MAX(1, p) : MAX(1, maxmp); // Needs workspace if 'C'
     // Input coefficients
     int min_lddcoe_f = MAX(1, porm);
     int min_lduco1_f = (rowcol_upper == 'R') ? MAX(1, p) : MAX(1, maxmp); // Needs workspace
     int min_lduco2_f = (rowcol_upper == 'R') ? MAX(1, m) : MAX(1, maxmp); // Needs workspace

     if (row_major) {
         // For row-major C, LD is number of columns
         int min_lda_rm_cols = n_calc; // Use max N for allocation check
         int min_ldb_rm_cols = maxmp; // Needs workspace
         int min_ldc_rm_cols = n_calc; // Use max N for allocation check
         int min_ldd_rm_cols = (rowcol_upper == 'R') ? m : maxmp; // Needs workspace if 'C'
         // Input coefficients (slices)
         int min_lddcoe_rm_rows = porm;
         int min_lduco1_rm_rows = (rowcol_upper == 'R') ? p : maxmp;
         int min_lduco2_rm_cols = (rowcol_upper == 'R') ? m : maxmp;

         if (lddcoe < min_lddcoe_rm_rows) { info = -6; goto cleanup; }
         if (lduco1 < min_lduco1_rm_rows) { info = -8; goto cleanup; }
         if (lduco2 < min_lduco2_rm_cols) { info = -9; goto cleanup; }
         if (lda < min_lda_rm_cols) { info = -12; goto cleanup; }
         if (ldb < min_ldb_rm_cols) { info = -14; goto cleanup; }
         if (ldc < min_ldc_rm_cols) { info = -16; goto cleanup; }
         if (ldd < min_ldd_rm_cols) { info = -18; goto cleanup; }
     } else {
         // For column-major C, LD is number of rows (Fortran style)
         if (lddcoe < min_lddcoe_f) { info = -6; goto cleanup; }
         if (lduco1 < min_lduco1_f) { info = -8; goto cleanup; }
         if (lduco2 < min_lduco2_f) { info = -9; goto cleanup; }
         if (lda < min_lda_f) { info = -12; goto cleanup; }
         if (ldb < min_ldb_f) { info = -14; goto cleanup; }
         if (ldc < min_ldc_f) { info = -16; goto cleanup; }
         if (ldd < min_ldd_f) { info = -18; goto cleanup; }
     }

     /* --- Prepare arrays for column-major format if using row-major --- */
     size_t elem_size = sizeof(double);
     // Calculate total sizes for coefficient arrays
     size_t dcoeff_total_size = (size_t)porm * kdcoef; // DCOEFF is porm x kdcoef
     size_t ucoeff_slice_size = (rowcol_upper == 'R') ? (size_t)p * m : (size_t)maxmp * maxmp; // Use workspace size for copy
     size_t ucoeff_total_size = ucoeff_slice_size * kdcoef;
     // Sizes for output state-space matrices (use calculated N)
     size_t a_rows_f = n_calc; size_t a_cols_f = n_calc; size_t a_size = a_rows_f * a_cols_f;
     size_t b_rows_f = n_calc; size_t b_cols_f = maxmp; size_t b_size = b_rows_f * b_cols_f; // Use max for workspace
     size_t c_rows_f = maxmp; size_t c_cols_f = n_calc; size_t c_size = c_rows_f * c_cols_f; // Use max for workspace
     size_t d_rows_f = (rowcol_upper == 'R') ? p : maxmp; size_t d_cols_f = (rowcol_upper == 'R') ? m : maxmp; size_t d_size = d_rows_f * d_cols_f; // Use max for workspace if 'C'

     // Allocate memory for copies/transpositions
     if (dcoeff_total_size > 0) { dcoeff_cm = (double*)malloc(dcoeff_total_size * elem_size); CHECK_ALLOC(dcoeff_cm); }
     if (ucoeff_total_size > 0) { ucoeff_cm = (double*)malloc(ucoeff_total_size * elem_size); CHECK_ALLOC(ucoeff_cm); }
     dcoeff_ptr = dcoeff_cm; // Use copy for Fortran call
     ucoeff_ptr = ucoeff_cm;  // Use copy for Fortran call

     if (row_major) {
         /* Transpose 2D input DCOEFF */
         if (dcoeff_total_size > 0 && dcoeff_cm) {
             slicot_transpose_to_fortran(dcoeff, dcoeff_cm, porm, kdcoef, elem_size);
         }

         /* Transpose each 2D slice of input UCOEFF */
         int u_in_rows = p; int u_in_cols = m; // Actual input U size
         size_t u_in_slice_size = (size_t)u_in_rows * u_in_cols;
         for (int k = 0; k < kdcoef; ++k) {
              if (u_in_slice_size > 0 && ucoeff_cm) {
                  // Transpose PxM input into potentially larger CM buffer slice
                  slicot_transpose_to_fortran(ucoeff + k * u_in_slice_size, ucoeff_cm + k * ucoeff_slice_size, u_in_rows, u_in_cols, elem_size);
              }
         }

         /* Allocate memory for column-major outputs */
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
         if (d_size > 0) { d_cm = (double*)malloc(d_size * elem_size); CHECK_ALLOC(d_cm); }
         a_ptr = a_cm; b_ptr = b_cm; c_ptr = c_cm; d_ptr = d_cm;

         /* Fortran leading dimensions (use calculated N and workspace sizes) */
         lddcoe_f = (porm > 0) ? porm : 1;
         lduco1_f = (rowcol_upper == 'R') ? ((p > 0) ? p : 1) : ((maxmp > 0) ? maxmp : 1);
         lduco2_f = (rowcol_upper == 'R') ? ((m > 0) ? m : 1) : ((maxmp > 0) ? maxmp : 1);
         lda_f = (n_calc > 0) ? n_calc : 1;
         ldb_f = (n_calc > 0) ? n_calc : 1;
         ldc_f = (maxmp > 0) ? maxmp : 1;
         ldd_f = (rowcol_upper == 'R') ? ((p > 0) ? p : 1) : ((maxmp > 0) ? maxmp : 1);

     } else {
         /* Column-major case - copy const inputs to mutable temps */
         if (dcoeff_total_size > 0 && dcoeff_cm) memcpy(dcoeff_cm, dcoeff, dcoeff_total_size * elem_size);
         if (ucoeff_total_size > 0 && ucoeff_cm) memcpy(ucoeff_cm, ucoeff, ucoeff_total_size * elem_size);
         // Use original output pointers
         a_ptr = a; b_ptr = b; c_ptr = c; d_ptr = d;
         // Use original LDs
         lddcoe_f = lddcoe; lduco1_f = lduco1; lduco2_f = lduco2;
         lda_f = lda; ldb_f = ldb; ldc_f = ldc; ldd_f = ldd;
     }


     /* --- Workspace Allocation --- */

     // Allocate IWORK
     iwork_size = n_calc + maxmp; // Size is N + MAX(M,P)
     if (iwork_size < 1) iwork_size = 1;
     iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
     CHECK_ALLOC(iwork);

     // Perform workspace query for DWORK
     ldwork = -1; // Query mode
     int nr_dummy; // Dummy output for query
     F77_FUNC(td04ad, TD04AD)(&rowcol_upper, &m, &p, index,
                              dcoeff_ptr, &lddcoe_f, ucoeff_ptr, &lduco1_f, &lduco2_f,
                              &nr_dummy, a_ptr, &lda_f, b_ptr, &ldb_f, c_ptr, &ldc_f, d_ptr, &ldd_f,
                              &tol, iwork, &dwork_query, &ldwork, &info,
                              rowcol_len);

     if (info < 0 && info != -22) { info = info; goto cleanup; } // Query failed due to invalid argument
     info = 0; // Reset info after query

     // Get the required dwork size from query result
     ldwork = (int)dwork_query;
     // Check against minimum documented size
     int min_ldwork = MAX(1, n_calc + MAX(n_calc, MAX(3*m, 3*p)));
     ldwork = MAX(ldwork, min_ldwork);

     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure


     /* --- Call the computational routine --- */
     F77_FUNC(td04ad, TD04AD)(&rowcol_upper, &m, &p, index,
                              dcoeff_ptr, &lddcoe_f, ucoeff_ptr, &lduco1_f, &lduco2_f,
                              nr, a_ptr, &lda_f, b_ptr, &ldb_f, c_ptr, &ldc_f, d_ptr, &ldd_f,
                              &tol, iwork, dwork, &ldwork, &info,
                              rowcol_len);

     /* --- Copy results back to row-major format if needed --- */
     if (row_major && info == 0) {
          int nr_val = *nr; // Get the computed minimal order
          // NR and IWORK are modified directly
          if (nr_val > 0) {
              if (a_size > 0 && a_cm) slicot_transpose_to_c(a_cm, a, nr_val, nr_val, elem_size);
              // Copy back only the relevant NR x M part of B
              if (b_size > 0 && m > 0 && b_cm) slicot_transpose_to_c(b_cm, b, nr_val, m, elem_size);
              // Copy back only the relevant P x NR part of C
              if (c_size > 0 && p > 0 && c_cm) slicot_transpose_to_c(c_cm, c, p, nr_val, elem_size);
          }
          // Copy back only the relevant P x M part of D
          if (d_size > 0 && p > 0 && m > 0 && d_cm) slicot_transpose_to_c(d_cm, d, p, m, elem_size);
     }
     // In column-major case, NR, A, B, C, D, IWORK modified in place.

  cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(iwork);
     free(dcoeff_cm); // Free mutable copies / transposed inputs
     free(ucoeff_cm);
     free(a_cm);      // Free output temps
     free(b_cm);
     free(c_cm);
     free(d_cm);

     return info;
 }