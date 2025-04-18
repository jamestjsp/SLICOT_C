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
 #include "tb04ad.h"
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
     double* b,              // DOUBLE PRECISION B(LDB,*) (in/out)
     const int* ldb,         // INTEGER LDB
     double* c,              // DOUBLE PRECISION C(LDC,*) (in/out)
     const int* ldc,         // INTEGER LDC
     double* d,              // DOUBLE PRECISION D(LDD,*) (in) - Fortran modifies workspace part if ROWCOL='C'
     const int* ldd,         // INTEGER LDD
     int* nr,                // INTEGER NR (output)
     int* index,             // INTEGER INDEX(*) (output)
     double* dcoeff,         // DOUBLE PRECISION DCOEFF(LDDCOE,*) (output)
     const int* lddcoe,      // INTEGER LDDCOE
     double* ucoeff,         // DOUBLE PRECISION UCOEFF(LDUCO1,LDUCO2,*) (output)
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
 
     const int rowcol_len = 1;
 
     char rowcol_upper = toupper(rowcol);
 
     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL;
     // Output arrays dcoeff, ucoeff are filled directly
 
     /* Determine dimensions based on ROWCOL */
     int porm = (rowcol_upper == 'R') ? p : m; // Size for INDEX, DCOEFF dim
     int porp = (rowcol_upper == 'R') ? m : p; // Size for UCOEFF dim
     int maxmp = MAX(m, p);                  // Max dimension for workspace arrays B, C, D
 
     /* --- Input Parameter Validation --- */
     if (n < 0) { info = -2; goto cleanup; }
     if (m < 0) { info = -3; goto cleanup; }
     if (p < 0) { info = -4; goto cleanup; }
     if (rowcol_upper != 'R' && rowcol_upper != 'C') { info = -1; goto cleanup; }
     // TOL1, TOL2 checks done by Fortran
 
     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, n);
     int min_ldb_f = MAX(1, n);
     int min_ldc_f = (rowcol_upper == 'R') ? MAX(1, p) : MAX(1, maxmp); // C needs workspace if 'C' and m > p
     int min_ldd_f = (rowcol_upper == 'R') ? MAX(1, p) : MAX(1, maxmp); // D needs workspace if 'C'
     int min_lddcoe_f = MAX(1, porm);
     int min_lduco1_f = MAX(1, porm);
     int min_lduco2_f = MAX(1, porp);
 
 
     if (row_major) {
         // For row-major C, LD is the number of columns
         int min_lda_rm_cols = n;
         int min_ldb_rm_cols = (rowcol_upper == 'R') ? m : maxmp; // B needs workspace if 'C' and p > m
         int min_ldc_rm_cols = n;
         int min_ldd_rm_cols = (rowcol_upper == 'R') ? m : maxmp; // D needs workspace if 'C'
         int min_lddcoe_rm_rows = porm; // 1st dim is rows
         int min_lduco1_rm_rows = porm; // 1st dim is rows
         int min_lduco2_rm_cols = porp; // 2nd dim is cols
 
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
 
     /* --- Workspace Allocation --- */
 
     // Allocate IWORK
     iwork_size = n + maxmp;
     iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
     CHECK_ALLOC(iwork);
 
     // Allocate DWORK based on query
     ldwork = -1; // Query mode
     F77_FUNC(tb04ad, TB04AD)(&rowcol_upper, &n, &m, &p,
                              NULL, &lda, NULL, &ldb, NULL, &ldc, NULL, &ldd,
                              nr, index, NULL, &lddcoe, NULL, &lduco1, &lduco2,
                              &tol1, &tol2, iwork, &dwork_query, &ldwork, &info,
                              rowcol_len);
 
     if (info < 0 && info != -24) { info = info; goto cleanup; } // Query failed due to invalid argument
     info = 0; // Reset info after query
 
     // Get the required dwork size from query result
     ldwork = (int)dwork_query;
     // Check against minimum documented size
     int min_ldwork_mp = (rowcol_upper == 'R') ? m : p;
     int min_ldwork_pm = (rowcol_upper == 'R') ? p : m;
     int term1 = n * (n + 1);
     int term2 = n * min_ldwork_mp + 2 * n + MAX(n, min_ldwork_mp);
     int term3 = 3 * min_ldwork_mp;
     int term4 = min_ldwork_pm;
     int min_ldwork = MAX(1, MAX(term1 + MAX(term2, term3), term4));
     ldwork = MAX(ldwork, min_ldwork);
 
     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure
 
     /* --- Prepare Arrays and Call Fortran Routine --- */
     size_t elem_size = sizeof(double);
 
     // Determine sizes for potential copies (max dimensions needed for workspace)
     size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
     size_t b_rows = n; size_t b_cols = (rowcol_upper == 'R') ? m : maxmp; size_t b_size = b_rows * b_cols;
     size_t c_rows = (rowcol_upper == 'R') ? p : maxmp; size_t c_cols = n; size_t c_size = c_rows * c_cols;
     size_t d_rows = (rowcol_upper == 'R') ? p : maxmp; size_t d_cols = (rowcol_upper == 'R') ? m : maxmp; size_t d_size = d_rows * d_cols;
 
 
     if (row_major) {
         /* --- Row-Major Case --- */
 
         /* Allocate memory for column-major copies of inputs A, B, C, D */
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
         if (d_size > 0) { d_cm = (double*)malloc(d_size * elem_size); CHECK_ALLOC(d_cm); }
 
         /* Transpose C inputs to Fortran copies */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         // Only transpose the relevant N x M part of B
         if (n > 0 && m > 0) slicot_transpose_to_fortran(b, b_cm, n, m, elem_size);
         // Only transpose the relevant P x N part of C
         if (p > 0 && n > 0) slicot_transpose_to_fortran(c, c_cm, p, n, elem_size);
         // Only transpose the relevant P x M part of D
         if (p > 0 && m > 0) slicot_transpose_to_fortran(d, d_cm, p, m, elem_size);
 
         /* Fortran leading dimensions (use max dimensions for workspace arrays) */
         int lda_f = (a_rows > 0) ? a_rows : 1;
         int ldb_f = (b_rows > 0) ? b_rows : 1;
         int ldc_f = (c_rows > 0) ? c_rows : 1;
         int ldd_f = (d_rows > 0) ? d_rows : 1;
         int lddcoe_f = (porm > 0) ? porm : 1;
         int lduco1_f = (porm > 0) ? porm : 1;
         int lduco2_f = (porp > 0) ? porp : 1;
 
 
         /* Call the Fortran routine */
         F77_FUNC(tb04ad, TB04AD)(&rowcol_upper, &n, &m, &p,
                                  a_cm, &lda_f, b_cm, &ldb_f, c_cm, &ldc_f, d_cm, &ldd_f,
                                  nr, index, dcoeff, &lddcoe_f, ucoeff, &lduco1_f, &lduco2_f,
                                  &tol1, &tol2, iwork, dwork, &ldwork, &info,
                                  rowcol_len);
 
         /* Copy back results from column-major temps to original row-major arrays */
         if (info == 0) { // Only copy back if successful
              int nr_val = *nr; // Get the computed minimal order
              if (nr_val > 0) {
                  if (a_size > 0) slicot_transpose_to_c(a_cm, a, nr_val, nr_val, elem_size);
                  if (b_size > 0 && m > 0) slicot_transpose_to_c(b_cm, b, nr_val, m, elem_size);
                  if (c_size > 0 && p > 0) slicot_transpose_to_c(c_cm, c, p, nr_val, elem_size);
              }
              // D is input only (or workspace), no copy back needed.
              // DCOEFF, UCOEFF are outputs computed directly into user arrays.
              // INDEX, NR, IWORK are outputs computed directly.
         }
         /* Temps freed in cleanup */
 
     } else {
         /* --- Column-Major Case --- */
         /* Call the Fortran routine directly with user-provided arrays */
         F77_FUNC(tb04ad, TB04AD)(&rowcol_upper, &n, &m, &p,
                                  a, &lda, b, &ldb, c, &ldc, d, &ldd,
                                  nr, index, dcoeff, &lddcoe, ucoeff, &lduco1, &lduco2,
                                  &tol1, &tol2, iwork, dwork, &ldwork, &info,
                                  rowcol_len);
         // A, B, C, NR, INDEX, DCOEFF, UCOEFF, IWORK modified in place.
     }
 
  cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(iwork);
     free(a_cm);
     free(b_cm);
     free(c_cm);
     free(d_cm);
 
     return info;
 }
 