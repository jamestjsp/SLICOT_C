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
 #include "tb03ad.h"
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
     double* b,              // DOUBLE PRECISION B(LDB,*) (in/out)
     const int* ldb,         // INTEGER LDB
     double* c,              // DOUBLE PRECISION C(LDC,*) (in/out)
     const int* ldc,         // INTEGER LDC
     double* d,              // DOUBLE PRECISION D(LDD,*) (in) - Fortran modifies workspace part
     const int* ldd,         // INTEGER LDD
     int* nr,                // INTEGER NR (output)
     int* index,             // INTEGER INDEX(*) (output)
     double* pcoeff,         // DOUBLE PRECISION PCOEFF(LDPCO1,LDPCO2,*) (output)
     const int* ldpco1,      // INTEGER LDPCO1
     const int* ldpco2,      // INTEGER LDPCO2
     double* qcoeff,         // DOUBLE PRECISION QCOEFF(LDQCO1,LDQCO2,*) (output)
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
 
     const int leri_len = 1, equil_len = 1;
 
     char leri_upper = toupper(leri);
     char equil_upper = toupper(equil);
 
     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL;
     double *pcoeff_cm = NULL, *qcoeff_cm = NULL, *vcoeff_cm = NULL;
 
     /* Determine dimensions based on LERI */
     int porm = (leri_upper == 'L') ? p : m; // Size for INDEX, PCOEFF dims
     int porp = (leri_upper == 'L') ? m : p; // Size for QCOEFF dims
     int maxmp = MAX(m, p);                  // Max dimension for workspace arrays B, C, D
 
     /* --- Input Parameter Validation --- */
     if (n < 0) { info = -3; goto cleanup; }
     if (m < 0) { info = -4; goto cleanup; }
     if (p < 0) { info = -5; goto cleanup; }
     if (leri_upper != 'L' && leri_upper != 'R') { info = -1; goto cleanup; }
     if (equil_upper != 'S' && equil_upper != 'N') { info = -2; goto cleanup; }
     // TOL check done by Fortran
 
     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, n);
     int min_ldb_f = MAX(1, n);
     int min_ldc_f = MAX(1, maxmp); // C needs workspace if m > p
     int min_ldd_f = MAX(1, maxmp); // D needs workspace
 
     int min_ldpco1_f = MAX(1, porm);
     int min_ldpco2_f = MAX(1, porm);
     int min_ldqco1_f = (leri_upper == 'L') ? MAX(1, p) : MAX(1, maxmp); // Q needs workspace if m > p and LERI='R'
     int min_ldqco2_f = (leri_upper == 'L') ? MAX(1, m) : MAX(1, maxmp); // Q needs workspace if p > m and LERI='R'
     int min_ldvco1_f = MAX(1, porm);
     int min_ldvco2_f = MAX(1, n); // Fortran doc says N, C wrapper uses NR later, but max N needed
 
     if (row_major) {
         // For row-major C, LD is the number of columns
         int min_lda_rm_cols = n;
         int min_ldb_rm_cols = maxmp; // B needs workspace if p > m
         int min_ldc_rm_cols = n;
         int min_ldd_rm_cols = maxmp; // D needs workspace
 
         int min_ldpco1_rm_rows = porm; // 1st dim is rows in C for 3D
         int min_ldpco2_rm_cols = porm; // 2nd dim is cols in C for 3D
         int min_ldqco1_rm_rows = (leri_upper == 'L') ? p : maxmp;
         int min_ldqco2_rm_cols = (leri_upper == 'L') ? m : maxmp;
         int min_ldvco1_rm_rows = porm;
         int min_ldvco2_rm_cols = n; // Max N needed
 
         if (lda < min_lda_rm_cols) { info = -7; goto cleanup; }
         if (ldb < min_ldb_rm_cols) { info = -9; goto cleanup; }
         if (ldc < min_ldc_rm_cols) { info = -11; goto cleanup; }
         if (ldd < min_ldd_rm_cols) { info = -13; goto cleanup; }
         // Check 3D array dimensions (interpreted as 2D slices)
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
 
     /* --- Workspace Allocation --- */
 
     // Allocate IWORK
     iwork_size = n + maxmp;
     iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
     CHECK_ALLOC(iwork);
 
     // Allocate DWORK based on query
     ldwork = -1; // Query mode
     F77_FUNC(tb03ad, TB03AD)(&leri_upper, &equil_upper, &n, &m, &p,
                              NULL, &lda, NULL, &ldb, NULL, &ldc, NULL, &ldd,
                              nr, index, NULL, &ldpco1, &ldpco2,
                              NULL, &ldqco1, &ldqco2, NULL, &ldvco1, &ldvco2,
                              &tol, iwork, &dwork_query, &ldwork, &info,
                              leri_len, equil_len);
 
     if (info < 0 && info != -28) { info = info; goto cleanup; } // Query failed due to invalid argument
     info = 0; // Reset info after query
 
     // Get the required dwork size from query result
     ldwork = (int)dwork_query;
     // Check against minimum documented size
     int min_ldwork_pm = (leri_upper == 'L') ? p : m;
     // FIX: Properly nest MAX calls for three arguments
     int min_ldwork = MAX(1, MAX(n + MAX(n, MAX(3*m, 3*p)), min_ldwork_pm * (min_ldwork_pm + 2)));
     ldwork = MAX(ldwork, min_ldwork);
 
     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure
 
     /* --- Prepare Arrays and Call Fortran Routine --- */
     size_t elem_size = sizeof(double);
 
     // Determine sizes for potential copies (max dimensions needed for workspace)
     size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
     size_t b_rows = n; size_t b_cols = maxmp; size_t b_size = b_rows * b_cols;
     size_t c_rows = maxmp; size_t c_cols = n; size_t c_size = c_rows * c_cols;
     size_t d_rows = maxmp; size_t d_cols = maxmp; size_t d_size = d_rows * d_cols;
 
     // Sizes for 3D arrays (number of elements per slice)
     // Note: kpcoef (max degree + 1) is unknown until after the call.
     // We cannot pre-allocate/transpose pcoeff, qcoeff, vcoeff easily.
     // The Fortran routine computes these. We only need to allocate space
     // for the row-major case if we were to copy *back* results, but the
     // C function passes the user's pointers directly for output.
     // We *do* need to transpose A, B, C, D inputs for row-major.
 
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
         int ldpco1_f = (porm > 0) ? porm : 1;
         int ldpco2_f = (porm > 0) ? porm : 1;
         int ldqco1_f = (leri_upper == 'L') ? ((p > 0) ? p : 1) : ((maxmp > 0) ? maxmp : 1);
         int ldqco2_f = (leri_upper == 'L') ? ((m > 0) ? m : 1) : ((maxmp > 0) ? maxmp : 1);
         int ldvco1_f = (porm > 0) ? porm : 1;
         int ldvco2_f = (n > 0) ? n : 1; // Use N for input LD, actual output size NR unknown
 
 
         /* Call the Fortran routine */
         F77_FUNC(tb03ad, TB03AD)(&leri_upper, &equil_upper, &n, &m, &p,
                                  a_cm, &lda_f, b_cm, &ldb_f, c_cm, &ldc_f, d_cm, &ldd_f,
                                  nr, index, pcoeff, &ldpco1_f, &ldpco2_f,
                                  qcoeff, &ldqco1_f, &ldqco2_f, vcoeff, &ldvco1_f, &ldvco2_f,
                                  &tol, iwork, dwork, &ldwork, &info,
                                  leri_len, equil_len);
 
         /* Copy back results from column-major temps to original row-major arrays */
         if (info == 0) { // Only copy back if successful
              int nr_val = *nr; // Get the computed minimal order
              if (nr_val > 0) {
                  if (a_size > 0) slicot_transpose_to_c(a_cm, a, nr_val, nr_val, elem_size);
                  if (b_size > 0 && m > 0) slicot_transpose_to_c(b_cm, b, nr_val, m, elem_size);
                  if (c_size > 0 && p > 0) slicot_transpose_to_c(c_cm, c, p, nr_val, elem_size);
              }
              // D is input only, no copy back needed.
              // PCOEFF, QCOEFF, VCOEFF are outputs computed directly into user arrays.
              // INDEX, NR, IWORK are outputs computed directly.
         }
         /* Temps freed in cleanup */
 
     } else {
         /* --- Column-Major Case --- */
         /* Call the Fortran routine directly with user-provided arrays */
         F77_FUNC(tb03ad, TB03AD)(&leri_upper, &equil_upper, &n, &m, &p,
                                  a, &lda, b, &ldb, c, &ldc, d, &ldd,
                                  nr, index, pcoeff, &ldpco1, &ldpco2,
                                  qcoeff, &ldqco1, &ldqco2, vcoeff, &ldvco1, &ldvco2,
                                  &tol, iwork, dwork, &ldwork, &info,
                                  leri_len, equil_len);
         // A, B, C, NR, INDEX, PCOEFF, QCOEFF, VCOEFF, IWORK modified in place.
     }
 
  cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(iwork);
     free(a_cm);
     free(b_cm);
     free(c_cm);
     free(d_cm);
     // No need to free pcoeff_cm etc. as they weren't allocated for copy-back
 
     return info;
 }
 