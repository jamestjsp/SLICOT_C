/**
 * @file tc04ad.c
 * @brief C wrapper implementation for SLICOT routine TC04AD
 *
 * This file provides a C wrapper implementation for the SLICOT routine TC04AD,
 * which finds a state-space representation (A,B,C,D) equivalent to
 * a given left or right polynomial matrix representation.
 */

 #include <stdlib.h>
 #include <ctype.h>
 #include <string.h> // For memcpy
 #include <stddef.h> // For size_t
 #include <stdio.h>  // For diagnostics if needed
 
 // Include the header file for this wrapper
 #include "tc04ad.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, MIN, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions
 
 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * Note the 3D arrays PCOEFF, QCOEFF are passed as flat pointers.
  * Hidden length for CHARACTER argument is added at the end.
  * PCOEFF and QCOEFF are const in C wrapper but non-const in Fortran if LERI='R'
  */
 extern void F77_FUNC(tc04ad, TC04AD)(
     const char* leri,       // CHARACTER*1 LERI
     const int* m,           // INTEGER M
     const int* p,           // INTEGER P
     const int* index,       // INTEGER INDEX(*) (in)
     double* pcoeff,         // DOUBLE PRECISION PCOEFF(LDPCO1,LDPCO2,*) (in) - Fortran modifies if LERI='R'
     const int* ldpco1,      // INTEGER LDPCO1
     const int* ldpco2,      // INTEGER LDPCO2
     double* qcoeff,         // DOUBLE PRECISION QCOEFF(LDQCO1,LDQCO2,*) (in) - Fortran modifies if LERI='R'
     const int* ldqco1,      // INTEGER LDQCO1
     const int* ldqco2,      // INTEGER LDQCO2
     int* n,                 // INTEGER N (output)
     double* rcond,          // DOUBLE PRECISION RCOND (output)
     double* a,              // DOUBLE PRECISION A(LDA,*) (output)
     const int* lda,         // INTEGER LDA
     double* b,              // DOUBLE PRECISION B(LDB,*) (output) - Needs workspace
     const int* ldb,         // INTEGER LDB
     double* c,              // DOUBLE PRECISION C(LDC,*) (output) - Needs workspace
     const int* ldc,         // INTEGER LDC
     double* d,              // DOUBLE PRECISION D(LDD,*) (output) - Needs workspace
     const int* ldd,         // INTEGER LDD
     int* iwork,             // INTEGER IWORK(*)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* info,              // INTEGER INFO (output)
     int leri_len            // Hidden length
 );
 
 
 /* C wrapper function definition */
 int slicot_tc04ad(char leri, int m, int p, const int* index,
                   const double* pcoeff, int ldpco1, int ldpco2,
                   const double* qcoeff, int ldqco1, int ldqco2,
                   int* n, double* rcond,
                   double* a, int lda, double* b, int ldb,
                   double* c, int ldc, double* d, int ldd,
                   int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = -1; /* Use -1 for workspace query */
     double dwork_query;
     double* dwork = NULL;
     int* iwork = NULL;
     int iwork_size = 0;
     int n_val = 0; // To store computed N
 
     const int leri_len = 1;
     char leri_upper = toupper(leri);
 
     /* Pointers for column-major copies if needed */
     double *pcoeff_cm = NULL, *qcoeff_cm = NULL; // For inputs
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL; // For outputs
 
     /* Determine dimensions based on LERI */
     int porm = (leri_upper == 'L') ? p : m; // Dimension for INDEX, PCOEFF
     int porp = (leri_upper == 'L') ? m : p; // Other dimension for QCOEFF
     int maxmp = MAX(m, p);
 
     /* --- Input Parameter Validation --- */
     if (m < 0) { info = -2; goto cleanup; }
     if (p < 0) { info = -3; goto cleanup; }
     if (leri_upper != 'L' && leri_upper != 'R') { info = -1; goto cleanup; }
     if (!index) { info = -4; goto cleanup; } // Check index pointer
 
     // Calculate kpcoef (max degree + 1) and N (sum of degrees)
     int kpcoef = 0;
     n_val = 0;
     if (porm > 0) {
         for (int i = 0; i < porm; ++i) {
             if (index[i] < 0) { info = -4; goto cleanup; } // Degrees must be non-negative
             kpcoef = MAX(kpcoef, index[i]);
             n_val += index[i];
         }
         kpcoef += 1;
     } else {
         kpcoef = 1; // If porm=0, n=0, need kpcoef=1 for array bounds
         n_val = 0;
     }
     if (kpcoef <= 0) kpcoef = 1; // Ensure kpcoef is at least 1
 
     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, n_val);
     int min_ldb_f = MAX(1, n_val);
     int min_ldc_f = MAX(1, maxmp); // C needs workspace
     int min_ldd_f = MAX(1, maxmp); // D needs workspace
     int min_ldpco1_f = MAX(1, porm);
     int min_ldpco2_f = MAX(1, porm);
     int min_ldqco1_f = (leri_upper == 'L') ? MAX(1, p) : MAX(1, maxmp);
     int min_ldqco2_f = (leri_upper == 'L') ? MAX(1, m) : MAX(1, maxmp);
 
 
     if (row_major) {
         // For row-major C, LD is number of columns
         int min_lda_rm_cols = n_val;
         int min_ldb_rm_cols = maxmp; // B needs workspace
         int min_ldc_rm_cols = n_val;
         int min_ldd_rm_cols = maxmp; // D needs workspace
         int min_ldpco1_rm_rows = porm;
         int min_ldpco2_rm_cols = porm;
         int min_ldqco1_rm_rows = (leri_upper == 'L') ? p : maxmp;
         int min_ldqco2_rm_cols = (leri_upper == 'L') ? m : maxmp;
 
         if (lda < min_lda_rm_cols) { info = -14; goto cleanup; }
         if (ldb < min_ldb_rm_cols) { info = -16; goto cleanup; }
         if (ldc < min_ldc_rm_cols) { info = -18; goto cleanup; }
         if (ldd < min_ldd_rm_cols) { info = -20; goto cleanup; }
         if (ldpco1 < min_ldpco1_rm_rows) { info = -6; goto cleanup; }
         if (ldpco2 < min_ldpco2_rm_cols) { info = -7; goto cleanup; }
         if (ldqco1 < min_ldqco1_rm_rows) { info = -9; goto cleanup; }
         if (ldqco2 < min_ldqco2_rm_cols) { info = -10; goto cleanup; }
     } else {
         // For column-major C, LD is number of rows (Fortran style)
         if (lda < min_lda_f) { info = -14; goto cleanup; }
         if (ldb < min_ldb_f) { info = -16; goto cleanup; }
         if (ldc < min_ldc_f) { info = -18; goto cleanup; }
         if (ldd < min_ldd_f) { info = -20; goto cleanup; }
         if (ldpco1 < min_ldpco1_f) { info = -6; goto cleanup; }
         if (ldpco2 < min_ldpco2_f) { info = -7; goto cleanup; }
         if (ldqco1 < min_ldqco1_f) { info = -9; goto cleanup; }
         if (ldqco2 < min_ldqco2_f) { info = -10; goto cleanup; }
     }
 
     /* --- Workspace Allocation --- */
 
     // Allocate IWORK
     iwork_size = 2 * maxmp;
     if (iwork_size < 1) iwork_size = 1;
     iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
     CHECK_ALLOC(iwork);
 
     // Allocate DWORK based on query
     ldwork = -1; // Query mode
     // Need to pass valid pointers even for query, use iwork
     int n_dummy = 0; // Use dummy output vars for query
     double rcond_dummy = 0.0;
     F77_FUNC(tc04ad, TC04AD)(&leri_upper, &m, &p, index,
                              NULL, &ldpco1, &ldpco2, NULL, &ldqco1, &ldqco2,
                              &n_dummy, &rcond_dummy, NULL, &lda, NULL, &ldb, NULL, &ldc, NULL, &ldd,
                              iwork, &dwork_query, &ldwork, &info,
                              leri_len);
 
     if (info < 0 && info != -23) { info = info; goto cleanup; } // Query failed due to invalid argument
     info = 0; // Reset info after query
 
     // Get the required dwork size from query result
     ldwork = (int)dwork_query;
     // Check against minimum documented size
     int min_ldwork = MAX(1, maxmp * (maxmp + 4));
     ldwork = MAX(ldwork, min_ldwork);
 
     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure
 
     /* --- Prepare Arrays and Call Fortran Routine --- */
     size_t elem_size = sizeof(double);
 
     // Calculate total sizes for coefficient arrays
     size_t pcoeff_slice_size = (size_t)porm * porm;
     size_t qcoeff_in_slice_size = (leri_upper == 'L') ? (size_t)p * m : (size_t)p * m; // Q is PxM or PxM? Example implies PxM for LERI='L'. Doc is confusing. Let's assume PxM.
     size_t pcoeff_total_size = pcoeff_slice_size * kpcoef;
     size_t qcoeff_total_size = qcoeff_in_slice_size * kpcoef;
 
     // Sizes for output state-space matrices
     size_t a_rows = n_val; size_t a_cols = n_val; size_t a_size = a_rows * a_cols;
     size_t b_rows = n_val; size_t b_cols = maxmp; size_t b_size = b_rows * b_cols; // Use max for workspace
     size_t c_rows = maxmp; size_t c_cols = n_val; size_t c_size = c_rows * c_cols; // Use max for workspace
     size_t d_rows = maxmp; size_t d_cols = maxmp; size_t d_size = d_rows * d_cols; // Use max for workspace
 
 
     // Need mutable copies of pcoeff, qcoeff for Fortran call, even if input logically
     // Allocate memory for mutable copies
     if (pcoeff_total_size > 0) { pcoeff_cm = (double*)malloc(pcoeff_total_size * elem_size); CHECK_ALLOC(pcoeff_cm); }
     if (qcoeff_total_size > 0) { qcoeff_cm = (double*)malloc(qcoeff_total_size * elem_size); CHECK_ALLOC(qcoeff_cm); }
 
 
     if (row_major) {
         /* --- Row-Major Case --- */
 
         /* Transpose each 2D slice of input coefficients */
         for (int k = 0; k < kpcoef; ++k) {
             size_t offset = k * pcoeff_slice_size;
             if (pcoeff_slice_size > 0) {
                 slicot_transpose_to_fortran(pcoeff + offset, pcoeff_cm + offset, porm, porm, elem_size);
             }
              offset = k * qcoeff_in_slice_size;
              if (qcoeff_in_slice_size > 0) {
                  // Assuming Q is PxM for transposition
                  slicot_transpose_to_fortran(qcoeff + offset, qcoeff_cm + offset, p, m, elem_size);
              }
         }
 
         /* Allocate memory for column-major outputs */
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
         if (d_size > 0) { d_cm = (double*)malloc(d_size * elem_size); CHECK_ALLOC(d_cm); }
 
         /* Fortran leading dimensions */
         int ldpco1_f = (porm > 0) ? porm : 1;
         int ldpco2_f = (porm > 0) ? porm : 1;
         int ldqco1_f = (p > 0) ? p : 1; // Assuming PxM
         int ldqco2_f = (p > 0) ? p : 1; // Assuming PxM? Example uses MAXMP. Let's use P.
         int lda_f = (a_rows > 0) ? a_rows : 1;
         int ldb_f = (b_rows > 0) ? b_rows : 1;
         int ldc_f = (c_rows > 0) ? c_rows : 1;
         int ldd_f = (d_rows > 0) ? d_rows : 1;
 
         /* Call the Fortran routine */
         F77_FUNC(tc04ad, TC04AD)(&leri_upper, &m, &p, index,
                                  pcoeff_cm, &ldpco1_f, &ldpco2_f, qcoeff_cm, &ldqco1_f, &ldqco2_f,
                                  n, rcond, a_cm, &lda_f, b_cm, &ldb_f, c_cm, &ldc_f, d_cm, &ldd_f,
                                  iwork, dwork, &ldwork, &info,
                                  leri_len);
 
         /* Copy back results from column-major temps to original row-major arrays */
         if (info == 0) { // Only copy back if successful
              // N and RCOND are modified directly
              if (a_size > 0) slicot_transpose_to_c(a_cm, a, a_rows, a_cols, elem_size);
              // Copy back only the relevant N x M part of B
              if (b_size > 0 && m > 0) slicot_transpose_to_c(b_cm, b, n_val, m, elem_size);
              // Copy back only the relevant P x N part of C
              if (c_size > 0 && p > 0) slicot_transpose_to_c(c_cm, c, p, n_val, elem_size);
              // Copy back only the relevant P x M part of D
              if (d_size > 0 && p > 0 && m > 0) slicot_transpose_to_c(d_cm, d, p, m, elem_size);
         }
         /* Temps freed in cleanup */
 
     } else {
         /* --- Column-Major Case --- */
         /* Copy const inputs to mutable temps */
         if (pcoeff_total_size > 0) memcpy(pcoeff_cm, pcoeff, pcoeff_total_size * elem_size);
         if (qcoeff_total_size > 0) memcpy(qcoeff_cm, qcoeff, qcoeff_total_size * elem_size);
 
         /* Call the Fortran routine directly with user-provided arrays */
         F77_FUNC(tc04ad, TC04AD)(&leri_upper, &m, &p, index,
                                  pcoeff_cm, &ldpco1, &ldpco2, qcoeff_cm, &ldqco1, &ldqco2,
                                  n, rcond, a, &lda, b, &ldb, c, &ldc, d, &ldd,
                                  iwork, dwork, &ldwork, &info,
                                  leri_len);
         // N, RCOND, A, B, C, D modified in place.
     }
 
  cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(iwork);
     free(pcoeff_cm); // Free mutable copies
     free(qcoeff_cm);
     free(a_cm);      // Free output temps
     free(b_cm);
     free(c_cm);
     free(d_cm);
 
     return info;
 }
 