/**
 * @file ab13dd.c
 * @brief C wrapper implementation for SLICOT routine AB13DD
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB13DD,
 * which computes the L-infinity norm of a standard or descriptor system.
 */

 #include <stdlib.h>
 #include <ctype.h>
 #include <stddef.h> // For size_t
 
 // Include the header file for this wrapper
 #include "ab13dd.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, MIN, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines, slicot_complex_double
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions
 
 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * Note use of slicot_complex_double for CWORK.
  */
 extern void F77_FUNC(ab13dd, AB13DD)(
     const char* dico,       // CHARACTER*1 DICO
     const char* jobe,       // CHARACTER*1 JOBE
     const char* equil,      // CHARACTER*1 EQUIL
     const char* jobd,       // CHARACTER*1 JOBD
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     const int* p,           // INTEGER P
     double* fpeak,          // DOUBLE PRECISION FPEAK(2) (in/out)
     const double* a,        // DOUBLE PRECISION A(LDA,*)
     const int* lda,         // INTEGER LDA
     const double* e,        // DOUBLE PRECISION E(LDE,*)
     const int* lde,         // INTEGER LDE
     const double* b,        // DOUBLE PRECISION B(LDB,*)
     const int* ldb,         // INTEGER LDB
     const double* c,        // DOUBLE PRECISION C(LDC,*)
     const int* ldc,         // INTEGER LDC
     const double* d,        // DOUBLE PRECISION D(LDD,*)
     const int* ldd,         // INTEGER LDD
     double* gpeak,          // DOUBLE PRECISION GPEAK(2) (output)
     const double* tol,      // DOUBLE PRECISION TOL
     int* iwork,             // INTEGER IWORK(*)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     slicot_complex_double* cwork, // COMPLEX*16 CWORK(*)
     const int* lcwork,      // INTEGER LCWORK
     int* info,              // INTEGER INFO (output)
     int dico_len,           // Hidden length
     int jobe_len,           // Hidden length
     int equil_len,          // Hidden length
     int jobd_len            // Hidden length
 );
 
 
 /* C wrapper function definition */
 int slicot_ab13dd(char dico, char jobe, char equil, char jobd,
                   int n, int m, int p, double* fpeak,
                   const double* a, int lda, const double* e, int lde,
                   const double* b, int ldb, const double* c, int ldc,
                   const double* d, int ldd, double* gpeak,
                   double tol, int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = -1; /* Use -1 for workspace query */
     int lcwork = -1; /* Use -1 for workspace query */
     double dwork_query;
     slicot_complex_double cwork_query;
     double* dwork = NULL;
     slicot_complex_double* cwork = NULL;
     int* iwork = NULL;
     int iwork_size = 0;
     const int dico_len = 1, jobe_len = 1, equil_len = 1, jobd_len = 1;
 
     char dico_upper = toupper(dico);
     char jobe_upper = toupper(jobe);
     char equil_upper = toupper(equil);
     char jobd_upper = toupper(jobd);
 
     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *e_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL;
 
     /* --- Input Parameter Validation --- */
 
     if (n < 0) { info = -5; goto cleanup; }
     if (m < 0) { info = -6; goto cleanup; }
     if (p < 0) { info = -7; goto cleanup; }
     if (dico_upper != 'C' && dico_upper != 'D') { info = -1; goto cleanup; }
     if (jobe_upper != 'G' && jobe_upper != 'I') { info = -2; goto cleanup; }
     if (equil_upper != 'S' && equil_upper != 'N') { info = -3; goto cleanup; }
     if (jobd_upper != 'D' && jobd_upper != 'Z') { info = -4; goto cleanup; }
     if (fpeak[0] < 0.0 || fpeak[1] < 0.0) { info = -8; goto cleanup; } // FPEAK elements must be non-negative
     if (tol < 0.0 || tol >= 1.0) { info = -19; goto cleanup; }
 
 
     // Check leading dimensions based on storage order and flags
     int min_lda_f = MAX(1, n);
     int min_lde_f = (jobe_upper == 'G') ? MAX(1, n) : 1;
     int min_ldb_f = MAX(1, n);
     int min_ldc_f = MAX(1, p);
     int min_ldd_f = (jobd_upper == 'D') ? MAX(1, p) : 1;
 
     if (row_major) {
         // For row-major C, LDA is the number of columns
         int min_lda_rm_cols = n;
         int min_lde_rm_cols = (jobe_upper == 'G') ? n : 1;
         int min_ldb_rm_cols = m;
         int min_ldc_rm_cols = n;
         int min_ldd_rm_cols = (jobd_upper == 'D') ? m : 1;
         if (lda < min_lda_rm_cols) { info = -10; goto cleanup; }
         if (jobe_upper == 'G' && lde < min_lde_rm_cols) { info = -12; goto cleanup; }
         if (ldb < min_ldb_rm_cols) { info = -14; goto cleanup; }
         if (ldc < min_ldc_rm_cols) { info = -16; goto cleanup; }
         if (jobd_upper == 'D' && ldd < min_ldd_rm_cols) { info = -18; goto cleanup; }
     } else {
         // For column-major C, LDA is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -10; goto cleanup; }
         if (lde < min_lde_f) { info = -12; goto cleanup; } // Check even if JOBE='I'
         if (ldb < min_ldb_f) { info = -14; goto cleanup; }
         if (ldc < min_ldc_f) { info = -16; goto cleanup; }
         if (ldd < min_ldd_f) { info = -18; goto cleanup; } // Check even if JOBD='Z'
     }
 
     /* --- Workspace Allocation --- */
 
     // Allocate IWORK (size N)
     iwork_size = MAX(1, n);
     iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
     CHECK_ALLOC(iwork);
 
     // Query DWORK and CWORK sizes
     ldwork = -1; // Query mode
     lcwork = -1; // Query mode
     F77_FUNC(ab13dd, AB13DD)(&dico_upper, &jobe_upper, &equil_upper, &jobd_upper,
                              &n, &m, &p, fpeak, a, &lda, e, &lde, b, &ldb,
                              c, &ldc, d, &ldd, gpeak, &tol, iwork,
                              &dwork_query, &ldwork, // Pass address for query result
                              &cwork_query, &lcwork, // Pass address for query result
                              &info,
                              dico_len, jobe_len, equil_len, jobd_len);
 
     if (info < 0) { goto cleanup; } // Query failed due to invalid argument
     info = 0; // Reset info after query
 
     // Get the required workspace sizes from query results
     ldwork = (int)dwork_query;
     lcwork = (int)SLICOT_COMPLEX_REAL(cwork_query); // Use macro to get real part
 
     // Check against minimum documented sizes
     int min_ldwork = 1; // Placeholder - formula is very complex
     int min_lcwork = 1;
     if (n > 0 && m > 0 && p > 0) { // Basic check if system non-trivial
         min_lcwork = MAX(1, (n + m) * (n + p) + 2 * MIN(p, m) + MAX(p, m));
     }
     // Add complex formula check for min_ldwork if needed, or rely on query
     ldwork = MAX(ldwork, min_ldwork);
     lcwork = MAX(lcwork, min_lcwork);
 
     // Allocate workspaces
     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork);
     cwork = (slicot_complex_double*)malloc((size_t)lcwork * sizeof(slicot_complex_double));
     CHECK_ALLOC(cwork);
 
     /* --- Prepare Arrays and Call Fortran Routine --- */
     size_t elem_size = sizeof(double);
 
     if (row_major) {
         /* --- Row-Major Case --- */
 
         /* Allocate memory for column-major copies */
         size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
         size_t e_rows = n; size_t e_cols = n; size_t e_size = e_rows * e_cols;
         size_t b_rows = n; size_t b_cols = m; size_t b_size = b_rows * b_cols;
         size_t c_rows = p; size_t c_cols = n; size_t c_size = c_rows * c_cols;
         size_t d_rows = p; size_t d_cols = m; size_t d_size = d_rows * d_cols;
 
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (jobe_upper == 'G' && e_size > 0) { e_cm = (double*)malloc(e_size * elem_size); CHECK_ALLOC(e_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
         if (jobd_upper == 'D' && d_size > 0) { d_cm = (double*)malloc(d_size * elem_size); CHECK_ALLOC(d_cm); }
 
         /* Transpose C (row-major) inputs to Fortran (column-major) copies */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         if (jobe_upper == 'G' && e_size > 0) slicot_transpose_to_fortran(e, e_cm, e_rows, e_cols, elem_size);
         if (b_size > 0) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, elem_size);
         if (c_size > 0) slicot_transpose_to_fortran(c, c_cm, c_rows, c_cols, elem_size);
         if (jobd_upper == 'D' && d_size > 0) slicot_transpose_to_fortran(d, d_cm, d_rows, d_cols, elem_size);
 
         /* Fortran leading dimensions */
         int lda_f = (a_rows > 0) ? a_rows : 1;
         int lde_f = (e_rows > 0) ? e_rows : 1;
         int ldb_f = (b_rows > 0) ? b_rows : 1;
         int ldc_f = (c_rows > 0) ? c_rows : 1;
         int ldd_f = (d_rows > 0) ? d_rows : 1;
 
         /* Call the Fortran routine */
         F77_FUNC(ab13dd, AB13DD)(&dico_upper, &jobe_upper, &equil_upper, &jobd_upper,
                                  &n, &m, &p, fpeak,
                                  a_cm, &lda_f,                 // Pass CM A
                                  (jobe_upper == 'G' ? e_cm : NULL), &lde_f, // Pass CM E or NULL
                                  b_cm, &ldb_f,                 // Pass CM B
                                  c_cm, &ldc_f,                 // Pass CM C
                                  (jobd_upper == 'D' ? d_cm : NULL), &ldd_f, // Pass CM D or NULL
                                  gpeak, &tol, iwork,
                                  dwork, &ldwork, cwork, &lcwork,
                                  &info,
                                  dico_len, jobe_len, equil_len, jobd_len);
 
         /* No copy-back needed for A, E, B, C, D as they are input only */
         /* FPEAK and GPEAK are modified directly */
         /* Temps will be freed in cleanup */
 
     } else {
         /* --- Column-Major Case --- */
 
         /* Call the Fortran routine directly with user-provided arrays */
         F77_FUNC(ab13dd, AB13DD)(&dico_upper, &jobe_upper, &equil_upper, &jobd_upper,
                                  &n, &m, &p, fpeak,
                                  a, &lda,                      // Pass original A
                                  (jobe_upper == 'G' ? e : NULL), &lde, // Pass original E or NULL
                                  b, &ldb,                      // Pass original B
                                  c, &ldc,                      // Pass original C
                                  (jobd_upper == 'D' ? d : NULL), &ldd, // Pass original D or NULL
                                  gpeak, &tol, iwork,
                                  dwork, &ldwork, cwork, &lcwork,
                                  &info,
                                  dico_len, jobe_len, equil_len, jobd_len);
         // FPEAK and GPEAK are modified directly.
     }
 
 cleanup:
     /* --- Cleanup --- */
     free(cwork);
     free(dwork);
     free(iwork);
     free(a_cm);
     free(e_cm);
     free(b_cm);
     free(c_cm);
     free(d_cm);
 
     return info;
 }
 