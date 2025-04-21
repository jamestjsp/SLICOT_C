/**
 * @file ab05nd.c
 * @brief C wrapper implementation for SLICOT routine AB05ND
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB05ND,
 * which computes the state-space model (A,B,C,D) for the feedback
 * inter-connection of two systems.
 * This wrapper assumes OVER = 'N' (no overlapping arrays).
 * Refactored to align with ab01nd.c structure.
 */

 #include <stdlib.h>
 #include <ctype.h>  // For toupper
 #include <stddef.h> // For size_t

 // Include the header file for this wrapper
 #include "ab05nd.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * Use const for input arrays as this wrapper assumes OVER = 'N'.
  * Hidden string length arguments are explicitly included.
  */
 extern void F77_FUNC(ab05nd, AB05ND)(
     const char* over,       // CHARACTER*1 OVER
     const int* n1,          // INTEGER N1
     const int* m1,          // INTEGER M1
     const int* p1,          // INTEGER P1
     const int* n2,          // INTEGER N2
     const double* alpha,    // DOUBLE PRECISION ALPHA
     const double* a1,       // DOUBLE PRECISION A1(LDA1,*)
     const int* lda1,        // INTEGER LDA1
     const double* b1,       // DOUBLE PRECISION B1(LDB1,*)
     const int* ldb1,        // INTEGER LDB1
     const double* c1,       // DOUBLE PRECISION C1(LDC1,*)
     const int* ldc1,        // INTEGER LDC1
     const double* d1,       // DOUBLE PRECISION D1(LDD1,*)
     const int* ldd1,        // INTEGER LDD1
     const double* a2,       // DOUBLE PRECISION A2(LDA2,*)
     const int* lda2,        // INTEGER LDA2
     const double* b2,       // DOUBLE PRECISION B2(LDB2,*)
     const int* ldb2,        // INTEGER LDB2
     const double* c2,       // DOUBLE PRECISION C2(LDC2,*)
     const int* ldc2,        // INTEGER LDC2
     const double* d2,       // DOUBLE PRECISION D2(LDD2,*)
     const int* ldd2,        // INTEGER LDD2
     int* n,                 // INTEGER N (output)
     double* a,              // DOUBLE PRECISION A(LDA,*) (output)
     const int* lda,         // INTEGER LDA
     double* b,              // DOUBLE PRECISION B(LDB,*) (output)
     const int* ldb,         // INTEGER LDB
     double* c,              // DOUBLE PRECISION C(LDC,*) (output)
     const int* ldc,         // INTEGER LDC
     double* d,              // DOUBLE PRECISION D(LDD,*) (output)
     const int* ldd,         // INTEGER LDD
     int* iwork,             // INTEGER IWORK(*)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* info,              // INTEGER INFO (output)
     int over_len            // Hidden length for over
 );


 /* C wrapper function definition */
 SLICOT_C_WRAPPER_API
 int slicot_ab05nd(char over,
                   int n1, int m1, int p1, int n2, double alpha,
                   const double* a1, int lda1, const double* b1, int ldb1,
                   const double* c1, int ldc1, const double* d1, int ldd1,
                   const double* a2, int lda2, const double* b2, int ldb2,
                   const double* c2, int ldc2, const double* d2, int ldd2,
                   int* n_out, // Renamed to avoid conflict with local n
                   double* a, int lda, double* b, int ldb,
                   double* c, int ldc, double* d, int ldd,
                   int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = 1;
     double* dwork = NULL;
     int* iwork = NULL;
     int iwork_size = 0;
     const int over_len = 1;
     int n = 0; // Calculated total state dimension (local variable)

     char over_upper = toupper(over);

     /* Pointers for column-major copies if needed */
     double *a1_cm = NULL, *b1_cm = NULL, *c1_cm = NULL, *d1_cm = NULL;
     double *a2_cm = NULL, *b2_cm = NULL, *c2_cm = NULL, *d2_cm = NULL;
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL;

     /* Pointers to pass to Fortran */
     const double *a1_ptr, *b1_ptr, *c1_ptr, *d1_ptr;
     const double *a2_ptr, *b2_ptr, *c2_ptr, *d2_ptr;
     double *a_ptr, *b_ptr, *c_ptr, *d_ptr;
     int lda1_f, ldb1_f, ldc1_f, ldd1_f;
     int lda2_f, ldb2_f, ldc2_f, ldd2_f;
     int lda_f, ldb_f, ldc_f, ldd_f;

     /* --- Input Parameter Validation --- */

     if (n1 < 0) { info = -2; goto cleanup; }
     if (m1 < 0) { info = -3; goto cleanup; }
     if (p1 < 0) { info = -4; goto cleanup; }
     if (n2 < 0) { info = -5; goto cleanup; }
     // No check for ALPHA needed based on docs

     // Force OVER = 'N' as 'O' is not supported by this wrapper interface
     if (over_upper == 'O') {
         // Optionally print a warning or return an error
         // For now, proceed assuming 'N' was intended for the wrapper call
         over_upper = 'N';
     }
     if (over_upper != 'N') { info = -1; goto cleanup; }

     // Calculate resulting dimension N
     n = n1 + n2;
     if (n_out) { *n_out = n; } // Set output N if pointer is provided

     // Check leading dimensions for inputs
     int min_lda1_f = MAX(1, n1); int min_ldb1_f = MAX(1, n1);
     int min_ldc1_f = MAX(1, p1); int min_ldd1_f = MAX(1, p1); // Rows=P1
     int min_lda2_f = MAX(1, n2); int min_ldb2_f = MAX(1, n2);
     int min_ldc2_f = MAX(1, m1); int min_ldd2_f = MAX(1, m1); // Rows=M1
     // Check leading dimensions for outputs
     int min_lda_f = MAX(1, n); int min_ldb_f = MAX(1, n);
     int min_ldc_f = MAX(1, p1); int min_ldd_f = MAX(1, p1); // Rows=P1

     if (row_major) {
         // Check C dimensions (number of columns)
         if (n1 > 0 && lda1 < n1) { info = -8; goto cleanup; }
         if (n1 > 0 && ldb1 < m1) { info = -10; goto cleanup; }
         if (p1 > 0 && ldc1 < n1) { info = -12; goto cleanup; } // C1 is P1-by-N1
         if (p1 > 0 && ldd1 < m1) { info = -14; goto cleanup; } // D1 is P1-by-M1
         if (n2 > 0 && lda2 < n2) { info = -16; goto cleanup; }
         if (n2 > 0 && ldb2 < p1) { info = -18; goto cleanup; } // B2 is N2-by-P1
         if (m1 > 0 && ldc2 < n2) { info = -20; goto cleanup; } // C2 is M1-by-N2
         if (m1 > 0 && ldd2 < p1) { info = -22; goto cleanup; } // D2 is M1-by-P1
         // Check output C dimensions
         if (n > 0 && lda < n) { info = -25; goto cleanup; }
         if (n > 0 && ldb < m1) { info = -27; goto cleanup; } // B is N-by-M1
         if (p1 > 0 && ldc < n) { info = -29; goto cleanup; } // C is P1-by-N
         if (p1 > 0 && ldd < m1) { info = -31; goto cleanup; } // D is P1-by-M1
     } else {
         // Check Fortran dimensions (number of rows)
         if (lda1 < min_lda1_f) { info = -8; goto cleanup; }
         if (ldb1 < min_ldb1_f) { info = -10; goto cleanup; }
         if (ldc1 < min_ldc1_f) { info = -12; goto cleanup; }
         if (ldd1 < min_ldd1_f) { info = -14; goto cleanup; }
         if (lda2 < min_lda2_f) { info = -16; goto cleanup; }
         if (ldb2 < min_ldb2_f) { info = -18; goto cleanup; }
         if (ldc2 < min_ldc2_f) { info = -20; goto cleanup; }
         if (ldd2 < min_ldd2_f) { info = -22; goto cleanup; }
         // Check output Fortran dimensions
         if (lda < min_lda_f) { info = -25; goto cleanup; }
         if (ldb < min_ldb_f) { info = -27; goto cleanup; }
         if (ldc < min_ldc_f) { info = -29; goto cleanup; }
         if (ldd < min_ldd_f) { info = -31; goto cleanup; }
     }

     /* --- Workspace Allocation --- */

     // Allocate IWORK (size P1)
     iwork_size = MAX(1, p1); // Ensure minimum size 1
     iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
     CHECK_ALLOC(iwork);

     // Allocate DWORK (calculate size for OVER='N')
     // LDWORK >= MAX(1, P1*P1, M1*M1, N1*P1)
     ldwork = 1; // Start with minimum
     if (p1 > 0) ldwork = MAX(ldwork, p1 * p1);
     if (m1 > 0) ldwork = MAX(ldwork, m1 * m1);
     if (n1 > 0 && p1 > 0) ldwork = MAX(ldwork, n1 * p1);
     // Note: If N1=0, P1=0, M1=0, ldwork remains 1.

     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure

     /* --- Prepare Arrays and Call Fortran Routine --- */

     if (row_major) {
         /* --- Row-Major Case --- */

         /* Allocate memory for column-major copies of inputs */
         size_t a1_size = (size_t)n1 * n1; if (n1 == 0) a1_size = 0;
         size_t b1_size = (size_t)n1 * m1; if (n1 == 0 || m1 == 0) b1_size = 0;
         size_t c1_size = (size_t)p1 * n1; if (p1 == 0 || n1 == 0) c1_size = 0; // P1 rows, N1 cols
         size_t d1_size = (size_t)p1 * m1; if (p1 == 0 || m1 == 0) d1_size = 0; // P1 rows, M1 cols
         size_t a2_size = (size_t)n2 * n2; if (n2 == 0) a2_size = 0;
         size_t b2_size = (size_t)n2 * p1; if (n2 == 0 || p1 == 0) b2_size = 0; // N2 rows, P1 cols
         size_t c2_size = (size_t)m1 * n2; if (m1 == 0 || n2 == 0) c2_size = 0; // M1 rows, N2 cols
         size_t d2_size = (size_t)m1 * p1; if (m1 == 0 || p1 == 0) d2_size = 0; // M1 rows, P1 cols

         if (a1_size > 0) { a1_cm = (double*)malloc(a1_size * sizeof(double)); CHECK_ALLOC(a1_cm); }
         if (b1_size > 0) { b1_cm = (double*)malloc(b1_size * sizeof(double)); CHECK_ALLOC(b1_cm); }
         if (c1_size > 0) { c1_cm = (double*)malloc(c1_size * sizeof(double)); CHECK_ALLOC(c1_cm); }
         if (d1_size > 0) { d1_cm = (double*)malloc(d1_size * sizeof(double)); CHECK_ALLOC(d1_cm); }
         if (a2_size > 0) { a2_cm = (double*)malloc(a2_size * sizeof(double)); CHECK_ALLOC(a2_cm); }
         if (b2_size > 0) { b2_cm = (double*)malloc(b2_size * sizeof(double)); CHECK_ALLOC(b2_cm); }
         if (c2_size > 0) { c2_cm = (double*)malloc(c2_size * sizeof(double)); CHECK_ALLOC(c2_cm); }
         if (d2_size > 0) { d2_cm = (double*)malloc(d2_size * sizeof(double)); CHECK_ALLOC(d2_cm); }

         /* Allocate memory for column-major copies of outputs */
         size_t a_size = (size_t)n * n; if (n == 0) a_size = 0;
         size_t b_size = (size_t)n * m1; if (n == 0 || m1 == 0) b_size = 0; // N rows, M1 cols
         size_t c_size = (size_t)p1 * n; if (p1 == 0 || n == 0) c_size = 0; // P1 rows, N cols
         size_t d_size = (size_t)p1 * m1; if (p1 == 0 || m1 == 0) d_size = 0; // P1 rows, M1 cols

         if (a_size > 0) { a_cm = (double*)malloc(a_size * sizeof(double)); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * sizeof(double)); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * sizeof(double)); CHECK_ALLOC(c_cm); }
         if (d_size > 0) { d_cm = (double*)malloc(d_size * sizeof(double)); CHECK_ALLOC(d_cm); }

         /* Transpose C (row-major) inputs to Fortran (column-major) copies */
         if (a1_size > 0) slicot_transpose_to_fortran(a1, a1_cm, n1, n1, sizeof(double));
         if (b1_size > 0) slicot_transpose_to_fortran(b1, b1_cm, n1, m1, sizeof(double));
         if (c1_size > 0) slicot_transpose_to_fortran(c1, c1_cm, p1, n1, sizeof(double)); // P1 rows, N1 cols
         if (d1_size > 0) slicot_transpose_to_fortran(d1, d1_cm, p1, m1, sizeof(double)); // P1 rows, M1 cols
         if (a2_size > 0) slicot_transpose_to_fortran(a2, a2_cm, n2, n2, sizeof(double));
         if (b2_size > 0) slicot_transpose_to_fortran(b2, b2_cm, n2, p1, sizeof(double)); // N2 rows, P1 cols
         if (c2_size > 0) slicot_transpose_to_fortran(c2, c2_cm, m1, n2, sizeof(double)); // M1 rows, N2 cols
         if (d2_size > 0) slicot_transpose_to_fortran(d2, d2_cm, m1, p1, sizeof(double)); // M1 rows, P1 cols

         /* Fortran leading dimensions (number of rows) */
         lda1_f = MAX(1, n1); ldb1_f = MAX(1, n1);
         ldc1_f = MAX(1, p1); ldd1_f = MAX(1, p1);
         lda2_f = MAX(1, n2); ldb2_f = MAX(1, n2);
         ldc2_f = MAX(1, m1); ldd2_f = MAX(1, m1);
         lda_f  = MAX(1, n);  ldb_f  = MAX(1, n);
         ldc_f  = MAX(1, p1); ldd_f  = MAX(1, p1);

         /* Set pointers for Fortran call */
         a1_ptr = a1_cm; b1_ptr = b1_cm; c1_ptr = c1_cm; d1_ptr = d1_cm;
         a2_ptr = a2_cm; b2_ptr = b2_cm; c2_ptr = c2_cm; d2_ptr = d2_cm;
         a_ptr = a_cm; b_ptr = b_cm; c_ptr = c_cm; d_ptr = d_cm;

     } else {
         /* --- Column-Major Case --- */
         lda1_f = lda1; ldb1_f = ldb1; ldc1_f = ldc1; ldd1_f = ldd1;
         lda2_f = lda2; ldb2_f = ldb2; ldc2_f = ldc2; ldd2_f = ldd2;
         lda_f = lda; ldb_f = ldb; ldc_f = ldc; ldd_f = ldd;

         a1_ptr = a1; b1_ptr = b1; c1_ptr = c1; d1_ptr = d1;
         a2_ptr = a2; b2_ptr = b2; c2_ptr = c2; d2_ptr = d2;
         a_ptr = a; b_ptr = b; c_ptr = c; d_ptr = d;
     }

     /* Call the Fortran routine */
     F77_FUNC(ab05nd, AB05ND)(&over_upper,
                              &n1, &m1, &p1, &n2, &alpha,
                              a1_ptr, &lda1_f, b1_ptr, &ldb1_f,
                              c1_ptr, &ldc1_f, d1_ptr, &ldd1_f,
                              a2_ptr, &lda2_f, b2_ptr, &ldb2_f,
                              c2_ptr, &ldc2_f, d2_ptr, &ldd2_f,
                              &n, // Pass calculated N
                              a_ptr, &lda_f, b_ptr, &ldb_f,
                              c_ptr, &ldc_f, d_ptr, &ldd_f,
                              iwork, dwork, &ldwork, &info,
                              over_len);

     /* Copy back results if row_major and successful */
     if (row_major && info == 0) {
         size_t a_size = (size_t)n * n; if (n == 0) a_size = 0;
         size_t b_size = (size_t)n * m1; if (n == 0 || m1 == 0) b_size = 0; // N rows, M1 cols
         size_t c_size = (size_t)p1 * n; if (p1 == 0 || n == 0) c_size = 0; // P1 rows, N cols
         size_t d_size = (size_t)p1 * m1; if (p1 == 0 || m1 == 0) d_size = 0; // P1 rows, M1 cols

         if (a_size > 0) slicot_transpose_to_c(a_cm, a, n, n, sizeof(double));
         if (b_size > 0) slicot_transpose_to_c(b_cm, b, n, m1, sizeof(double)); // N rows, M1 cols
         if (c_size > 0) slicot_transpose_to_c(c_cm, c, p1, n, sizeof(double)); // P1 rows, N cols
         if (d_size > 0) slicot_transpose_to_c(d_cm, d, p1, m1, sizeof(double)); // P1 rows, M1 cols
     }

 cleanup:
     /* --- Cleanup --- */
     // Free allocated memory (free(NULL) is safe)
     free(dwork);
     free(iwork);
     // Free column-major copies if they were allocated
     free(a1_cm); free(b1_cm); free(c1_cm); free(d1_cm);
     free(a2_cm); free(b2_cm); free(c2_cm); free(d2_cm);
     free(a_cm); free(b_cm); free(c_cm); free(d_cm);

     /* Return the info code from the Fortran routine or SLICOT_MEMORY_ERROR */
     return info;
 }
