/**
 * @file ab05od.c
 * @brief C wrapper implementation for SLICOT routine AB05OD
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB05OD,
 * which computes the state-space model (A,B,C,D) for the rowwise
 * concatenation (parallel inter-connection with separate inputs)
 * of two systems.
 * This wrapper assumes OVER = 'N' (no overlapping arrays).
 */

 #include <stdlib.h>
 #include <ctype.h>  // For toupper
 #include <stddef.h> // For size_t
 
 // Include the header file for this wrapper
 #include "ab05od.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions
 
 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * Use const for input arrays as this wrapper assumes OVER = 'N'.
  * Hidden string length arguments are explicitly included.
  */
 extern void F77_FUNC(ab05od, AB05OD)(
     const char* over,       // CHARACTER*1 OVER
     const int* n1,          // INTEGER N1
     const int* m1,          // INTEGER M1
     const int* p1,          // INTEGER P1
     const int* n2,          // INTEGER N2
     const int* m2,          // INTEGER M2
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
     int* m,                 // INTEGER M (output)
     double* a,              // DOUBLE PRECISION A(LDA,*) (output)
     const int* lda,         // INTEGER LDA
     double* b,              // DOUBLE PRECISION B(LDB,*) (output)
     const int* ldb,         // INTEGER LDB
     double* c,              // DOUBLE PRECISION C(LDC,*) (output)
     const int* ldc,         // INTEGER LDC
     double* d,              // DOUBLE PRECISION D(LDD,*) (output)
     const int* ldd,         // INTEGER LDD
     int* info,              // INTEGER INFO (output)
     int over_len            // Hidden length for over
 );
 
 
 /* C wrapper function definition */
 int slicot_ab05od(char over,
                   int n1, int m1, int p1, int n2, int m2, double alpha,
                   const double* a1, int lda1, const double* b1, int ldb1,
                   const double* c1, int ldc1, const double* d1, int ldd1,
                   const double* a2, int lda2, const double* b2, int ldb2,
                   const double* c2, int ldc2, const double* d2, int ldd2,
                   int* n, int* m,
                   double* a, int lda, double* b, int ldb,
                   double* c, int ldc, double* d, int ldd,
                   int row_major)
 {
     /* Local variables */
     int info = 0;
     const int over_len = 1;
     int n_calc = 0;       // Calculated total state dimension
     int m_calc = 0;       // Calculated total input dimension
 
     char over_upper = toupper(over); // Although we assume 'N'
 
     /* Pointers for column-major copies if needed */
     double *a1_cm = NULL, *b1_cm = NULL, *c1_cm = NULL, *d1_cm = NULL;
     double *a2_cm = NULL, *b2_cm = NULL, *c2_cm = NULL, *d2_cm = NULL;
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL;
 
     /* --- Input Parameter Validation --- */
 
     if (n1 < 0) { info = -2; goto cleanup; }
     if (m1 < 0) { info = -3; goto cleanup; }
     if (p1 < 0) { info = -4; goto cleanup; }
     if (n2 < 0) { info = -5; goto cleanup; }
     if (m2 < 0) { info = -6; goto cleanup; }
     // No check for ALPHA needed based on docs
 
     // Force OVER = 'N' as 'O' is not supported by this wrapper interface
     if (over_upper == 'O') {
         // Optionally print a warning or return an error
         // For now, proceed assuming 'N' was intended for the wrapper call
         over_upper = 'N';
     }
     if (over_upper != 'N') { info = -1; goto cleanup; } // Should not happen now
 
     // Calculate resulting dimensions N and M
     n_calc = n1 + n2;
     m_calc = m1 + m2;
     if (n) { *n = n_calc; } // Set output N
     if (m) { *m = m_calc; } // Set output M
 
     // Check leading dimensions for inputs
     int min_lda1_f = MAX(1, n1); int min_ldb1_f = MAX(1, n1);
     int min_ldc1_f = (n1 > 0) ? MAX(1, p1) : 1; int min_ldd1_f = MAX(1, p1);
     int min_lda2_f = MAX(1, n2); int min_ldb2_f = MAX(1, n2);
     int min_ldc2_f = (n2 > 0) ? MAX(1, p1) : 1; int min_ldd2_f = MAX(1, p1); // Rows=P1
     // Check leading dimensions for outputs
     int min_lda_f = MAX(1, n_calc); int min_ldb_f = MAX(1, n_calc);
     int min_ldc_f = (n_calc > 0) ? MAX(1, p1) : 1; int min_ldd_f = MAX(1, p1); // Rows=P1
 
     if (row_major) {
         // Check C dimensions (number of columns)
         if (lda1 < n1) { info = -9; goto cleanup; }
         if (ldb1 < m1) { info = -11; goto cleanup; }
         if (ldc1 < n1) { info = -13; goto cleanup; }
         if (ldd1 < m1) { info = -15; goto cleanup; }
         if (lda2 < n2) { info = -17; goto cleanup; }
         if (ldb2 < m2) { info = -19; goto cleanup; }
         if (ldc2 < n2) { info = -21; goto cleanup; }
         if (ldd2 < m2) { info = -23; goto cleanup; }
         // Check output C dimensions
         if (lda < n_calc) { info = -26; goto cleanup; }
         if (ldb < m_calc) { info = -28; goto cleanup; } // Cols=M
         if (ldc < n_calc) { info = -30; goto cleanup; }
         if (ldd < m_calc) { info = -32; goto cleanup; } // Cols=M
     } else {
         // Check Fortran dimensions (number of rows)
         if (lda1 < min_lda1_f) { info = -9; goto cleanup; }
         if (ldb1 < min_ldb1_f) { info = -11; goto cleanup; }
         if (ldc1 < min_ldc1_f) { info = -13; goto cleanup; }
         if (ldd1 < min_ldd1_f) { info = -15; goto cleanup; }
         if (lda2 < min_lda2_f) { info = -17; goto cleanup; }
         if (ldb2 < min_ldb2_f) { info = -19; goto cleanup; }
         if (ldc2 < min_ldc2_f) { info = -21; goto cleanup; }
         if (ldd2 < min_ldd2_f) { info = -23; goto cleanup; }
         // Check output Fortran dimensions
         if (lda < min_lda_f) { info = -26; goto cleanup; }
         if (ldb < min_ldb_f) { info = -28; goto cleanup; }
         if (ldc < min_ldc_f) { info = -30; goto cleanup; }
         if (ldd < min_ldd_f) { info = -32; goto cleanup; }
     }
 
     /* --- Workspace Allocation (None needed for AB05OD) --- */
 
     /* --- Prepare Arrays and Call Fortran Routine --- */
 
     if (row_major) {
         /* --- Row-Major Case --- */
 
         /* Allocate memory for column-major copies of inputs */
         size_t a1_size = (size_t)n1 * n1; if (n1 == 0) a1_size = 0;
         size_t b1_size = (size_t)n1 * m1; if (n1 == 0 || m1 == 0) b1_size = 0;
         size_t c1_size = (size_t)p1 * n1; if (p1 == 0 || n1 == 0) c1_size = 0;
         size_t d1_size = (size_t)p1 * m1; if (p1 == 0 || m1 == 0) d1_size = 0;
         size_t a2_size = (size_t)n2 * n2; if (n2 == 0) a2_size = 0;
         size_t b2_size = (size_t)n2 * m2; if (n2 == 0 || m2 == 0) b2_size = 0;
         size_t c2_size = (size_t)p1 * n2; if (p1 == 0 || n2 == 0) c2_size = 0; // Rows=P1
         size_t d2_size = (size_t)p1 * m2; if (p1 == 0 || m2 == 0) d2_size = 0; // Rows=P1
 
         if (a1_size > 0) { a1_cm = (double*)malloc(a1_size * sizeof(double)); CHECK_ALLOC(a1_cm); }
         if (b1_size > 0) { b1_cm = (double*)malloc(b1_size * sizeof(double)); CHECK_ALLOC(b1_cm); }
         if (c1_size > 0) { c1_cm = (double*)malloc(c1_size * sizeof(double)); CHECK_ALLOC(c1_cm); }
         if (d1_size > 0) { d1_cm = (double*)malloc(d1_size * sizeof(double)); CHECK_ALLOC(d1_cm); }
         if (a2_size > 0) { a2_cm = (double*)malloc(a2_size * sizeof(double)); CHECK_ALLOC(a2_cm); }
         if (b2_size > 0) { b2_cm = (double*)malloc(b2_size * sizeof(double)); CHECK_ALLOC(b2_cm); }
         if (c2_size > 0) { c2_cm = (double*)malloc(c2_size * sizeof(double)); CHECK_ALLOC(c2_cm); }
         if (d2_size > 0) { d2_cm = (double*)malloc(d2_size * sizeof(double)); CHECK_ALLOC(d2_cm); }
 
         /* Allocate memory for column-major copies of outputs */
         size_t a_size = (size_t)n_calc * n_calc; if (n_calc == 0) a_size = 0;
         size_t b_size = (size_t)n_calc * m_calc; if (n_calc == 0 || m_calc == 0) b_size = 0;
         size_t c_size = (size_t)p1 * n_calc; if (p1 == 0 || n_calc == 0) c_size = 0; // Rows=P1
         size_t d_size = (size_t)p1 * m_calc; if (p1 == 0 || m_calc == 0) d_size = 0; // Rows=P1
 
         if (a_size > 0) { a_cm = (double*)malloc(a_size * sizeof(double)); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * sizeof(double)); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * sizeof(double)); CHECK_ALLOC(c_cm); }
         if (d_size > 0) { d_cm = (double*)malloc(d_size * sizeof(double)); CHECK_ALLOC(d_cm); }
 
         /* Transpose C (row-major) inputs to Fortran (column-major) copies */
         if (a1_size > 0) slicot_transpose_to_fortran(a1, a1_cm, n1, n1, sizeof(double));
         if (b1_size > 0) slicot_transpose_to_fortran(b1, b1_cm, n1, m1, sizeof(double));
         if (c1_size > 0) slicot_transpose_to_fortran(c1, c1_cm, p1, n1, sizeof(double));
         if (d1_size > 0) slicot_transpose_to_fortran(d1, d1_cm, p1, m1, sizeof(double));
         if (a2_size > 0) slicot_transpose_to_fortran(a2, a2_cm, n2, n2, sizeof(double));
         if (b2_size > 0) slicot_transpose_to_fortran(b2, b2_cm, n2, m2, sizeof(double));
         if (c2_size > 0) slicot_transpose_to_fortran(c2, c2_cm, p1, n2, sizeof(double)); // Rows=P1
         if (d2_size > 0) slicot_transpose_to_fortran(d2, d2_cm, p1, m2, sizeof(double)); // Rows=P1
 
         /* Fortran leading dimensions */
         int lda1_f = (n1 > 0) ? n1 : 1; int ldb1_f = (n1 > 0) ? n1 : 1;
         int ldc1_f = (p1 > 0) ? p1 : 1; int ldd1_f = (p1 > 0) ? p1 : 1;
         int lda2_f = (n2 > 0) ? n2 : 1; int ldb2_f = (n2 > 0) ? n2 : 1;
         int ldc2_f = (p1 > 0) ? p1 : 1; int ldd2_f = (p1 > 0) ? p1 : 1; // Rows=P1
         int lda_f = (n_calc > 0) ? n_calc : 1; int ldb_f = (n_calc > 0) ? n_calc : 1;
         int ldc_f = (p1 > 0) ? p1 : 1; int ldd_f = (p1 > 0) ? p1 : 1; // Rows=P1
 
         /* Call the Fortran routine */
         F77_FUNC(ab05od, AB05OD)(&over_upper,
                                  &n1, &m1, &p1, &n2, &m2, &alpha,
                                  a1_cm, &lda1_f, b1_cm, &ldb1_f,
                                  c1_cm, &ldc1_f, d1_cm, &ldd1_f,
                                  a2_cm, &lda2_f, b2_cm, &ldb2_f,
                                  c2_cm, &ldc2_f, d2_cm, &ldd2_f,
                                  &n_calc, &m_calc, // Pass calculated N, M
                                  a_cm, &lda_f, b_cm, &ldb_f,
                                  c_cm, &ldc_f, d_cm, &ldd_f,
                                  &info,
                                  over_len);
 
         /* Copy back results from column-major temps to original row-major arrays */
         if (info == 0) {
             if (a_size > 0) slicot_transpose_to_c(a_cm, a, n_calc, n_calc, sizeof(double));
             if (b_size > 0) slicot_transpose_to_c(b_cm, b, n_calc, m_calc, sizeof(double));
             if (c_size > 0) slicot_transpose_to_c(c_cm, c, p1, n_calc, sizeof(double)); // Rows=P1
             if (d_size > 0) slicot_transpose_to_c(d_cm, d, p1, m_calc, sizeof(double)); // Rows=P1
         }
         /* Column-major copies will be freed in cleanup */
 
     } else {
         /* --- Column-Major Case --- */
 
         /* Call the Fortran routine directly with user-provided arrays */
         F77_FUNC(ab05od, AB05OD)(&over_upper,
                                  &n1, &m1, &p1, &n2, &m2, &alpha,
                                  a1, &lda1, b1, &ldb1,
                                  c1, &ldc1, d1, &ldd1,
                                  a2, &lda2, b2, &ldb2,
                                  c2, &ldc2, d2, &ldd2,
                                  &n_calc, &m_calc, // Pass calculated N, M
                                  a, &lda, b, &ldb,
                                  c, &ldc, d, &ldd,
                                  &info,
                                  over_len);
         // Output arrays a, b, c, d are filled in place.
     }
 
 cleanup:
     /* --- Cleanup --- */
     // Free allocated memory (free(NULL) is safe)
     // No workspace allocated by this wrapper
     // Free column-major copies if they were allocated
     free(a1_cm); free(b1_cm); free(c1_cm); free(d1_cm);
     free(a2_cm); free(b2_cm); free(c2_cm); free(d2_cm);
     free(a_cm); free(b_cm); free(c_cm); free(d_cm);
 
     /* Return the info code from the Fortran routine or SLICOT_MEMORY_ERROR */
     return info;
 }
 