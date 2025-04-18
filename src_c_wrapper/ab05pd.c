/**
 * @file ab05pd.c
 * @brief C wrapper implementation for SLICOT routine AB05PD
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB05PD,
 * which computes the state-space model (A,B,C,D) for the parallel
 * inter-connection (sum) G = G1 + alpha*G2 of two systems.
 * This wrapper assumes OVER = 'N' (no overlapping arrays).
 */

 #include <stdlib.h>
 #include <ctype.h>  // For toupper
 #include <stddef.h> // For size_t
 
 // Include the header file for this wrapper
 #include "ab05pd.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions
 
 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * Use const for input arrays as this wrapper assumes OVER = 'N'.
  * Hidden string length arguments are explicitly included.
  */
 extern void F77_FUNC(ab05pd, AB05PD)(
     const char* over,       // CHARACTER*1 OVER
     const int* n1,          // INTEGER N1
     const int* m,           // INTEGER M
     const int* p,           // INTEGER P
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
     int* info,              // INTEGER INFO (output)
     int over_len            // Hidden length for over
 );
 
 
 /* C wrapper function definition */
 SLICOT_C_WRAPPER_API
 int slicot_ab05pd(char over,
                   int n1, int m_in, int p_in, int n2, double alpha,
                   const double* a1, int lda1, const double* b1, int ldb1,
                   const double* c1, int ldc1, const double* d1, int ldd1,
                   const double* a2, int lda2, const double* b2, int ldb2,
                   const double* c2, int ldc2, const double* d2, int ldd2,
                   int* n,
                   double* a, int lda, double* b, int ldb,
                   double* c, int ldc, double* d, int ldd,
                   int row_major)
 {
     /* Local variables */
     int info = 0;
     const int over_len = 1;
     int n_calc = 0;       // Calculated total state dimension
 
     char over_upper = toupper(over); // Although we assume 'N'
 
     /* Pointers for column-major copies if needed */
     double *a1_cm = NULL, *b1_cm = NULL, *c1_cm = NULL, *d1_cm = NULL;
     double *a2_cm = NULL, *b2_cm = NULL, *c2_cm = NULL, *d2_cm = NULL;
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL;
 
     /* --- Input Parameter Validation --- */
 
     if (n1 < 0) { info = -2; goto cleanup; }
     if (m_in < 0) { info = -3; goto cleanup; }
     if (p_in < 0) { info = -4; goto cleanup; }
     if (n2 < 0) { info = -5; goto cleanup; }
     // No check for ALPHA needed based on docs
 
     // Force OVER = 'N' as 'O' is not supported by this wrapper interface
     if (over_upper == 'O') {
         // Optionally print a warning or return an error
         // For now, proceed assuming 'N' was intended for the wrapper call
         over_upper = 'N';
     }
     if (over_upper != 'N') { info = -1; goto cleanup; } // Should not happen now
 
     // Calculate resulting dimension N
     n_calc = n1 + n2;
     if (n) { *n = n_calc; } // Set output N
     // M (number of inputs) is an input parameter m_in
 
     // Check leading dimensions for inputs
     int min_lda1_f = MAX(1, n1); int min_ldb1_f = MAX(1, n1);
     int min_ldc1_f = (n1 > 0) ? MAX(1, p_in) : 1; int min_ldd1_f = MAX(1, p_in);
     int min_lda2_f = MAX(1, n2); int min_ldb2_f = MAX(1, n2);
     int min_ldc2_f = (n2 > 0) ? MAX(1, p_in) : 1; int min_ldd2_f = MAX(1, p_in); // Rows=P
     // Check leading dimensions for outputs
     int min_lda_f = MAX(1, n_calc); int min_ldb_f = MAX(1, n_calc);
     int min_ldc_f = (n_calc > 0) ? MAX(1, p_in) : 1; int min_ldd_f = MAX(1, p_in); // Rows=P
 
     if (row_major) {
         // Check C dimensions (number of columns)
         if (lda1 < n1) { info = -8; goto cleanup; }
         if (ldb1 < m_in) { info = -10; goto cleanup; }
         if (ldc1 < n1) { info = -12; goto cleanup; }
         if (ldd1 < m_in) { info = -14; goto cleanup; }
         if (lda2 < n2) { info = -16; goto cleanup; }
         if (ldb2 < m_in) { info = -18; goto cleanup; }
         if (ldc2 < n2) { info = -20; goto cleanup; }
         if (ldd2 < m_in) { info = -22; goto cleanup; }
         // Check output C dimensions
         if (lda < n_calc) { info = -25; goto cleanup; }
         if (ldb < m_in) { info = -27; goto cleanup; } // Cols=M
         if (ldc < n_calc) { info = -29; goto cleanup; }
         if (ldd < m_in) { info = -31; goto cleanup; } // Cols=M
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
 
     /* --- Workspace Allocation (None needed for AB05PD) --- */
 
     /* --- Prepare Arrays and Call Fortran Routine --- */
 
     if (row_major) {
         /* --- Row-Major Case --- */
 
         /* Allocate memory for column-major copies of inputs */
         size_t a1_size = (size_t)n1 * n1; if (n1 == 0) a1_size = 0;
         size_t b1_size = (size_t)n1 * m_in; if (n1 == 0 || m_in == 0) b1_size = 0;
         size_t c1_size = (size_t)p_in * n1; if (p_in == 0 || n1 == 0) c1_size = 0; // Rows=P
         size_t d1_size = (size_t)p_in * m_in; if (p_in == 0 || m_in == 0) d1_size = 0; // Rows=P
         size_t a2_size = (size_t)n2 * n2; if (n2 == 0) a2_size = 0;
         size_t b2_size = (size_t)n2 * m_in; if (n2 == 0 || m_in == 0) b2_size = 0;
         size_t c2_size = (size_t)p_in * n2; if (p_in == 0 || n2 == 0) c2_size = 0; // Rows=P
         size_t d2_size = (size_t)p_in * m_in; if (p_in == 0 || m_in == 0) d2_size = 0; // Rows=P
 
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
         size_t b_size = (size_t)n_calc * m_in; if (n_calc == 0 || m_in == 0) b_size = 0;
         size_t c_size = (size_t)p_in * n_calc; if (p_in == 0 || n_calc == 0) c_size = 0; // Rows=P
         size_t d_size = (size_t)p_in * m_in; if (p_in == 0 || m_in == 0) d_size = 0; // Rows=P
 
         if (a_size > 0) { a_cm = (double*)malloc(a_size * sizeof(double)); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * sizeof(double)); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * sizeof(double)); CHECK_ALLOC(c_cm); }
         if (d_size > 0) { d_cm = (double*)malloc(d_size * sizeof(double)); CHECK_ALLOC(d_cm); }
 
         /* Transpose C (row-major) inputs to Fortran (column-major) copies */
         if (a1_size > 0) slicot_transpose_to_fortran(a1, a1_cm, n1, n1, sizeof(double));
         if (b1_size > 0) slicot_transpose_to_fortran(b1, b1_cm, n1, m_in, sizeof(double));
         if (c1_size > 0) slicot_transpose_to_fortran(c1, c1_cm, p_in, n1, sizeof(double)); // Rows=P
         if (d1_size > 0) slicot_transpose_to_fortran(d1, d1_cm, p_in, m_in, sizeof(double)); // Rows=P
         if (a2_size > 0) slicot_transpose_to_fortran(a2, a2_cm, n2, n2, sizeof(double));
         if (b2_size > 0) slicot_transpose_to_fortran(b2, b2_cm, n2, m_in, sizeof(double));
         if (c2_size > 0) slicot_transpose_to_fortran(c2, c2_cm, p_in, n2, sizeof(double)); // Rows=P
         if (d2_size > 0) slicot_transpose_to_fortran(d2, d2_cm, p_in, m_in, sizeof(double)); // Rows=P
 
         /* Fortran leading dimensions */
         int lda1_f = (n1 > 0) ? n1 : 1; int ldb1_f = (n1 > 0) ? n1 : 1;
         int ldc1_f = (p_in > 0) ? p_in : 1; int ldd1_f = (p_in > 0) ? p_in : 1; // Rows=P
         int lda2_f = (n2 > 0) ? n2 : 1; int ldb2_f = (n2 > 0) ? n2 : 1;
         int ldc2_f = (p_in > 0) ? p_in : 1; int ldd2_f = (p_in > 0) ? p_in : 1; // Rows=P
         int lda_f = (n_calc > 0) ? n_calc : 1; int ldb_f = (n_calc > 0) ? n_calc : 1;
         int ldc_f = (p_in > 0) ? p_in : 1; int ldd_f = (p_in > 0) ? p_in : 1; // Rows=P
 
         /* Call the Fortran routine */
         F77_FUNC(ab05pd, AB05PD)(&over_upper,
                                  &n1, &m_in, &p_in, &n2, &alpha,
                                  a1_cm, &lda1_f, b1_cm, &ldb1_f,
                                  c1_cm, &ldc1_f, d1_cm, &ldd1_f,
                                  a2_cm, &lda2_f, b2_cm, &ldb2_f,
                                  c2_cm, &ldc2_f, d2_cm, &ldd2_f,
                                  &n_calc, // Pass calculated N
                                  a_cm, &lda_f, b_cm, &ldb_f,
                                  c_cm, &ldc_f, d_cm, &ldd_f,
                                  &info,
                                  over_len);
 
         /* Copy back results from column-major temps to original row-major arrays */
         if (info == 0) {
             if (a_size > 0) slicot_transpose_to_c(a_cm, a, n_calc, n_calc, sizeof(double));
             if (b_size > 0) slicot_transpose_to_c(b_cm, b, n_calc, m_in, sizeof(double));
             if (c_size > 0) slicot_transpose_to_c(c_cm, c, p_in, n_calc, sizeof(double)); // Rows=P
             if (d_size > 0) slicot_transpose_to_c(d_cm, d, p_in, m_in, sizeof(double)); // Rows=P
         }
         /* Column-major copies will be freed in cleanup */
 
     } else {
         /* --- Column-Major Case --- */
 
         /* Call the Fortran routine directly with user-provided arrays */
         F77_FUNC(ab05pd, AB05PD)(&over_upper,
                                  &n1, &m_in, &p_in, &n2, &alpha,
                                  a1, &lda1, b1, &ldb1,
                                  c1, &ldc1, d1, &ldd1,
                                  a2, &lda2, b2, &ldb2,
                                  c2, &ldc2, d2, &ldd2,
                                  &n_calc, // Pass calculated N
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
 