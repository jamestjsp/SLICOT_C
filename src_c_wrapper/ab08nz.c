/**
 * @file ab08nz.c
 * @brief C wrapper implementation for SLICOT routine AB08NZ
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB08NZ,
 * which constructs a regular pencil for a given complex system such that its
 * generalized eigenvalues are invariant zeros of the system.
 * Refactored to align with ab01nd.c structure.
 */

 #include <stdlib.h>
 #include <ctype.h>  // For toupper
 #include <stddef.h> // For size_t
 #include <complex.h> // For creal

 // Include the header file for this wrapper
 #include "ab08nz.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, MIN, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines, slicot_complex_double, SLICOT_COMPLEX_REAL
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * Note the use of slicot_complex_double for COMPLEX*16 arguments.
  */
 extern void F77_FUNC(ab08nz, AB08NZ)(
     const char* equil,       // CHARACTER*1 EQUIL
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     const int* p,           // INTEGER P
     slicot_complex_double* a, // COMPLEX*16 A(LDA,*) (in/out)
     const int* lda,         // INTEGER LDA
     slicot_complex_double* b, // COMPLEX*16 B(LDB,*) (in/out)
     const int* ldb,         // INTEGER LDB
     slicot_complex_double* c, // COMPLEX*16 C(LDC,*) (in/out)
     const int* ldc,         // INTEGER LDC
     slicot_complex_double* d, // COMPLEX*16 D(LDD,*) (in/out)
     const int* ldd,         // INTEGER LDD
     int* nu,                // INTEGER NU (output)
     int* rank,              // INTEGER RANK (output)
     int* dinfz,             // INTEGER DINFZ (output)
     int* nkror,             // INTEGER NKROR (output)
     int* nkrol,             // INTEGER NKROL (output)
     int* infz,              // INTEGER INFZ(*) (output)
     int* kronr,             // INTEGER KRONR(*) (output)
     int* kronl,             // INTEGER KRONL(*) (output)
     slicot_complex_double* af,// COMPLEX*16 AF(LDAF,*) (output)
     const int* ldaf,        // INTEGER LDAF
     slicot_complex_double* bf,// COMPLEX*16 BF(LDBF,*) (output)
     const int* ldbf,        // INTEGER LDBF
     const double* tol,      // DOUBLE PRECISION TOL
     int* iwork,             // INTEGER IWORK(*)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     slicot_complex_double* zwork,// COMPLEX*16 ZWORK(*)
     const int* lzwork,      // INTEGER LZWORK
     int* info,              // INTEGER INFO (output)
     int equil_len           // Hidden length for equil
 );


 /* C wrapper function definition */
 SLICOT_C_WRAPPER_API
 int slicot_ab08nz(char equil, int n, int m, int p,
                   slicot_complex_double* a, int lda,
                   slicot_complex_double* b, int ldb,
                   slicot_complex_double* c, int ldc,
                   slicot_complex_double* d, int ldd,
                   int* nu, int* rank,
                   int* dinfz, int* nkror, int* nkrol,
                   int* infz, int* kronr, int* kronl,
                   slicot_complex_double* af, int ldaf,
                   slicot_complex_double* bf, int ldbf,
                   double tol, int row_major)
 {
     /* Local variables */
     int info = 0;
     int lzwork = -1; /* Use -1 for workspace query */
     slicot_complex_double zwork_query; // To receive optimal size
     slicot_complex_double* zwork = NULL; // Complex workspace
     double* dwork = NULL; // Real workspace
     int dwork_size = 0;
     int* iwork = NULL; // Integer workspace
     int iwork_size = 0;
     const int equil_len = 1;

     char equil_upper = toupper(equil);

     /* Pointers for column-major copies if needed */
     slicot_complex_double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL;
     slicot_complex_double *af_cm = NULL, *bf_cm = NULL;

     /* Pointers to pass to Fortran */
     slicot_complex_double *a_ptr, *b_ptr, *c_ptr, *d_ptr;
     slicot_complex_double *af_ptr, *bf_ptr;
     int lda_f, ldb_f, ldc_f, ldd_f;
     int ldaf_f, ldbf_f;

     /* --- Input Parameter Validation --- */

     if (n < 0) { info = -2; goto cleanup; }
     if (m < 0) { info = -3; goto cleanup; }
     if (p < 0) { info = -4; goto cleanup; }
     if (equil_upper != 'S' && equil_upper != 'N') { info = -1; goto cleanup; }
     // Optional: Check TOL range if necessary
     // if (tol < 0.0) { info = -23; goto cleanup; }

     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, n);
     int min_ldb_f = MAX(1, n);
     int min_ldc_f = MAX(1, p);
     int min_ldd_f = MAX(1, p);
     int min_ldaf_f = MAX(1, n + m);
     int min_ldbf_f = MAX(1, n + p);

     if (row_major) {
         // For row-major C, LDA is the number of columns
         int min_lda_rm_cols = n;
         int min_ldb_rm_cols = m;
         int min_ldc_rm_cols = n;
         int min_ldd_rm_cols = m;
         // LDAF/LDBF check depends on NU, done later
         if (n > 0 && lda < min_lda_rm_cols) { info = -6; goto cleanup; }
         if (n > 0 && ldb < min_ldb_rm_cols) { info = -8; goto cleanup; }
         if (p > 0 && ldc < min_ldc_rm_cols) { info = -10; goto cleanup; }
         if (p > 0 && ldd < min_ldd_rm_cols) { info = -12; goto cleanup; }
         // Check LDAF/LDBF against minimum Fortran requirements (rows)
         if (ldaf < min_ldaf_f) { info = -20; goto cleanup; }
         if (ldbf < min_ldbf_f) { info = -22; goto cleanup; }
     } else {
         // For column-major C, LDA is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -6; goto cleanup; }
         if (ldb < min_ldb_f) { info = -8; goto cleanup; }
         if (ldc < min_ldc_f) { info = -10; goto cleanup; }
         if (ldd < min_ldd_f) { info = -12; goto cleanup; }
         if (ldaf < min_ldaf_f) { info = -20; goto cleanup; }
         if (ldbf < min_ldbf_f) { info = -22; goto cleanup; }
     }

     /* --- Workspace Allocation --- */

     // Allocate IWORK (size MAX(M,P))
     iwork_size = MAX(1, MAX(m, p));
     iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
     CHECK_ALLOC(iwork);

     // Allocate DWORK (size MAX(N, 2*MAX(P,M)))
     dwork_size = MAX(1, MAX(n, 2 * MAX(p, m)));
     dwork = (double*)malloc((size_t)dwork_size * sizeof(double));
     CHECK_ALLOC(dwork);

     // Query and allocate ZWORK (Complex workspace)
     lzwork = -1; // Query mode
     // Use dummy LDs for query if dimensions are 0
     int lda_q = row_major ? MAX(1, n) : lda;
     int ldb_q = row_major ? MAX(1, n) : ldb;
     int ldc_q = row_major ? MAX(1, p) : ldc;
     int ldd_q = row_major ? MAX(1, p) : ldd;
     int ldaf_q = ldaf; // Pass user-provided LDs
     int ldbf_q = ldbf;

     F77_FUNC(ab08nz, AB08NZ)(&equil_upper, &n, &m, &p,
                              NULL, &lda_q, NULL, &ldb_q, NULL, &ldc_q, NULL, &ldd_q,
                              nu, rank, dinfz, nkror, nkrol,
                              NULL, NULL, NULL, // NULL integer arrays
                              NULL, &ldaf_q, NULL, &ldbf_q, // NULL output arrays
                              &tol, iwork, dwork,
                              &zwork_query, &lzwork, &info, // Pass address for query result
                              equil_len);

     if (info != 0) { goto cleanup; } // Query failed

     lzwork = (int)SLICOT_COMPLEX_REAL(zwork_query); // Extract real part for size
     // Check against minimum size formula if needed
     int min_lzwork = 1; // Basic minimum
     // Example parsing (needs refinement for full formula)
     min_lzwork = MAX(min_lzwork, MIN(p, m) + MAX(3 * m - 1, n));
     min_lzwork = MAX(min_lzwork, MIN(p, n) + MAX(3 * p - 1, MAX(n + p, n + m)));
     min_lzwork = MAX(min_lzwork, MIN(m, n) + MAX(3 * m - 1, n + m));
     lzwork = MAX(lzwork, min_lzwork);

     zwork = (slicot_complex_double*)malloc((size_t)lzwork * sizeof(slicot_complex_double));
     CHECK_ALLOC(zwork);

     /* --- Prepare Arrays and Call Fortran Routine --- */
     size_t elem_size = sizeof(slicot_complex_double);

     if (row_major) {
         /* --- Row-Major Case --- */

         /* Allocate memory for column-major copies */
         size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
         size_t b_rows = n; size_t b_cols = m; size_t b_size = b_rows * b_cols;
         size_t c_rows = p; size_t c_cols = n; size_t c_size = c_rows * c_cols;
         size_t d_rows = p; size_t d_cols = m; size_t d_size = d_rows * d_cols;
         // Fortran AF is LDAF x (N+MIN(P,M)), BF is LDBF x (N+M)
         size_t af_fort_rows = ldaf; size_t af_fort_cols = n + MIN(p, m); size_t af_size = af_fort_rows * af_fort_cols;
         size_t bf_fort_rows = ldbf; size_t bf_fort_cols = n + m; size_t bf_size = bf_fort_rows * bf_fort_cols;

         if (a_size > 0) { a_cm = (slicot_complex_double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (slicot_complex_double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (slicot_complex_double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
         if (d_size > 0) { d_cm = (slicot_complex_double*)malloc(d_size * elem_size); CHECK_ALLOC(d_cm); }
         if (af_size > 0) { af_cm = (slicot_complex_double*)malloc(af_size * elem_size); CHECK_ALLOC(af_cm); }
         if (bf_size > 0) { bf_cm = (slicot_complex_double*)malloc(bf_size * elem_size); CHECK_ALLOC(bf_cm); }

         /* Transpose row-major inputs to column-major copies */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         if (b_size > 0) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, elem_size);
         if (c_size > 0) slicot_transpose_to_fortran(c, c_cm, c_rows, c_cols, elem_size);
         if (d_size > 0) slicot_transpose_to_fortran(d, d_cm, d_rows, d_cols, elem_size);

         /* Fortran leading dimensions */
         lda_f = MAX(1, a_rows);
         ldb_f = MAX(1, b_rows);
         ldc_f = MAX(1, c_rows);
         ldd_f = MAX(1, d_rows);
         ldaf_f = ldaf; // Use passed C ldaf as Fortran LDAF
         ldbf_f = ldbf; // Use passed C ldbf as Fortran LDBF

         /* Set pointers for Fortran call */
         a_ptr = a_cm; b_ptr = b_cm; c_ptr = c_cm; d_ptr = d_cm;
         af_ptr = af_cm; bf_ptr = bf_cm;

     } else {
         /* --- Column-Major Case --- */
         lda_f = lda; ldb_f = ldb; ldc_f = ldc; ldd_f = ldd;
         ldaf_f = ldaf; ldbf_f = ldbf;
         a_ptr = a; b_ptr = b; c_ptr = c; d_ptr = d;
         af_ptr = af; bf_ptr = bf;
     }

     /* Call the computational routine */
     F77_FUNC(ab08nz, AB08NZ)(&equil_upper, &n, &m, &p,
                              a_ptr, &lda_f, b_ptr, &ldb_f, c_ptr, &ldc_f, d_ptr, &ldd_f,
                              nu, rank, dinfz, nkror, nkrol,
                              infz, kronr, kronl,
                              af_ptr, &ldaf_f, bf_ptr, &ldbf_f,
                              &tol, iwork, dwork, zwork, &lzwork, &info,
                              equil_len);

     /* Copy back results if row_major */
     if (row_major && info == 0) {
          int nu_val = *nu;
          // Check target C array dimensions (ldaf/ldbf = #cols >= nu)
          if (ldaf < nu_val) { info = -20; goto cleanup; } // Invalid LDAF for computed NU
          if (ldbf < nu_val) { info = -22; goto cleanup; } // Invalid LDBF for computed NU

          // Copy back AF and BF
          if (nu_val > 0) {
              // Fortran AF is LDAF x NU, BF is LDBF x NU
              size_t af_fort_rows = ldaf; size_t af_fort_cols = nu_val;
              size_t bf_fort_rows = ldbf; size_t bf_fort_cols = nu_val;
              if (af_fort_rows * af_fort_cols > 0) {
                  slicot_transpose_to_c(af_cm, af, af_fort_rows, af_fort_cols, elem_size);
              }
              if (bf_fort_rows * bf_fort_cols > 0) {
                  slicot_transpose_to_c(bf_cm, bf, bf_fort_rows, bf_fort_cols, elem_size);
              }
          }
          /* Update original input matrices if changed by the routine (EQUIL='S') */
          if (equil_upper == 'S') {
              size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
              size_t b_rows = n; size_t b_cols = m; size_t b_size = b_rows * b_cols;
              size_t c_rows = p; size_t c_cols = n; size_t c_size = c_rows * c_cols;
              size_t d_rows = p; size_t d_cols = m; size_t d_size = d_rows * d_cols;

              if (a_size > 0) slicot_transpose_to_c(a_cm, a, a_rows, a_cols, elem_size);
              if (b_size > 0) slicot_transpose_to_c(b_cm, b, b_rows, b_cols, elem_size);
              if (c_size > 0) slicot_transpose_to_c(c_cm, c, c_rows, c_cols, elem_size);
              if (d_size > 0) slicot_transpose_to_c(d_cm, d, d_rows, d_cols, elem_size);
          }
     }

 cleanup:
     /* --- Cleanup --- */
     free(zwork);
     free(dwork);
     free(iwork);
     free(a_cm);
     free(b_cm);
     free(c_cm);
     free(d_cm);
     free(af_cm);
     free(bf_cm);

     return info;
 }
