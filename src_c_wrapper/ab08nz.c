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
 // Ensure slicot_utils.h provides: MAX, MIN, CHECK_ALLOC, SLICOT_MEMORY_ERROR,
 // slicot_transpose_to_fortran, slicot_transpose_to_c, slicot_complex_double, SLICOT_COMPLEX_REAL
 #include "slicot_utils.h"
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
 SLICOT_EXPORT
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

     /* Check leading dimensions based on storage order */
     if (row_major) {
         /* For row-major C, make sure arrays have enough columns */
         if (n > 0 && lda < n) { info = -6; goto cleanup; } // A: n rows x n cols
         if (n > 0 && ldb < m) { info = -8; goto cleanup; } // B: n rows x m cols
         if (p > 0 && ldc < n) { info = -10; goto cleanup; } // C: p rows x n cols
         if (p > 0 && ldd < m) { info = -12; goto cleanup; } // D: p rows x m cols
         
         /* For output arrays, C dimensions need to be sufficient for Fortran result */
         if (ldaf < MIN(n + m, n + MIN(p, m))) { info = -20; goto cleanup; } // Need enough columns
         if (ldbf < n + m) { info = -22; goto cleanup; } // Need enough columns
     } else {
         /* For column-major (Fortran style) */
         if (lda < MAX(1, n)) { info = -6; goto cleanup; } // A: need >= n rows
         if (ldb < MAX(1, n)) { info = -8; goto cleanup; } // B: need >= n rows
         if (ldc < MAX(1, p)) { info = -10; goto cleanup; } // C: need >= p rows
         if (ldd < MAX(1, p)) { info = -12; goto cleanup; } // D: need >= p rows
         if (ldaf < MAX(1, n + m)) { info = -20; goto cleanup; } // AF: need >= n+m rows
         if (ldbf < MAX(1, n + p)) { info = -22; goto cleanup; } // BF: need >= n+p rows
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
     int lda_q = MAX(1, n);  // Ensure valid LDA for query - needs to be >= n
     int ldb_q = MAX(1, n);  // Ensure valid LDB for query - needs to be >= n
     int ldc_q = MAX(1, p);  // Ensure valid LDC for query - needs to be >= p
     int ldd_q = MAX(1, p);  // Ensure valid LDD for query - needs to be >= p
     int ldaf_q = MAX(1, n + m); // Ensure valid LDAF for query
     int ldbf_q = MAX(1, n + p); // Ensure valid LDBF for query

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
     // Check against minimum size formula if needed (based on AB08MZ doc)
     int min_lzwork = 1;
     if (n > 0 || m > 0 || p > 0) {
        min_lzwork = (n+p)*(n+m) + MAX(1, MIN(p,m) + MAX(3*m-1,n));
        min_lzwork = MAX(min_lzwork, (n+p)*(n+m) + MIN(p,n) + MAX(3*p-1, MAX(n+p, n+m)));
     }
     lzwork = MAX(lzwork, min_lzwork);

     zwork = (slicot_complex_double*)malloc((size_t)lzwork * sizeof(slicot_complex_double));
     CHECK_ALLOC(zwork);

     /* --- Prepare Arrays and Call Fortran Routine --- */
     size_t elem_size = sizeof(slicot_complex_double);

     if (row_major) {
         /* --- Row-Major Case --- */

         /* Allocate memory for column-major copies */
         size_t a_size = n * n;     // A: NxN matrix
         size_t b_size = n * m;     // B: NxM matrix
         size_t c_size = p * n;     // C: PxN matrix
         size_t d_size = p * m;     // D: PxM matrix
         
         /* Fortran dimensions for AF and BF */
         ldaf_f = MAX(1, n + m);    // AF needs at least N+M rows in Fortran
         ldbf_f = MAX(1, n + p);    // BF needs at least N+P rows in Fortran
         
         /* Determine sizes for Fortran arrays */
         size_t af_size = ldaf_f * (n + MIN(p, m));  // AF is LDAF_f x (N+MIN(P,M))
         size_t bf_size = ldbf_f * (n + m);          // BF is LDBF_f x (N+M)

         /* Allocate memory for input array copies */
         if (a_size > 0) { a_cm = (slicot_complex_double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (slicot_complex_double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (slicot_complex_double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
         if (d_size > 0) { d_cm = (slicot_complex_double*)malloc(d_size * elem_size); CHECK_ALLOC(d_cm); }
         
         /* Allocate memory for output arrays */
         if (af_size > 0) { af_cm = (slicot_complex_double*)malloc(af_size * elem_size); CHECK_ALLOC(af_cm); }
         if (bf_size > 0) { bf_cm = (slicot_complex_double*)malloc(bf_size * elem_size); CHECK_ALLOC(bf_cm); }

         /* Transpose row-major inputs to column-major copies */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, n, n, elem_size);
         if (b_size > 0) slicot_transpose_to_fortran(b, b_cm, n, m, elem_size);
         if (c_size > 0) slicot_transpose_to_fortran(c, c_cm, p, n, elem_size);
         if (d_size > 0) slicot_transpose_to_fortran(d, d_cm, p, m, elem_size);

         /* Fortran leading dimensions */
         lda_f = MAX(1, n);    // A: needs N rows for Fortran
         ldb_f = MAX(1, n);    // B: needs N rows for Fortran
         ldc_f = MAX(1, p);    // C: needs P rows for Fortran
         ldd_f = MAX(1, p);    // D: needs P rows for Fortran
         /* ldaf_f and ldbf_f already set above */

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

          // Copy back AF and BF (only the NU x NU part)
          if (nu_val > 0) {
              // Use slicot_transpose_to_c_with_ld to properly handle leading dimensions
              if (af_cm) { // Check pointer validity
                   slicot_transpose_to_c_with_ld(af_cm, af, nu_val, nu_val, 
                                                ldaf_f, ldaf, elem_size);
              }
              if (bf_cm) { // Check pointer validity
                   slicot_transpose_to_c_with_ld(bf_cm, bf, nu_val, nu_val, 
                                                ldbf_f, ldbf, elem_size);
              }
          }
          
          /* Update original input matrices if changed by the routine (EQUIL='S') */
          if (equil_upper == 'S') {
              // Dimensions of the original C matrices
              size_t a_size = n * n;
              size_t b_size = n * m;
              size_t c_size = p * n;
              size_t d_size = p * m;

              // Transpose back from temporary column-major storage
              if (a_cm && a_size > 0) slicot_transpose_to_c(a_cm, a, n, n, elem_size);
              if (b_cm && b_size > 0) slicot_transpose_to_c(b_cm, b, n, m, elem_size);
              if (c_cm && c_size > 0) slicot_transpose_to_c(c_cm, c, p, n, elem_size);
              if (d_cm && d_size > 0) slicot_transpose_to_c(d_cm, d, p, m, elem_size);
          }
     }

 cleanup:
     /* --- Cleanup --- */
     free(zwork);
     free(dwork);
     free(iwork);
     // Free temporary column-major arrays if allocated
     free(a_cm);
     free(b_cm);
     free(c_cm);
     free(d_cm);
     free(af_cm);
     free(bf_cm);

     return info;
 }
