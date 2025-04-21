/**
 * @file sb10dd.c
 * @brief C wrapper implementation for SLICOT routine SB10DD
 *
 * This file provides a C wrapper implementation for the SLICOT routine SB10DD,
 * which computes an H-infinity (sub)optimal controller for a
 * discrete-time system.
 */

 #include <stdlib.h>
 #include <stddef.h> // For size_t

 // Include the header file for this wrapper
 #include "sb10dd.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines, slicot_transpose_symmetric_to_c
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * Input matrices are const. Output matrices AK..DK, X, Z are non-const.
  * Output scalar array is RCOND.
  */
 extern void F77_FUNC(sb10dd, SB10DD)(
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     const int* np,          // INTEGER NP
     const int* ncon,        // INTEGER NCON
     const int* nmeas,       // INTEGER NMEAS
     const double* gamma,    // DOUBLE PRECISION GAMMA (in)
     const double* a,        // DOUBLE PRECISION A(LDA,*) (in)
     const int* lda,         // INTEGER LDA
     const double* b,        // DOUBLE PRECISION B(LDB,*) (in)
     const int* ldb,         // INTEGER LDB
     const double* c,        // DOUBLE PRECISION C(LDC,*) (in)
     const int* ldc,         // INTEGER LDC
     const double* d,        // DOUBLE PRECISION D(LDD,*) (in)
     const int* ldd,         // INTEGER LDD
     double* ak,             // DOUBLE PRECISION AK(LDAK,*) (output)
     const int* ldak,        // INTEGER LDAK
     double* bk,             // DOUBLE PRECISION BK(LDBK,*) (output)
     const int* ldbk,        // INTEGER LDBK
     double* ck,             // DOUBLE PRECISION CK(LDCK,*) (output)
     const int* ldck,        // INTEGER LDCK
     double* dk,             // DOUBLE PRECISION DK(LDDK,*) (output)
     const int* lddk,        // INTEGER LDDK
     double* x,              // DOUBLE PRECISION X(LDX,*) (output, symmetric)
     const int* ldx,         // INTEGER LDX
     double* z,              // DOUBLE PRECISION Z(LDZ,*) (output, symmetric)
     const int* ldz,         // INTEGER LDZ
     double* rcond,          // DOUBLE PRECISION RCOND(8) (output)
     const double* tol,      // DOUBLE PRECISION TOL (in)
     int* iwork,             // INTEGER IWORK(*)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* bwork,             // LOGICAL BWORK(*) -> int*
     int* info               // INTEGER INFO (output)
 );


 /* C wrapper function definition */
 SLICOT_C_WRAPPER_API
 int slicot_sb10dd(int n, int m, int np, int ncon, int nmeas,
                   double gamma, const double* a, int lda,
                   const double* b, int ldb, const double* c, int ldc,
                   const double* d, int ldd, double* ak, int ldak,
                   double* bk, int ldbk, double* ck, int ldck,
                   double* dk, int lddk, double* x, int ldx,
                   double* z, int ldz, double* rcond, double tol,
                   int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = -1; /* Use -1 for workspace query */
     int liwork = 0;  /* Fixed size */
     int lbwork = 0;  /* Fixed size */
     double dwork_query;
     double* dwork = NULL;
     int* iwork = NULL;
     int* bwork = NULL; // Map LOGICAL to int

     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL;
     double *ak_cm = NULL, *bk_cm = NULL, *ck_cm = NULL, *dk_cm = NULL;
     double *x_cm = NULL, *z_cm = NULL;
     const double *a_ptr, *b_ptr, *c_ptr, *d_ptr; // Input pointers
     double *ak_ptr, *bk_ptr, *ck_ptr, *dk_ptr; // Output pointers (K)
     double *x_ptr, *z_ptr; // Output pointers (X, Z)
     int lda_f, ldb_f, ldc_f, ldd_f, ldak_f, ldbk_f, ldck_f, lddk_f;
     int ldx_f, ldz_f;


     /* --- Input Parameter Validation --- */

     if (n < 0) { info = -1; goto cleanup; }
     if (m < 0) { info = -2; goto cleanup; }
     if (np < 0) { info = -3; goto cleanup; }
     if (ncon < 0 || ncon > m) { info = -4; goto cleanup; }
     if (nmeas < 0 || nmeas > np) { info = -5; goto cleanup; }
     if (gamma <= 0.0) { info = -6; goto cleanup; }
     // TOL checked by Fortran

     // Calculate partition dimensions needed for validation/workspace
     int m1 = m - ncon;      // Columns of B1, D11, D21
     int np2 = nmeas;        // Rows of C2, D21, D22
     int np1 = np - np2;     // Rows of C1, D11, D12
     int m2 = ncon;          // Columns of B2, D12, D22

     // Further dimension checks related to partitioning
     if (np1 < 0) { info = -5; goto cleanup; } // Implied by nmeas <= np
     if (m1 < 0)  { info = -4; goto cleanup; } // Implied by ncon <= m


     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, n);
     int min_ldb_f = MAX(1, n);
     int min_ldc_f = MAX(1, np);
     int min_ldd_f = MAX(1, np);
     int min_ldak_f = MAX(1, n);
     int min_ldbk_f = MAX(1, n);
     int min_ldck_f = MAX(1, ncon);
     int min_lddk_f = MAX(1, ncon);
     int min_ldx_f = MAX(1, n);
     int min_ldz_f = MAX(1, n);


     if (row_major) {
         // For row-major C, LD is the number of columns
         int min_lda_rm_cols = n;
         int min_ldb_rm_cols = m;
         int min_ldc_rm_cols = n;
         int min_ldd_rm_cols = m;
         int min_ldak_rm_cols = n;
         int min_ldbk_rm_cols = nmeas;
         int min_ldck_rm_cols = n;
         int min_lddk_rm_cols = nmeas;
         int min_ldx_rm_cols = n;
         int min_ldz_rm_cols = n;

         if (lda < min_lda_rm_cols) { info = -8; goto cleanup; }
         if (ldb < min_ldb_rm_cols) { info = -10; goto cleanup; }
         if (ldc < min_ldc_rm_cols) { info = -12; goto cleanup; }
         if (ldd < min_ldd_rm_cols) { info = -14; goto cleanup; }
         if (ldak < min_ldak_rm_cols) { info = -16; goto cleanup; }
         if (ldbk < min_ldbk_rm_cols) { info = -18; goto cleanup; }
         if (ldck < min_ldck_rm_cols) { info = -20; goto cleanup; }
         if (lddk < min_lddk_rm_cols) { info = -22; goto cleanup; }
         if (ldx < min_ldx_rm_cols) { info = -24; goto cleanup; }
         if (ldz < min_ldz_rm_cols) { info = -26; goto cleanup; }
     } else {
         // For column-major C, LD is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -8; goto cleanup; }
         if (ldb < min_ldb_f) { info = -10; goto cleanup; }
         if (ldc < min_ldc_f) { info = -12; goto cleanup; }
         if (ldd < min_ldd_f) { info = -14; goto cleanup; }
         if (ldak < min_ldak_f) { info = -16; goto cleanup; }
         if (ldbk < min_ldbk_f) { info = -18; goto cleanup; }
         if (ldck < min_ldck_f) { info = -20; goto cleanup; }
         if (lddk < min_lddk_f) { info = -22; goto cleanup; }
         if (ldx < min_ldx_f) { info = -24; goto cleanup; }
         if (ldz < min_ldz_f) { info = -26; goto cleanup; }
     }

     /* --- Prepare arrays for column-major format if using row-major --- */
     size_t elem_size = sizeof(double);
     if (row_major) {
         /* Allocate memory for column-major copies */
         size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
         size_t b_rows = n; size_t b_cols = m; size_t b_size = b_rows * b_cols;
         size_t c_rows = np; size_t c_cols = n; size_t c_size = c_rows * c_cols;
         size_t d_rows = np; size_t d_cols = m; size_t d_size = d_rows * d_cols;
         size_t ak_rows = n; size_t ak_cols = n; size_t ak_size = ak_rows * ak_cols;
         size_t bk_rows = n; size_t bk_cols = nmeas; size_t bk_size = bk_rows * bk_cols;
         size_t ck_rows = ncon; size_t ck_cols = n; size_t ck_size = ck_rows * ck_cols;
         size_t dk_rows = ncon; size_t dk_cols = nmeas; size_t dk_size = dk_rows * dk_cols;
         size_t x_rows = n; size_t x_cols = n; size_t x_size = x_rows * x_cols;
         size_t z_rows = n; size_t z_cols = n; size_t z_size = z_rows * z_cols;


         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
         if (d_size > 0) { d_cm = (double*)malloc(d_size * elem_size); CHECK_ALLOC(d_cm); }
         if (ak_size > 0) { ak_cm = (double*)malloc(ak_size * elem_size); CHECK_ALLOC(ak_cm); } // Output
         if (bk_size > 0) { bk_cm = (double*)malloc(bk_size * elem_size); CHECK_ALLOC(bk_cm); } // Output
         if (ck_size > 0) { ck_cm = (double*)malloc(ck_size * elem_size); CHECK_ALLOC(ck_cm); } // Output
         if (dk_size > 0) { dk_cm = (double*)malloc(dk_size * elem_size); CHECK_ALLOC(dk_cm); } // Output
         if (x_size > 0) { x_cm = (double*)malloc(x_size * elem_size); CHECK_ALLOC(x_cm); }   // Output
         if (z_size > 0) { z_cm = (double*)malloc(z_size * elem_size); CHECK_ALLOC(z_cm); }   // Output


         /* Transpose C inputs to Fortran copies */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, n, n, elem_size);
         if (b_size > 0) slicot_transpose_to_fortran(b, b_cm, n, m, elem_size);
         if (c_size > 0) slicot_transpose_to_fortran(c, c_cm, np, n, elem_size);
         if (d_size > 0) slicot_transpose_to_fortran(d, d_cm, np, m, elem_size);


         /* Fortran leading dimensions */
         lda_f = (a_rows > 0) ? a_rows : 1;
         ldb_f = (b_rows > 0) ? b_rows : 1;
         ldc_f = (c_rows > 0) ? c_rows : 1;
         ldd_f = (d_rows > 0) ? d_rows : 1;
         ldak_f = (ak_rows > 0) ? ak_rows : 1;
         ldbk_f = (bk_rows > 0) ? bk_rows : 1;
         ldck_f = (ck_rows > 0) ? ck_rows : 1;
         lddk_f = (dk_rows > 0) ? dk_rows : 1;
         ldx_f = (x_rows > 0) ? x_rows : 1;
         ldz_f = (z_rows > 0) ? z_rows : 1;

         /* Set pointers */
         a_ptr = a_cm; b_ptr = b_cm; c_ptr = c_cm; d_ptr = d_cm;
         ak_ptr = ak_cm; bk_ptr = bk_cm; ck_ptr = ck_cm; dk_ptr = dk_cm;
         x_ptr = x_cm; z_ptr = z_cm;

     } else {
         /* Column-major case - use original arrays */
         lda_f = lda; ldb_f = ldb; ldc_f = ldc; ldd_f = ldd;
         ldak_f = ldak; ldbk_f = ldbk; ldck_f = ldck; lddk_f = lddk;
         ldx_f = ldx; ldz_f = ldz;
         a_ptr = (double*)a; b_ptr = (double*)b; c_ptr = (double*)c; d_ptr = (double*)d; // Cast away const for Fortran call
         ak_ptr = ak; bk_ptr = bk; ck_ptr = ck; dk_ptr = dk;
         x_ptr = x; z_ptr = z;
     }


     /* --- Workspace Allocation --- */

     // Determine fixed workspace sizes
     liwork = MAX(1, MAX(2 * MAX(m2, n), MAX(m, MAX(m2 + np2, n * n))));
     lbwork = MAX(1, 2 * n);

     // Allocate IWORK and BWORK
     iwork = (int*)malloc((size_t)liwork * sizeof(int));
     CHECK_ALLOC(iwork);
     bwork = (int*)malloc((size_t)lbwork * sizeof(int)); // Use int for LOGICAL
     CHECK_ALLOC(bwork);

     // Perform workspace query for DWORK
     ldwork = -1; // Query mode
     F77_FUNC(sb10dd, SB10DD)(&n, &m, &np, &ncon, &nmeas, &gamma,
                              a_ptr, &lda_f, b_ptr, &ldb_f, c_ptr, &ldc_f, d_ptr, &ldd_f,
                              ak_ptr, &ldak_f, bk_ptr, &ldbk_f, ck_ptr, &ldck_f, dk_ptr, &lddk_f,
                              x_ptr, &ldx_f, z_ptr, &ldz_f, rcond, &tol,
                              iwork, &dwork_query, &ldwork, bwork, &info);


     if (info < 0 && info != -30) { goto cleanup; } // Query failed due to invalid argument (allow INFO=-30 from query)
     info = 0; // Reset info after query

     // Get the required dwork size from query result
     ldwork = (int)dwork_query;
     // Check against minimum documented size (complex formula, rely on query + basic checks)
     int min_ldwork = 1;
     // Add basic checks based on components if needed, e.g., MAX(..., 16*n, ...)
     ldwork = MAX(ldwork, min_ldwork);

     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure


     /* --- Call the computational routine --- */
     F77_FUNC(sb10dd, SB10DD)(&n, &m, &np, &ncon, &nmeas, &gamma,
                              a_ptr, &lda_f, b_ptr, &ldb_f, c_ptr, &ldc_f, d_ptr, &ldd_f,
                              ak_ptr, &ldak_f, bk_ptr, &ldbk_f, ck_ptr, &ldck_f, dk_ptr, &lddk_f,
                              x_ptr, &ldx_f, z_ptr, &ldz_f, rcond, &tol,
                              iwork, dwork, &ldwork, bwork, &info);

     /* --- Copy results back to row-major format if needed --- */
     if (row_major && info == 0) {
         // Recalculate sizes needed for copy-back
         size_t ak_rows = n; size_t ak_cols = n; size_t ak_size = ak_rows * ak_cols;
         size_t bk_rows = n; size_t bk_cols = nmeas; size_t bk_size = bk_rows * bk_cols;
         size_t ck_rows = ncon; size_t ck_cols = n; size_t ck_size = ck_rows * ck_cols;
         size_t dk_rows = ncon; size_t dk_cols = nmeas; size_t dk_size = dk_rows * dk_cols;
         size_t x_rows = n; size_t x_cols = n; size_t x_size = x_rows * x_cols;
         size_t z_rows = n; size_t z_cols = n; size_t z_size = z_rows * z_cols;


         if (ak_size > 0) slicot_transpose_to_c(ak_cm, ak, ak_rows, ak_cols, elem_size);
         if (bk_size > 0) slicot_transpose_to_c(bk_cm, bk, bk_rows, bk_cols, elem_size);
         if (ck_size > 0) slicot_transpose_to_c(ck_cm, ck, ck_rows, ck_cols, elem_size);
         if (dk_size > 0) slicot_transpose_to_c(dk_cm, dk, dk_rows, dk_cols, elem_size);
         if (x_size > 0) slicot_transpose_symmetric_to_c(x_cm, x, x_rows, 'U', elem_size); // X is symmetric
         if (z_size > 0) slicot_transpose_symmetric_to_c(z_cm, z, z_rows, 'U', elem_size); // Z is symmetric
         // RCOND modified directly
     }
     // In column-major case, output arrays modified in place.

 cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(iwork);
     free(bwork);
     free(a_cm); free(b_cm); free(c_cm); free(d_cm);
     free(ak_cm); free(bk_cm); free(ck_cm); free(dk_cm);
     free(x_cm); free(z_cm);

     return info;
 }