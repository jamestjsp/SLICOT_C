/**
 * @file sb10ad.c
 * @brief C wrapper implementation for SLICOT routine SB10AD
 *
 * This file provides a C wrapper implementation for the SLICOT routine SB10AD,
 * which computes an H-infinity optimal controller for a continuous-time
 * system using modified Glover's and Doyle's formulas.
 */

 #include <stdlib.h>
 #include <stddef.h> // For size_t

 // Include the header file for this wrapper
 #include "sb10ad.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * Input matrices are const. GAMMA is in/out.
  * Output matrices are AK, BK, CK, DK, AC, BC, CC, DC.
  * Output scalar array is RCOND.
  */
 extern void F77_FUNC(sb10ad, SB10AD)(
     const int* job,         // INTEGER JOB
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     const int* np,          // INTEGER NP
     const int* ncon,        // INTEGER NCON
     const int* nmeas,       // INTEGER NMEAS
     double* gamma,          // DOUBLE PRECISION GAMMA (in/out)
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
     double* ac,             // DOUBLE PRECISION AC(LDAC,*) (output)
     const int* ldac,        // INTEGER LDAC
     double* bc,             // DOUBLE PRECISION BC(LDBC,*) (output)
     const int* ldbc,        // INTEGER LDBC
     double* cc,             // DOUBLE PRECISION CC(LDCC,*) (output)
     const int* ldcc,        // INTEGER LDCC
     double* dc,             // DOUBLE PRECISION DC(LDDC,*) (output)
     const int* lddc,        // INTEGER LDDC
     double* rcond,          // DOUBLE PRECISION RCOND(4) (output)
     const double* gtol,     // DOUBLE PRECISION GTOL
     const double* actol,    // DOUBLE PRECISION ACTOL
     int* iwork,             // INTEGER IWORK(*)
     const int* liwork,      // INTEGER LIWORK
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* bwork,             // LOGICAL BWORK(*) -> int*
     const int* lbwork,      // INTEGER LBWORK
     int* info               // INTEGER INFO (output)
 );


 /* C wrapper function definition */
 SLICOT_C_WRAPPER_API
 int slicot_sb10ad(int job, int n, int m, int np, int ncon, int nmeas,
                   double* gamma, const double* a, int lda,
                   const double* b, int ldb, const double* c, int ldc,
                   const double* d, int ldd, double* ak, int ldak,
                   double* bk, int ldbk, double* ck, int ldck,
                   double* dk, int lddk, double* ac, int ldac,
                   double* bc, int ldbc, double* cc, int ldcc,
                   double* dc, int lddc, double* rcond, double gtol,
                   double actol, int row_major)
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
     double *ac_cm = NULL, *bc_cm = NULL, *cc_cm = NULL, *dc_cm = NULL;
     const double *a_ptr, *b_ptr, *c_ptr, *d_ptr; // Input pointers
     double *ak_ptr, *bk_ptr, *ck_ptr, *dk_ptr; // Output pointers (K)
     double *ac_ptr, *bc_ptr, *cc_ptr, *dc_ptr; // Output pointers (CL)
     int lda_f, ldb_f, ldc_f, ldd_f, ldak_f, ldbk_f, ldck_f, lddk_f;
     int ldac_f, ldbc_f, ldcc_f, lddc_f;

     /* --- Input Parameter Validation --- */

     if (job < 1 || job > 4) { info = -1; goto cleanup; }
     if (n < 0) { info = -2; goto cleanup; }
     if (m < 0) { info = -3; goto cleanup; }
     if (np < 0) { info = -4; goto cleanup; }
     if (ncon < 0 || ncon > m) { info = -5; goto cleanup; }
     if (nmeas < 0 || nmeas > np) { info = -6; goto cleanup; }
     if (gamma == NULL || *gamma < 0.0) { info = -7; goto cleanup; }
     // GTOL, ACTOL checked by Fortran


     // Calculate partition dimensions needed for validation/workspace
     int m1 = m - ncon;      // Columns of B1, D11, D21
     int np2 = nmeas;        // Rows of C2, D21, D22
     int np1 = np - np2;     // Rows of C1, D11, D12
     int m2 = ncon;          // Columns of B2, D12, D22

     // Further dimension checks related to partitioning
     if (np1 < 0) { info = -6; goto cleanup; } // Implied by nmeas <= np
     if (m1 < 0)  { info = -5; goto cleanup; } // Implied by ncon <= m
     // These checks might be redundant based on earlier checks, but good for clarity
     if (np - nmeas < 0) { info = -6; goto cleanup; } // Check C1 rows >= 0
     if (m - ncon < 0)  { info = -5; goto cleanup; } // Check B1 cols >= 0


     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, n);
     int min_ldb_f = MAX(1, n);
     int min_ldc_f = MAX(1, np);
     int min_ldd_f = MAX(1, np);
     int min_ldak_f = MAX(1, n);
     int min_ldbk_f = MAX(1, n);
     int min_ldck_f = MAX(1, ncon);
     int min_lddk_f = MAX(1, ncon);
     int min_ldac_f = MAX(1, 2 * n);
     int min_ldbc_f = MAX(1, 2 * n);
     int min_ldcc_f = MAX(1, np1);
     int min_lddc_f = MAX(1, np1);

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
         int min_ldac_rm_cols = 2 * n;
         int min_ldbc_rm_cols = m1;
         int min_ldcc_rm_cols = 2 * n;
         int min_lddc_rm_cols = m1;

         if (lda < min_lda_rm_cols) { info = -9; goto cleanup; }
         if (ldb < min_ldb_rm_cols) { info = -11; goto cleanup; }
         if (ldc < min_ldc_rm_cols) { info = -13; goto cleanup; }
         if (ldd < min_ldd_rm_cols) { info = -15; goto cleanup; }
         if (ldak < min_ldak_rm_cols) { info = -17; goto cleanup; }
         if (ldbk < min_ldbk_rm_cols) { info = -19; goto cleanup; }
         if (ldck < min_ldck_rm_cols) { info = -21; goto cleanup; }
         if (lddk < min_lddk_rm_cols) { info = -23; goto cleanup; }
         if (ldac < min_ldac_rm_cols) { info = -25; goto cleanup; }
         if (ldbc < min_ldbc_rm_cols) { info = -27; goto cleanup; }
         if (ldcc < min_ldcc_rm_cols) { info = -29; goto cleanup; }
         if (lddc < min_lddc_rm_cols) { info = -31; goto cleanup; }
     } else {
         // For column-major C, LD is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -9; goto cleanup; }
         if (ldb < min_ldb_f) { info = -11; goto cleanup; }
         if (ldc < min_ldc_f) { info = -13; goto cleanup; }
         if (ldd < min_ldd_f) { info = -15; goto cleanup; }
         if (ldak < min_ldak_f) { info = -17; goto cleanup; }
         if (ldbk < min_ldbk_f) { info = -19; goto cleanup; }
         if (ldck < min_ldck_f) { info = -21; goto cleanup; }
         if (lddk < min_lddk_f) { info = -23; goto cleanup; }
         if (ldac < min_ldac_f) { info = -25; goto cleanup; }
         if (ldbc < min_ldbc_f) { info = -27; goto cleanup; }
         if (ldcc < min_ldcc_f) { info = -29; goto cleanup; }
         if (lddc < min_lddc_f) { info = -31; goto cleanup; }
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
         size_t ac_rows = 2*n; size_t ac_cols = 2*n; size_t ac_size = ac_rows * ac_cols;
         size_t bc_rows = 2*n; size_t bc_cols = m1; size_t bc_size = bc_rows * bc_cols;
         size_t cc_rows = np1; size_t cc_cols = 2*n; size_t cc_size = cc_rows * cc_cols;
         size_t dc_rows = np1; size_t dc_cols = m1; size_t dc_size = dc_rows * dc_cols;

         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
         if (d_size > 0) { d_cm = (double*)malloc(d_size * elem_size); CHECK_ALLOC(d_cm); }
         if (ak_size > 0) { ak_cm = (double*)malloc(ak_size * elem_size); CHECK_ALLOC(ak_cm); } // Output
         if (bk_size > 0) { bk_cm = (double*)malloc(bk_size * elem_size); CHECK_ALLOC(bk_cm); } // Output
         if (ck_size > 0) { ck_cm = (double*)malloc(ck_size * elem_size); CHECK_ALLOC(ck_cm); } // Output
         if (dk_size > 0) { dk_cm = (double*)malloc(dk_size * elem_size); CHECK_ALLOC(dk_cm); } // Output
         if (ac_size > 0) { ac_cm = (double*)malloc(ac_size * elem_size); CHECK_ALLOC(ac_cm); } // Output
         if (bc_size > 0) { bc_cm = (double*)malloc(bc_size * elem_size); CHECK_ALLOC(bc_cm); } // Output
         if (cc_size > 0) { cc_cm = (double*)malloc(cc_size * elem_size); CHECK_ALLOC(cc_cm); } // Output
         if (dc_size > 0) { dc_cm = (double*)malloc(dc_size * elem_size); CHECK_ALLOC(dc_cm); } // Output

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
         ldac_f = (ac_rows > 0) ? ac_rows : 1;
         ldbc_f = (bc_rows > 0) ? bc_rows : 1;
         ldcc_f = (cc_rows > 0) ? cc_rows : 1;
         lddc_f = (dc_rows > 0) ? dc_rows : 1;

         /* Set pointers */
         a_ptr = a_cm; b_ptr = b_cm; c_ptr = c_cm; d_ptr = d_cm;
         ak_ptr = ak_cm; bk_ptr = bk_cm; ck_ptr = ck_cm; dk_ptr = dk_cm;
         ac_ptr = ac_cm; bc_ptr = bc_cm; cc_ptr = cc_cm; dc_ptr = dc_cm;

     } else {
         /* Column-major case - use original arrays */
         lda_f = lda; ldb_f = ldb; ldc_f = ldc; ldd_f = ldd;
         ldak_f = ldak; ldbk_f = ldbk; ldck_f = ldck; lddk_f = lddk;
         ldac_f = ldac; ldbc_f = ldbc; ldcc_f = ldcc; lddc_f = lddc;
         a_ptr = (double*)a; b_ptr = (double*)b; c_ptr = (double*)c; d_ptr = (double*)d; // Cast away const for Fortran call
         ak_ptr = ak; bk_ptr = bk; ck_ptr = ck; dk_ptr = dk;
         ac_ptr = ac; bc_ptr = bc; cc_ptr = cc; dc_ptr = dc;
     }

     /* --- Workspace Allocation --- */

     // Determine fixed workspace sizes
     liwork = MAX(1, MAX(2 * MAX(n, MAX(m1, MAX(np1, MAX(m2, np2)))), n * n));
     lbwork = MAX(1, 2 * n);

     // Allocate IWORK and BWORK
     iwork = (int*)malloc((size_t)liwork * sizeof(int));
     CHECK_ALLOC(iwork);
     bwork = (int*)malloc((size_t)lbwork * sizeof(int)); // Use int for LOGICAL
     CHECK_ALLOC(bwork);

     // Perform workspace query for DWORK
     ldwork = -1; // Query mode
     F77_FUNC(sb10ad, SB10AD)(&job, &n, &m, &np, &ncon, &nmeas, gamma,
                              a_ptr, &lda_f, b_ptr, &ldb_f, c_ptr, &ldc_f, d_ptr, &ldd_f,
                              ak_ptr, &ldak_f, bk_ptr, &ldbk_f, ck_ptr, &ldck_f, dk_ptr, &lddk_f,
                              ac_ptr, &ldac_f, bc_ptr, &ldbc_f, cc_ptr, &ldcc_f, dc_ptr, &lddc_f,
                              rcond, &gtol, &actol,
                              iwork, &liwork, &dwork_query, &ldwork,
                              bwork, &lbwork, &info);

     if (info < 0 && info != -37) { goto cleanup; } // Query failed due to invalid argument (allow INFO=-37 from query)
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
     F77_FUNC(sb10ad, SB10AD)(&job, &n, &m, &np, &ncon, &nmeas, gamma,
                              a_ptr, &lda_f, b_ptr, &ldb_f, c_ptr, &ldc_f, d_ptr, &ldd_f,
                              ak_ptr, &ldak_f, bk_ptr, &ldbk_f, ck_ptr, &ldck_f, dk_ptr, &lddk_f,
                              ac_ptr, &ldac_f, bc_ptr, &ldbc_f, cc_ptr, &ldcc_f, dc_ptr, &lddc_f,
                              rcond, &gtol, &actol,
                              iwork, &liwork, dwork, &ldwork,
                              bwork, &lbwork, &info);

     /* --- Copy results back to row-major format if needed --- */
     if (row_major && info == 0) {
         // Recalculate sizes needed for copy-back
         size_t ak_rows = n; size_t ak_cols = n; size_t ak_size = ak_rows * ak_cols;
         size_t bk_rows = n; size_t bk_cols = nmeas; size_t bk_size = bk_rows * bk_cols;
         size_t ck_rows = ncon; size_t ck_cols = n; size_t ck_size = ck_rows * ck_cols;
         size_t dk_rows = ncon; size_t dk_cols = nmeas; size_t dk_size = dk_rows * dk_cols;
         size_t ac_rows = 2*n; size_t ac_cols = 2*n; size_t ac_size = ac_rows * ac_cols;
         size_t bc_rows = 2*n; size_t bc_cols = m1; size_t bc_size = bc_rows * bc_cols;
         size_t cc_rows = np1; size_t cc_cols = 2*n; size_t cc_size = cc_rows * cc_cols;
         size_t dc_rows = np1; size_t dc_cols = m1; size_t dc_size = dc_rows * dc_cols;

         if (ak_size > 0) slicot_transpose_to_c(ak_cm, ak, ak_rows, ak_cols, elem_size);
         if (bk_size > 0) slicot_transpose_to_c(bk_cm, bk, bk_rows, bk_cols, elem_size);
         if (ck_size > 0) slicot_transpose_to_c(ck_cm, ck, ck_rows, ck_cols, elem_size);
         if (dk_size > 0) slicot_transpose_to_c(dk_cm, dk, dk_rows, dk_cols, elem_size);
         if (ac_size > 0) slicot_transpose_to_c(ac_cm, ac, ac_rows, ac_cols, elem_size);
         if (bc_size > 0) slicot_transpose_to_c(bc_cm, bc, bc_rows, bc_cols, elem_size);
         if (cc_size > 0) slicot_transpose_to_c(cc_cm, cc, cc_rows, cc_cols, elem_size);
         if (dc_size > 0) slicot_transpose_to_c(dc_cm, dc, dc_rows, dc_cols, elem_size);
         // GAMMA, RCOND modified directly
     }
     // In column-major case, output arrays modified in place.

 cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(iwork);
     free(bwork);
     free(a_cm); free(b_cm); free(c_cm); free(d_cm);
     free(ak_cm); free(bk_cm); free(ck_cm); free(dk_cm);
     free(ac_cm); free(bc_cm); free(cc_cm); free(dc_cm);

     return info;
 }