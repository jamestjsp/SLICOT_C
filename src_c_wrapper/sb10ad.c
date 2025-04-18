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
     const double* a,        // DOUBLE PRECISION A(LDA,*)
     const int* lda,         // INTEGER LDA
     const double* b,        // DOUBLE PRECISION B(LDB,*)
     const int* ldb,         // INTEGER LDB
     const double* c,        // DOUBLE PRECISION C(LDC,*)
     const int* ldc,         // INTEGER LDC
     const double* d,        // DOUBLE PRECISION D(LDD,*)
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
     int liwork = 0;
     int lbwork = 0;
     double dwork_query;
     double* dwork = NULL;
     int* iwork = NULL;
     int* bwork = NULL; // Map LOGICAL to int
 
     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL;
     double *ak_cm = NULL, *bk_cm = NULL, *ck_cm = NULL, *dk_cm = NULL;
     double *ac_cm = NULL, *bc_cm = NULL, *cc_cm = NULL, *dc_cm = NULL;
 
     /* --- Input Parameter Validation --- */
 
     if (n < 0) { info = -2; goto cleanup; }
     if (m < 0) { info = -3; goto cleanup; }
     if (np < 0) { info = -4; goto cleanup; }
     if (ncon < 0 || ncon > m) { info = -5; goto cleanup; }
     if (nmeas < 0 || nmeas > np) { info = -6; goto cleanup; }
     if (gamma == NULL || *gamma < 0.0) { info = -7; goto cleanup; }
     if (job < 1 || job > 4) { info = -1; goto cleanup; }
     // Further dimension checks related to partitioning
     if (np - nmeas < ncon) { info = -5; goto cleanup; } // Check C1 rows >= NCON
     if (m - ncon < nmeas) { info = -6; goto cleanup; } // Check B1 cols >= NMEAS
 
     // Check leading dimensions based on storage order
     int m1 = m - ncon;      // Columns of B1, D11, D21
     int np2 = nmeas;        // Rows of C2, D21, D22
     int np1 = np - np2;     // Rows of C1, D11, D12
     int m2 = ncon;          // Columns of B2, D12, D22
 
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
 
     /* --- Workspace Allocation --- */
 
     // Determine fixed workspace sizes
     liwork = MAX(1, MAX(2 * MAX(n, MAX(m1, MAX(np1, MAX(m2, np2)))), n * n));
     lbwork = MAX(1, 2 * n);
 
     // Allocate IWORK and BWORK
     iwork = (int*)malloc((size_t)liwork * sizeof(int));
     CHECK_ALLOC(iwork);
     bwork = (int*)malloc((size_t)lbwork * sizeof(int)); // Use int for LOGICAL
     CHECK_ALLOC(bwork);
 
     // Allocate DWORK based on query
     ldwork = -1; // Query mode
     F77_FUNC(sb10ad, SB10AD)(&job, &n, &m, &np, &ncon, &nmeas, gamma,
                              a, &lda, b, &ldb, c, &ldc, d, &ldd,
                              ak, &ldak, bk, &ldbk, ck, &ldck, dk, &lddk,
                              ac, &ldac, bc, &ldbc, cc, &ldcc, dc, &lddc,
                              rcond, &gtol, &actol,
                              iwork, &liwork, &dwork_query, &ldwork,
                              bwork, &lbwork, &info);
 
     if (info < 0) { goto cleanup; } // Query failed due to invalid argument
     info = 0; // Reset info after query
 
     // Get the required dwork size from query result
     ldwork = (int)dwork_query;
     // Check against minimum documented size (complex formula, rely on query + basic checks)
     int min_ldwork = 1;
     // Add basic checks based on components if needed, e.g., MAX(..., 16*n, ...)
     ldwork = MAX(ldwork, min_ldwork);
 
     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure
 
     /* --- Prepare Arrays and Call Fortran Routine --- */
     size_t elem_size = sizeof(double);
 
     // Determine sizes for potential copies
     size_t a_size = (size_t)n * n;
     size_t b_size = (size_t)n * m;
     size_t c_size = (size_t)np * n;
     size_t d_size = (size_t)np * m;
     size_t ak_size = (size_t)n * n;
     size_t bk_size = (size_t)n * nmeas;
     size_t ck_size = (size_t)ncon * n;
     size_t dk_size = (size_t)ncon * nmeas;
     size_t ac_size = (size_t)(2*n) * (2*n);
     size_t bc_size = (size_t)(2*n) * m1;
     size_t cc_size = (size_t)np1 * (2*n);
     size_t dc_size = (size_t)np1 * m1;
 
 
     if (row_major) {
         /* --- Row-Major Case --- */
 
         /* Allocate memory for column-major copies */
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
         if (d_size > 0) { d_cm = (double*)malloc(d_size * elem_size); CHECK_ALLOC(d_cm); }
         if (ak_size > 0) { ak_cm = (double*)malloc(ak_size * elem_size); CHECK_ALLOC(ak_cm); }
         if (bk_size > 0) { bk_cm = (double*)malloc(bk_size * elem_size); CHECK_ALLOC(bk_cm); }
         if (ck_size > 0) { ck_cm = (double*)malloc(ck_size * elem_size); CHECK_ALLOC(ck_cm); }
         if (dk_size > 0) { dk_cm = (double*)malloc(dk_size * elem_size); CHECK_ALLOC(dk_cm); }
         if (ac_size > 0) { ac_cm = (double*)malloc(ac_size * elem_size); CHECK_ALLOC(ac_cm); }
         if (bc_size > 0) { bc_cm = (double*)malloc(bc_size * elem_size); CHECK_ALLOC(bc_cm); }
         if (cc_size > 0) { cc_cm = (double*)malloc(cc_size * elem_size); CHECK_ALLOC(cc_cm); }
         if (dc_size > 0) { dc_cm = (double*)malloc(dc_size * elem_size); CHECK_ALLOC(dc_cm); }
 
 
         /* Transpose C inputs to Fortran copies */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, n, n, elem_size);
         if (b_size > 0) slicot_transpose_to_fortran(b, b_cm, n, m, elem_size);
         if (c_size > 0) slicot_transpose_to_fortran(c, c_cm, np, n, elem_size);
         if (d_size > 0) slicot_transpose_to_fortran(d, d_cm, np, m, elem_size);
 
         /* Fortran leading dimensions (use input LDs for call) */
         // Note: Fortran LDs are passed directly as they define the stride in the original data
 
         /* Call the Fortran routine */
         F77_FUNC(sb10ad, SB10AD)(&job, &n, &m, &np, &ncon, &nmeas, gamma,
                                  a_cm, &lda, b_cm, &ldb, c_cm, &ldc, d_cm, &ldd,
                                  ak_cm, &ldak, bk_cm, &ldbk, ck_cm, &ldck, dk_cm, &lddk,
                                  ac_cm, &ldac, bc_cm, &ldbc, cc_cm, &ldcc, dc_cm, &lddc,
                                  rcond, &gtol, &actol,
                                  iwork, &liwork, dwork, &ldwork,
                                  bwork, &lbwork, &info);
 
         /* Copy back results from column-major temps to original row-major arrays */
         if (info == 0) {
             if (ak_size > 0) slicot_transpose_to_c(ak_cm, ak, n, n, elem_size);
             if (bk_size > 0) slicot_transpose_to_c(bk_cm, bk, n, nmeas, elem_size);
             if (ck_size > 0) slicot_transpose_to_c(ck_cm, ck, ncon, n, elem_size);
             if (dk_size > 0) slicot_transpose_to_c(dk_cm, dk, ncon, nmeas, elem_size);
             if (ac_size > 0) slicot_transpose_to_c(ac_cm, ac, 2*n, 2*n, elem_size);
             if (bc_size > 0) slicot_transpose_to_c(bc_cm, bc, 2*n, m1, elem_size);
             if (cc_size > 0) slicot_transpose_to_c(cc_cm, cc, np1, 2*n, elem_size);
             if (dc_size > 0) slicot_transpose_to_c(dc_cm, dc, np1, m1, elem_size);
             // GAMMA, RCOND modified directly
         }
         /* Temps freed in cleanup */
 
     } else {
         /* --- Column-Major Case --- */
 
         /* Call the Fortran routine directly with user-provided arrays */
         F77_FUNC(sb10ad, SB10AD)(&job, &n, &m, &np, &ncon, &nmeas, gamma,
                                  a, &lda, b, &ldb, c, &ldc, d, &ldd,
                                  ak, &ldak, bk, &ldbk, ck, &ldck, dk, &lddk,
                                  ac, &ldac, bc, &ldbc, cc, &ldcc, dc, &lddc,
                                  rcond, &gtol, &actol,
                                  iwork, &liwork, dwork, &ldwork,
                                  bwork, &lbwork, &info);
         // Output arrays modified in place
     }
 
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
 