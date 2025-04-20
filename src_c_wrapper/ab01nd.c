/**
 * @file ab01nd.c
 * @brief C wrapper implementation for SLICOT routine AB01ND
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB01ND,
 * which finds a controllable realization for a linear time-invariant
 * multi-input system using orthogonal transformations.
 */

 #include <stdlib.h>
 #include <string.h> // For memset
 #include <ctype.h>  // For toupper
 
 // Include the header file for this wrapper
 #include "ab01nd.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines, set_identity
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions
 
 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * This handles potential name mangling issues between C and Fortran compilers.
  * All arguments are passed by reference (pointers).
  * The hidden string length argument for 'jobz' is explicitly included at the end.
  */
 extern void F77_FUNC(ab01nd, AB01ND)(
     const char* jobz,       // CHARACTER*1 JOBZ
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     double* a,              // DOUBLE PRECISION A(LDA,*) (in/out)
     const int* lda,         // INTEGER LDA
     double* b,              // DOUBLE PRECISION B(LDB,*) (in/out)
     const int* ldb,         // INTEGER LDB
     int* ncont,             // INTEGER NCONT (output)
     int* indcon,            // INTEGER INDCON (output)
     int* nblk,              // INTEGER NBLK(*) (output)
     double* z,              // DOUBLE PRECISION Z(LDZ,*) (output)
     const int* ldz,         // INTEGER LDZ
     double* tau,            // DOUBLE PRECISION TAU(*) (output)
     const double* tol,      // DOUBLE PRECISION TOL
     int* iwork,             // INTEGER IWORK(*)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* info,              // INTEGER INFO (output)
     int jobz_len            // Hidden length argument for jobz (integer)
 );
 
 
 /* C wrapper function definition */
 SLICOT_C_WRAPPER_API
 int slicot_ab01nd(char jobz, int n, int m,
                   double* a, int lda,
                   double* b, int ldb,
                   int* ncont, int* indcon, int* nblk,
                   double* z, int ldz,
                   double* tau, double tol,
                   int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork;
     double* dwork = NULL;
     int* iwork = NULL;
     int jobz_len = 1; // Fortran expects 1-based length for strings

     char jobz_upper = toupper(jobz);

     /* Pointers for column-major copies if needed */
     double *a_cm = NULL;
     double *b_cm = NULL;
     double *z_cm = NULL; // Needed if row_major and jobz != 'N'

     /* --- Input Parameter Validation --- */

     if (n < 0) { info = -2; goto cleanup; }
     if (m < 0) { info = -3; goto cleanup; }

     if (jobz_upper != 'N' && jobz_upper != 'F' && jobz_upper != 'I') {
         info = -1; goto cleanup;
     }

     // Check leading dimensions based on storage order and JOBZ
     int min_lda_f = MAX(1, n);
     int min_ldb_f = MAX(1, n);
     int min_ldz_f = (jobz_upper == 'N') ? 1 : MAX(1, n);

     if (row_major) {
         // For row-major C, LDA/LDB/LDZ is the number of columns
         int min_lda_rm_cols = n;
         int min_ldb_rm_cols = m;
         int min_ldz_rm_cols = (jobz_upper == 'N') ? 1 : n; // Need n columns for Z
         if (lda < min_lda_rm_cols) { info = -5; goto cleanup; }
         if (ldb < min_ldb_rm_cols) { info = -7; goto cleanup; }
         if (jobz_upper != 'N' && ldz < min_ldz_rm_cols) { info = -12; goto cleanup; }
     } else {
         // For column-major C, LDA/LDB/LDZ is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -5; goto cleanup; }
         if (ldb < min_ldb_f) { info = -7; goto cleanup; }
         if (ldz < min_ldz_f) { info = -12; goto cleanup; } // Check even if JOBZ='N' as per docs
     }

     /* --- Prepare arrays for column-major format if using row-major --- */
     int lda_f, ldb_f, ldz_f;
     double *a_ptr, *b_ptr, *z_ptr;
     
     if (row_major) {
         /* Allocate and convert A from row-major to column-major */
         if (n > 0) {
             a_cm = (double*)malloc((size_t)n * n * sizeof(double));
             CHECK_ALLOC(a_cm);
             slicot_transpose_to_fortran(a, a_cm, n, n, sizeof(double));
         }
         
         /* Allocate and convert B from row-major to column-major */
         if (n > 0 && m > 0) {
             b_cm = (double*)malloc((size_t)n * m * sizeof(double));
             CHECK_ALLOC(b_cm);
             slicot_transpose_to_fortran(b, b_cm, n, m, sizeof(double));
         }
         
         /* Allocate Z if needed and initialize to identity if JOBZ='I' */
         if (jobz_upper != 'N' && n > 0) {
             z_cm = (double*)malloc((size_t)n * n * sizeof(double));
             CHECK_ALLOC(z_cm);
             
             if (jobz_upper == 'I') {
                 set_identity(n, z_cm, n, 0); // 0 for column-major
             }
         }
         
         lda_f = n;            // Column-major LDA for Fortran
         ldb_f = n;            // Column-major LDB for Fortran
         ldz_f = (jobz_upper == 'N') ? 1 : n; // Column-major LDZ for Fortran
         a_ptr = a_cm;
         b_ptr = b_cm;
         z_ptr = (jobz_upper == 'N') ? NULL : z_cm;
     } else {
         /* Column-major case - use original arrays */
         lda_f = lda;
         ldb_f = ldb;
         ldz_f = ldz;
         a_ptr = a;
         b_ptr = b;
         z_ptr = (jobz_upper == 'N') ? NULL : z;
         
         /* Initialize Z to identity if JOBZ='I' in column-major format */
         if (jobz_upper == 'I' && n > 0) {
             set_identity(n, z, ldz, 0); // 0 for column-major
         }
     }

     /* --- Workspace allocation --- */
     /* 
      * Skip the workspace query and directly allocate the minimum required workspace.
      * According to documentation, LDWORK >= MAX(1, N, 3*M) is required.
      */
     ldwork = MAX(1, MAX(n, 3 * m));
     
     /* Allocate a bit more than the minimum for better performance */
     ldwork = MAX(ldwork, MAX(n, 3 * m) * 2);
     
     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork);
     
     /* Allocate integer workspace */
     iwork = (int*)malloc((size_t)m * sizeof(int));
     CHECK_ALLOC(iwork);
     
     /* --- Call the computational routine --- */
     F77_FUNC(ab01nd, AB01ND)(&jobz_upper, &n, &m, 
                            a_ptr, &lda_f, 
                            b_ptr, &ldb_f, 
                            ncont, indcon, nblk, 
                            z_ptr, &ldz_f, 
                            tau, &tol, 
                            iwork, dwork, &ldwork, &info, 
                            jobz_len);
     
     /* --- Copy results back to row-major format if needed --- */
     if (row_major && info == 0) {
         if (n > 0) {
             slicot_transpose_to_c(a_cm, a, n, n, sizeof(double));
         }
         
         if (n > 0 && m > 0) {
             slicot_transpose_to_c(b_cm, b, n, m, sizeof(double));
         }
         
         if (jobz_upper != 'N' && n > 0) {
             slicot_transpose_to_c(z_cm, z, n, n, sizeof(double));
         }
     }
     
cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(iwork);
     free(a_cm);
     free(b_cm);
     free(z_cm);
     
     return info;
 }
