/**
 * @file tg01fd.c
 * @brief C wrapper implementation for SLICOT routine TG01FD
 *
 * This file provides a C wrapper implementation for the SLICOT routine TG01FD,
 * which computes an orthogonal reduction of a descriptor system
 * (A-lambda E,B,C) to a SVD-like coordinate form.
 */

 #include <stdlib.h>
 #include <ctype.h>
 #include <string.h> // For memcpy
 #include <stddef.h> // For size_t

 // Include the header file for this wrapper
 // #include "tg01fd.h" // Assuming a header file exists
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, MIN, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines, set_identity
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * Hidden lengths for CHARACTER arguments are added at the end.
  */
 extern void F77_FUNC(tg01fd, TG01FD)(
     const char* compq,      // CHARACTER*1 COMPQ
     const char* compz,      // CHARACTER*1 COMPZ
     const char* joba,       // CHARACTER*1 JOBA
     const int* l,           // INTEGER L
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     const int* p,           // INTEGER P
     double* a,              // DOUBLE PRECISION A(LDA,*) (in/out)
     const int* lda,         // INTEGER LDA
     double* e,              // DOUBLE PRECISION E(LDE,*) (in/out)
     const int* lde,         // INTEGER LDE
     double* b,              // DOUBLE PRECISION B(LDB,*) (in/out)
     const int* ldb,         // INTEGER LDB
     double* c,              // DOUBLE PRECISION C(LDC,*) (in/out)
     const int* ldc,         // INTEGER LDC
     double* q,              // DOUBLE PRECISION Q(LDQ,*) (in/out)
     const int* ldq,         // INTEGER LDQ
     double* z,              // DOUBLE PRECISION Z(LDZ,*) (in/out)
     const int* ldz,         // INTEGER LDZ
     int* ranke,             // INTEGER RANKE (output)
     int* rnka22,            // INTEGER RNKA22 (output)
     const double* tol,      // DOUBLE PRECISION TOL
     int* iwork,             // INTEGER IWORK(*)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* info,              // INTEGER INFO (output)
     int compq_len,          // Hidden length
     int compz_len,          // Hidden length
     int joba_len            // Hidden length
 );


 /* C wrapper function definition */
 SLICOT_C_WRAPPER_API
 int slicot_tg01fd(char compq, char compz, char joba, int l, int n, int m, int p,
                   double* a, int lda, double* e, int lde,
                   double* b, int ldb, double* c, int ldc,
                   double* q, int ldq, double* z, int ldz,
                   int* ranke, int* rnka22, double tol, int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = -1; /* Use -1 for workspace query */
     double dwork_query;
     double* dwork = NULL;
     int* iwork = NULL;
     int iwork_size = 0;

     const int compq_len = 1, compz_len = 1, joba_len = 1;
     char compq_upper = toupper(compq);
     char compz_upper = toupper(compz);
     char joba_upper = toupper(joba);

     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *e_cm = NULL, *b_cm = NULL, *c_cm = NULL;
     double *q_cm = NULL, *z_cm = NULL;
     double *a_ptr, *e_ptr, *b_ptr, *c_ptr, *q_ptr, *z_ptr;
     int lda_f, lde_f, ldb_f, ldc_f, ldq_f, ldz_f;

     /* --- Input Parameter Validation --- */
     if (compq_upper != 'N' && compq_upper != 'I' && compq_upper != 'U') { info = -1; goto cleanup; }
     if (compz_upper != 'N' && compz_upper != 'I' && compz_upper != 'U') { info = -2; goto cleanup; }
     if (joba_upper != 'N' && joba_upper != 'R' && joba_upper != 'T') { info = -3; goto cleanup; }
     if (l < 0) { info = -4; goto cleanup; }
     if (n < 0) { info = -5; goto cleanup; }
     if (m < 0) { info = -6; goto cleanup; }
     if (p < 0) { info = -7; goto cleanup; }
     // Tol checked by Fortran, but user check is good practice:
     // if (tol <= 0.0) { info = -21; goto cleanup; } // Allow TOL=0? Check docs. SLICOT usually requires TOL >= 0.

     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, l);
     int min_lde_f = MAX(1, l);
     int min_ldb_f = (m > 0) ? MAX(1, l) : 1; // B is L x M
     int min_ldc_f = MAX(1, p);             // C is P x N
     int min_ldq_f = (compq_upper == 'N') ? 1 : MAX(1, l); // Q is L x L
     int min_ldz_f = (compz_upper == 'N') ? 1 : MAX(1, n); // Z is N x N


     if (row_major) {
         // For row-major C, LD is number of columns
         int min_lda_rm_cols = n;
         int min_lde_rm_cols = n;
         int min_ldb_rm_cols = m;
         int min_ldc_rm_cols = n;
         int min_ldq_rm_cols = (compq_upper == 'N') ? 1 : l;
         int min_ldz_rm_cols = (compz_upper == 'N') ? 1 : n;
         if (lda < min_lda_rm_cols) { info = -9; goto cleanup; }
         if (lde < min_lde_rm_cols) { info = -11; goto cleanup; }
         if (m > 0 && ldb < min_ldb_rm_cols) { info = -13; goto cleanup; } // Check only if m > 0
         if (ldc < min_ldc_rm_cols) { info = -15; goto cleanup; }
         if (ldq < min_ldq_rm_cols) { info = -17; goto cleanup; }
         if (ldz < min_ldz_rm_cols) { info = -19; goto cleanup; }
     } else {
         // For column-major C, LD is number of rows (Fortran style)
         if (lda < min_lda_f) { info = -9; goto cleanup; }
         if (lde < min_lde_f) { info = -11; goto cleanup; }
         if (ldb < min_ldb_f) { info = -13; goto cleanup; } // Check even if m=0
         if (ldc < min_ldc_f) { info = -15; goto cleanup; }
         if (ldq < min_ldq_f) { info = -17; goto cleanup; }
         if (ldz < min_ldz_f) { info = -19; goto cleanup; }
     }

     /* --- Prepare arrays for column-major format if using row-major --- */
     size_t elem_size = sizeof(double);
     if (row_major) {
         // Determine sizes for potential copies
         size_t a_rows = l; size_t a_cols = n; size_t a_size = a_rows * a_cols;
         size_t e_rows = l; size_t e_cols = n; size_t e_size = e_rows * e_cols;
         size_t b_rows = l; size_t b_cols = m; size_t b_size = b_rows * b_cols;
         size_t c_rows = p; size_t c_cols = n; size_t c_size = c_rows * c_cols;
         size_t q_rows = l; size_t q_cols = l; size_t q_size = q_rows * q_cols;
         size_t z_rows = n; size_t z_cols = n; size_t z_size = z_rows * z_cols;

         /* Allocate memory for column-major copies */
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (e_size > 0) { e_cm = (double*)malloc(e_size * elem_size); CHECK_ALLOC(e_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
         if (compq_upper != 'N' && q_size > 0) { q_cm = (double*)malloc(q_size * elem_size); CHECK_ALLOC(q_cm); }
         if (compz_upper != 'N' && z_size > 0) { z_cm = (double*)malloc(z_size * elem_size); CHECK_ALLOC(z_cm); }

         /* Transpose C inputs to Fortran copies */
         if (a_cm) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         if (e_cm) slicot_transpose_to_fortran(e, e_cm, e_rows, e_cols, elem_size);
         if (b_cm) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, elem_size);
         if (c_cm) slicot_transpose_to_fortran(c, c_cm, c_rows, c_cols, elem_size);
         if (compq_upper == 'U' && q_cm) { // Copy Q if updating
             slicot_transpose_to_fortran(q, q_cm, q_rows, q_cols, elem_size);
         } else if (compq_upper == 'I' && q_cm) { // Initialize Q if requested
             set_identity(l, q_cm, l, 0); // Initialize CM Q to identity
         }
         if (compz_upper == 'U' && z_cm) { // Copy Z if updating
             slicot_transpose_to_fortran(z, z_cm, z_rows, z_cols, elem_size);
         } else if (compz_upper == 'I' && z_cm) { // Initialize Z if requested
             set_identity(n, z_cm, n, 0); // Initialize CM Z to identity
         }

         /* Fortran leading dimensions */
         lda_f = (l > 0) ? l : 1;
         lde_f = (l > 0) ? l : 1;
         ldb_f = (l > 0) ? l : 1;
         ldc_f = (p > 0) ? p : 1;
         ldq_f = (l > 0) ? l : 1;
         ldz_f = (n > 0) ? n : 1;

         /* Set pointers */
         a_ptr = a_cm; e_ptr = e_cm; b_ptr = b_cm; c_ptr = c_cm;
         q_ptr = q_cm; z_ptr = z_cm;

     } else {
         /* Column-major case - use original arrays */
         lda_f = lda; lde_f = lde; ldb_f = ldb; ldc_f = ldc; ldq_f = ldq; ldz_f = ldz;
         a_ptr = a; e_ptr = e; b_ptr = b; c_ptr = c; q_ptr = q; z_ptr = z;
         /* Initialize Q/Z if requested */
         if (compq_upper == 'I') set_identity(l, q_ptr, ldq_f, 0);
         if (compz_upper == 'I') set_identity(n, z_ptr, ldz_f, 0);
     }


     /* --- Workspace Allocation --- */

     // Allocate IWORK
     iwork_size = n; // Size is N
     if (iwork_size < 1) iwork_size = 1;
     iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
     CHECK_ALLOC(iwork);

     // Perform workspace query for DWORK
     ldwork = -1; // Query mode
     int ranke_dummy, rnka22_dummy; // Dummy outputs for query
     F77_FUNC(tg01fd, TG01FD)(&compq_upper, &compz_upper, &joba_upper, &l, &n, &m, &p,
                              a_ptr, &lda_f, e_ptr, &lde_f, b_ptr, &ldb_f, c_ptr, &ldc_f,
                              q_ptr, &ldq_f, z_ptr, &ldz_f,
                              &ranke_dummy, &rnka22_dummy, &tol,
                              iwork, &dwork_query, &ldwork, &info,
                              compq_len, compz_len, joba_len);

     if (info < 0 && info != -24) { info = info; goto cleanup; } // Query failed due to invalid argument
     info = 0; // Reset info after query

     // Get the required dwork size from query result
     ldwork = (int)dwork_query;
     // Check against minimum documented size
     int min_ldwork = MAX(1, MAX(n+p, MIN(l,n) + MAX(3*n-1, MAX(m, l))));
     ldwork = MAX(ldwork, min_ldwork);

     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure

     /* --- Call the computational routine --- */
     // RANKE, RNKA22 are scalar outputs, no transposition needed
     F77_FUNC(tg01fd, TG01FD)(&compq_upper, &compz_upper, &joba_upper, &l, &n, &m, &p,
                              a_ptr, &lda_f, e_ptr, &lde_f, b_ptr, &ldb_f, c_ptr, &ldc_f,
                              q_ptr, &ldq_f, z_ptr, &ldz_f,
                              ranke, rnka22, &tol,
                              iwork, dwork, &ldwork, &info,
                              compq_len, compz_len, joba_len);

     /* --- Copy results back to row-major format if needed --- */
     if (row_major && info == 0) {
         size_t a_rows = l; size_t a_cols = n; size_t a_size = a_rows * a_cols;
         size_t e_rows = l; size_t e_cols = n; size_t e_size = e_rows * e_cols;
         size_t b_rows = l; size_t b_cols = m; size_t b_size = b_rows * b_cols;
         size_t c_rows = p; size_t c_cols = n; size_t c_size = c_rows * c_cols;
         size_t q_rows = l; size_t q_cols = l; size_t q_size = q_rows * q_cols;
         size_t z_rows = n; size_t z_cols = n; size_t z_size = z_rows * z_cols;

         if (a_size > 0 && a_cm) slicot_transpose_to_c(a_cm, a, a_rows, a_cols, elem_size);
         if (e_size > 0 && e_cm) slicot_transpose_to_c(e_cm, e, e_rows, e_cols, elem_size);
         if (b_size > 0 && b_cm) slicot_transpose_to_c(b_cm, b, b_rows, b_cols, elem_size);
         if (c_size > 0 && c_cm) slicot_transpose_to_c(c_cm, c, c_rows, c_cols, elem_size);
         if (compq_upper != 'N' && q_size > 0 && q_cm) slicot_transpose_to_c(q_cm, q, q_rows, q_cols, elem_size);
         if (compz_upper != 'N' && z_size > 0 && z_cm) slicot_transpose_to_c(z_cm, z, z_rows, z_cols, elem_size);
         // RANKE, RNKA22 modified directly
     }
     // In column-major case, A, E, B, C, Q, Z, RANKE, RNKA22 modified in place.

  cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(iwork);
     free(a_cm); free(e_cm); free(b_cm); free(c_cm);
     free(q_cm); free(z_cm);

     return info;
 }