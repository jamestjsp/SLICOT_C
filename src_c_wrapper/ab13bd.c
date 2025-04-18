/**
 * @file ab13bd.c
 * @brief C wrapper implementation for SLICOT function AB13BD
 *
 * This file provides a C wrapper implementation for the SLICOT function AB13BD,
 * which computes the H2 or L2 norm of a system (A,B,C,D).
 */

 #include <stdlib.h>
 #include <ctype.h>
 #include <stddef.h> // For size_t
 
 // Include the header file for this wrapper
 #include "ab13bd.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions
 
 /*
  * Declare the external Fortran FUNCTION using the F77_FUNC macro.
  * Note A, B, C, D are input/output. NQ, IWARN, INFO are output.
  */
 extern double F77_FUNC(ab13bd, AB13BD)(
     const char* dico,       // CHARACTER*1 DICO
     const char* jobn,       // CHARACTER*1 JOBN
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     const int* p,           // INTEGER P
     double* a,              // DOUBLE PRECISION A(LDA,*) (in/out)
     const int* lda,         // INTEGER LDA
     double* b,              // DOUBLE PRECISION B(LDB,*) (in/out)
     const int* ldb,         // INTEGER LDB
     double* c,              // DOUBLE PRECISION C(LDC,*) (in/out)
     const int* ldc,         // INTEGER LDC
     double* d,              // DOUBLE PRECISION D(LDD,*) (in/out)
     const int* ldd,         // INTEGER LDD
     int* nq,                // INTEGER NQ (output)
     const double* tol,      // DOUBLE PRECISION TOL
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* iwarn,             // INTEGER IWARN (output)
     int* info,              // INTEGER INFO (output)
     int dico_len,           // Hidden length
     int jobn_len            // Hidden length
 );
 
 
 /* C wrapper function definition */
 SLICOT_C_WRAPPER_API
 double slicot_ab13bd(char dico, char jobn, int n, int m, int p,
                      double* a, int lda, double* b, int ldb,
                      double* c, int ldc, double* d, int ldd,
                      int* nq, double tol, int* iwarn, int* info_ptr,
                      int row_major)
 {
     /* Local variables */
     int info = 0; // Local info code
     int ldwork = -1; /* Use -1 for workspace query */
     double dwork_query;
     double* dwork = NULL;
     const int dico_len = 1, jobn_len = 1;
     double norm_value = 0.0; // Function return value
 
     char dico_upper = toupper(dico);
     char jobn_upper = toupper(jobn);
 
     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL;
 
     /* --- Input Parameter Validation --- */
 
     if (n < 0) { info = -3; goto cleanup; }
     if (m < 0) { info = -4; goto cleanup; }
     if (p < 0) { info = -5; goto cleanup; }
     if (dico_upper != 'C' && dico_upper != 'D') { info = -1; goto cleanup; }
     if (jobn_upper != 'H' && jobn_upper != 'L') { info = -2; goto cleanup; }
 
     // Check leading dimensions based on storage order
     int min_lda_f = MAX(1, n);
     int min_ldb_f = MAX(1, n);
     int min_ldc_f = MAX(1, p);
     int min_ldd_f = MAX(1, p);
 
     if (row_major) {
         // For row-major C, LDA is the number of columns
         int min_lda_rm_cols = n;
         int min_ldb_rm_cols = m;
         int min_ldc_rm_cols = n;
         int min_ldd_rm_cols = m;
         if (lda < min_lda_rm_cols) { info = -7; goto cleanup; }
         if (ldb < min_ldb_rm_cols) { info = -9; goto cleanup; }
         if (ldc < min_ldc_rm_cols) { info = -11; goto cleanup; }
         if (ldd < min_ldd_rm_cols) { info = -13; goto cleanup; }
     } else {
         // For column-major C, LDA is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -7; goto cleanup; }
         if (ldb < min_ldb_f) { info = -9; goto cleanup; }
         if (ldc < min_ldc_f) { info = -11; goto cleanup; }
         if (ldd < min_ldd_f) { info = -13; goto cleanup; }
     }
 
     /* --- Workspace Allocation --- */
 
     // Allocate DWORK based on query
     ldwork = -1; // Query mode
     // Call function just for query - result is irrelevant here
     (void) F77_FUNC(ab13bd, AB13BD)(&dico_upper, &jobn_upper, &n, &m, &p,
                                     a, &lda, b, &ldb, c, &ldc, d, &ldd,
                                     nq, &tol, &dwork_query, &ldwork,
                                     iwarn, &info,
                                     dico_len, jobn_len);
 
     if (info < 0) { goto cleanup; } // Query failed due to invalid argument
     info = 0; // Reset info after query
 
     // Get the required dwork size from query result
     ldwork = (int)dwork_query;
     // Check against minimum documented size:
     // MAX( 1, M*(N+M) + MAX( N*(N+5), M*(M+2), 4*P ), N*( MAX( N, P ) + 4 ) + MIN( N, P ) )
     int min_ldwork = 1;
     int term1 = m * (n + m);
     int term2 = MAX(n * (n + 5), MAX(m * (m + 2), 4 * p));
     int term3 = n * (MAX(n, p) + 4) + MIN(n, p);
     min_ldwork = MAX(min_ldwork, MAX(term1 + term2, term3));
     ldwork = MAX(ldwork, min_ldwork);
 
     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure
 
     /* --- Prepare Arrays and Call Fortran Routine --- */
     size_t elem_size = sizeof(double);
 
     if (row_major) {
         /* --- Row-Major Case --- */
 
         /* Allocate memory for column-major copies */
         size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
         size_t b_rows = n; size_t b_cols = m; size_t b_size = b_rows * b_cols;
         size_t c_rows = p; size_t c_cols = n; size_t c_size = c_rows * c_cols;
         size_t d_rows = p; size_t d_cols = m; size_t d_size = d_rows * d_cols;
 
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
         if (d_size > 0) { d_cm = (double*)malloc(d_size * elem_size); CHECK_ALLOC(d_cm); }
 
         /* Transpose C (row-major) inputs to Fortran (column-major) copies */
         if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         if (b_size > 0) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, elem_size);
         if (c_size > 0) slicot_transpose_to_fortran(c, c_cm, c_rows, c_cols, elem_size);
         if (d_size > 0) slicot_transpose_to_fortran(d, d_cm, d_rows, d_cols, elem_size);
 
         /* Fortran leading dimensions */
         int lda_f = (a_rows > 0) ? a_rows : 1;
         int ldb_f = (b_rows > 0) ? b_rows : 1;
         int ldc_f = (c_rows > 0) ? c_rows : 1;
         int ldd_f = (d_rows > 0) ? d_rows : 1;
 
         /* Call the Fortran FUNCTION */
         norm_value = F77_FUNC(ab13bd, AB13BD)(&dico_upper, &jobn_upper, &n, &m, &p,
                                               a_cm, &lda_f, b_cm, &ldb_f,
                                               c_cm, &ldc_f, d_cm, &ldd_f,
                                               nq, &tol, dwork, &ldwork,
                                               iwarn, &info,
                                               dico_len, jobn_len);
 
         /* Copy back modified arrays from column-major temps to original row-major arrays */
         // Copy back even if info > 0, as results might still be useful
         if (info >= 0) {
             // NQ is output, use it for dimensions where applicable
             int nq_val = (nq != NULL) ? *nq : n; // Use n if nq is NULL (should not happen)
             // A is NQxNQ, B is NQxM, C is PxNQ, D is PxM after transformation
             if (a_size > 0) slicot_transpose_to_c(a_cm, a, nq_val, nq_val, elem_size);
             if (b_size > 0) slicot_transpose_to_c(b_cm, b, nq_val, m, elem_size);
             if (c_size > 0) slicot_transpose_to_c(c_cm, c, p, nq_val, elem_size);
             if (d_size > 0) slicot_transpose_to_c(d_cm, d, p, m, elem_size);
         }
         /* Column-major copies will be freed in cleanup */
 
     } else {
         /* --- Column-Major Case --- */
 
         /* Call the Fortran FUNCTION directly with user-provided arrays */
         norm_value = F77_FUNC(ab13bd, AB13BD)(&dico_upper, &jobn_upper, &n, &m, &p,
                                               a, &lda, b, &ldb, c, &ldc, d, &ldd,
                                               nq, &tol, dwork, &ldwork,
                                               iwarn, &info,
                                               dico_len, jobn_len);
         // A, B, C, D are modified in place by the Fortran call.
     }
 
 cleanup:
     /* --- Cleanup --- */
     free(dwork);
     // No iwork for this routine
     free(a_cm);
     free(b_cm);
     free(c_cm);
     free(d_cm);
 
     // Store info code in user-provided pointer
     if (info_ptr) *info_ptr = info;
 
     // Return the computed norm (or 0.0 if error occurred before call)
     return (info >= 0) ? norm_value : 0.0;
 }
 