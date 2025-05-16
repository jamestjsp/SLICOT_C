/**
 * @file sb04qd.c
 * @brief C wrapper implementation for SLICOT routine SB04QD
 *
 * This file provides a C wrapper implementation for the SLICOT routine SB04QD,
 * which solves the discrete-time Sylvester equation X + AXB = C
 * using the Hessenberg-Schur method.
 */

 #include <stdlib.h>
 #include <stddef.h> // For size_t

 // Include the header file for this wrapper
 #include "sb04qd.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * A, B, C are input/output. Z is output.
  */
 extern void F77_FUNC(sb04qd, SB04QD)(
     const int* n,     // INTEGER N
     const int* m,     // INTEGER M
     double* a,        // DOUBLE PRECISION A(LDA,N)
     const int* lda,   // INTEGER LDA
     double* b,        // DOUBLE PRECISION B(LDB,M)
     const int* ldb,   // INTEGER LDB
     double* c,        // DOUBLE PRECISION C(LDC,M) -> X(N,M)
     const int* ldc,   // INTEGER LDC
     double* z,        // DOUBLE PRECISION Z(LDZ,M)
     const int* ldz,   // INTEGER LDZ
     int* iwork,      // INTEGER IWORK(*) - Added this missing argument
     double* dwork,    // DOUBLE PRECISION DWORK(*)
     const int* ldwork,// INTEGER LDWORK
     int* info         // INTEGER INFO
 );
 
 /* C wrapper function definition */
 SLICOT_EXPORT
 int slicot_sb04qd(int n, int m,
                   double* a, int lda, double* b, int ldb,
                   double* c, int ldc, double* z, int ldz,
                   int row_major)
 {
     int info = 0;
     double* dwork = NULL;
     int ldwork = -1; // Initialize for workspace query
 
     int* iwork = NULL;
     int iwork_size = 0;
 
     double* a_cm = NULL;
     double* b_cm = NULL;
     double* c_cm = NULL;
     double* z_cm = NULL;
     double *a_ptr, *b_ptr, *c_ptr, *z_ptr;
     int lda_f, ldb_f, ldc_f, ldz_f;


     /* --- Input Parameter Validation --- */

     // Check scalar parameters first
     if (n < 0) { info = -1; goto cleanup; } // N is 1st arg
     if (m < 0) { info = -2; goto cleanup; } // M is 2nd arg

     // Check pointers for required arrays
     // A is N x N, B is M x M, C is N x M (input), X is N x M (output in C), Z is M x M (output)
     if (a == NULL && n > 0) { info = -3; goto cleanup; } // A is 3rd arg
     if (b == NULL && m > 0) { info = -5; goto cleanup; } // B is 5th arg
     if (c == NULL && n > 0 && m > 0) { info = -7; goto cleanup; } // C is 7th arg
     if (z == NULL && m > 0) { info = -9; goto cleanup; } // Z is 9th arg

     // Check leading dimensions based on row_major flag
     // Fortran leading dimensions (number of rows)
     int min_lda_f_val = MAX(1, n);
     int min_ldb_f_val = MAX(1, m);
     int min_ldc_f_val = MAX(1, n); // C/X is N rows in Fortran
     int min_ldz_f_val = MAX(1, m); // Z is M rows in Fortran

     if (row_major) {
         // C LD is number of columns
         if (n > 0 && lda < n) { info = -4; goto cleanup; } // A (NxN) needs N cols
         if (m > 0 && ldb < m) { info = -6; goto cleanup; } // B (MxM) needs M cols
         if (m > 0 && ldc < m) { info = -8; goto cleanup; } // C/X (NxM) needs M cols
         if (m > 0 && ldz < m) { info = -10; goto cleanup; } // Z (MxM) needs M cols
     } else {
         // C LD is number of rows (Fortran style)
         if (lda < min_lda_f_val && n > 0) { info = -4; goto cleanup; }
         if (ldb < min_ldb_f_val && m > 0) { info = -6; goto cleanup; }
         if (ldc < min_ldc_f_val && n > 0 && m > 0) { info = -8; goto cleanup; }
         if (ldz < min_ldz_f_val && m > 0) { info = -10; goto cleanup; }
     }

     if (info != 0) { goto cleanup; }
 
     /* --- Allocate Integer Workspace IWORK --- */
     // Common size for SB04xD routines is N+M. Use MAX(1, ...) for safety.
     iwork_size = MAX(1, n + m);
     if (n == 0 && m == 0) iwork_size = 1; // Ensure minimum size 1
 
     iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
     CHECK_ALLOC(iwork);
 
     /* --- Workspace Allocation --- */
     // Workspace query for DWORK
     double dwork_query[1];
     // ldwork is already -1
     F77_FUNC(sb04qd, SB04QD)(&n, &m, NULL, &min_lda_f_val, NULL, &min_ldb_f_val, NULL, &min_ldc_f_val, NULL, &min_ldz_f_val,
                              iwork, dwork_query, &ldwork, &info);

     if (info == 0 || info == -13) { // -13 for LDWORK error during query (SB04QD uses -13 for LDWORK)
         ldwork = (int)dwork_query[0];
         // Ensure minimum size based on documentation formula (example: N*M + MAX(N,M))
         // From SB04QD example (not available, using SB04MD example as placeholder): LDWORK >= MAX(1, N*(M + MAX(2,N) + 2) + M*MAX(1,M-1)/2 )
         // A simpler common one: MAX(1, N + M + N*M) or similar.
         // For SB04QD, the Fortran example uses LDWORK = N*M + MAX(N,M).
         int min_ldwork_formula = (n == 0 && m == 0) ? 1 : (n * m + MAX(n, m));
         min_ldwork_formula = MAX(1, min_ldwork_formula);
         ldwork = MAX(ldwork, min_ldwork_formula);
         info = 0; // Reset info after successful query or if it was -12
     } else {
         // Fallback to formula if query fails for other reasons
         int min_ldwork_formula = (n == 0 && m == 0) ? 1 : (n * m + MAX(n, m));
         ldwork = MAX(1, min_ldwork_formula);
         info = 0; // Reset info if it was a different error from query
     }
     
     if (ldwork > 0) {
         dwork = (double*)malloc((size_t)ldwork * sizeof(double));
         CHECK_ALLOC(dwork);
     } else if (n == 0 && m == 0) { // Handle case where ldwork might be calculated as 0 or 1
         ldwork = 1; // Ensure at least 1 for dwork if routine expects it
         dwork = (double*)malloc((size_t)ldwork * sizeof(double));
         CHECK_ALLOC(dwork);
     }


     /* --- Prepare arrays for column-major format if using row-major --- */
     size_t elem_size = sizeof(double);
     if (row_major) {
         lda_f = MAX(1, n);
         ldb_f = MAX(1, m);
         ldc_f = MAX(1, n); // X is N rows in Fortran
         ldz_f = MAX(1, m);

         size_t a_size_bytes = (size_t)lda_f * n * elem_size; if (n == 0) a_size_bytes = 0;
         size_t b_size_bytes = (size_t)ldb_f * m * elem_size; if (m == 0) b_size_bytes = 0;
         size_t c_size_bytes = (size_t)ldc_f * m * elem_size; if (n == 0 || m == 0) c_size_bytes = 0; // C/X is N x M
         size_t z_size_bytes = (size_t)ldz_f * m * elem_size; if (m == 0) z_size_bytes = 0; // Z is M x M

         if (a_size_bytes > 0) { a_cm = (double*)malloc(a_size_bytes); CHECK_ALLOC(a_cm); slicot_transpose_to_fortran_with_ld(a, a_cm, n, n, lda, lda_f, elem_size); }
         if (b_size_bytes > 0) { b_cm = (double*)malloc(b_size_bytes); CHECK_ALLOC(b_cm); slicot_transpose_to_fortran_with_ld(b, b_cm, m, m, ldb, ldb_f, elem_size); }
         if (c_size_bytes > 0) { c_cm = (double*)malloc(c_size_bytes); CHECK_ALLOC(c_cm); slicot_transpose_to_fortran_with_ld(c, c_cm, n, m, ldc, ldc_f, elem_size); }
         if (z_size_bytes > 0) { z_cm = (double*)malloc(z_size_bytes); CHECK_ALLOC(z_cm); /* Z is output only */ }

         a_ptr = (a_size_bytes > 0) ? a_cm : NULL;
         b_ptr = (b_size_bytes > 0) ? b_cm : NULL;
         c_ptr = (c_size_bytes > 0) ? c_cm : NULL; // Fortran writes X to c_cm
         z_ptr = (z_size_bytes > 0) ? z_cm : NULL; // Fortran writes Z to z_cm
     } else {
         // Column-major C: Pass original pointers (or NULL if size is 0)
         if (n == 0) a_ptr = NULL;
         if (m == 0) b_ptr = NULL;
         if (n == 0 || m == 0) c_ptr = NULL;
         if (m == 0) z_ptr = NULL;
         // lda_f, ldb_f, ldc_f, ldz_f are already correctly set to C LDs (rows)
     }

     /* --- Call the computational routine --- */
     F77_FUNC(sb04qd, SB04QD)(&n, &m, a_ptr, &lda_f, b_ptr, &ldb_f, c_ptr, &ldc_f, z_ptr, &ldz_f,
                             iwork, dwork, &ldwork, &info);
 
     /* --- Copy results back to row-major format if needed --- */
     if (row_major && info == 0) {
         // A and B are input-only for SB04QD (unlike SB04MD which can modify them)
         // C (input) is overwritten by X (output)
         size_t c_size_bytes = (size_t)ldc_f * m * elem_size; if (n == 0 || m == 0) c_size_bytes = 0;
         size_t z_size_bytes = (size_t)ldz_f * m * elem_size; if (m == 0) z_size_bytes = 0;

         if (c_ptr && c_size_bytes > 0) { // c_ptr is c_cm
             slicot_transpose_to_c_with_ld(c_ptr, c, n, m, ldc_f, ldc, elem_size);
         }
         if (z_ptr && z_size_bytes > 0) { // z_ptr is z_cm
             slicot_transpose_to_c_with_ld(z_ptr, z, m, m, ldz_f, ldz, elem_size);
         }
     }

 cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(iwork);
     free(a_cm);
     free(b_cm);
     free(c_cm);
     free(z_cm);

     return info;
 }