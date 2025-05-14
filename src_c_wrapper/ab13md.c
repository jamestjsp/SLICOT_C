/**
 * @file ab13md.c
 * @brief C wrapper implementation for SLICOT routine AB13MD
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB13MD,
 * which computes an upper bound on the structured singular value for a
 * square complex matrix Z with a given block uncertainty structure.
 */

 #include <stdlib.h>
 #include <ctype.h>  // For toupper
 #include <stddef.h> // For size_t
 #include <string.h> // For memset
 
 // Include the header file for this wrapper
 #include "ab13md.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" 
 #include "slicot_f77.h"   
 
 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  */
 extern void F77_FUNC(ab13md, AB13MD)(
     const char* fact,       // CHARACTER*1 FACT
     const int* n,           // INTEGER N
     const slicot_complex_double* z, // COMPLEX*16 Z(LDZ,*)
     const int* ldz,         // INTEGER LDZ
     const int* m,           // INTEGER M
     const int* nblock,      // INTEGER NBLOCK(*)
     const int* itype,       // INTEGER ITYPE(*)
     double* x,              // DOUBLE PRECISION X(*) (in/out)
     double* bound,          // DOUBLE PRECISION BOUND (output)
     double* d,              // DOUBLE PRECISION D(*) (output)
     double* g,              // DOUBLE PRECISION G(*) (output)
     int* iwork,             // INTEGER IWORK(*)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     slicot_complex_double* zwork, // COMPLEX*16 ZWORK(*)
     const int* lzwork,      // INTEGER LZWORK
     int* info,              // INTEGER INFO (output)
     int fact_len            // Hidden length
 );
 
 
 /* C wrapper function definition */
 SLICOT_EXPORT
 int slicot_ab13md(char fact, int n, const slicot_complex_double* z, int ldz,
                   int m, const int* nblock, const int* itype,
                   double* x, double* bound, double* d, double* g,
                   int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = 1;
     int lzwork = 1;
     double* dwork = NULL;
     slicot_complex_double* zwork = NULL;
     int* iwork = NULL;
     int iwork_size = 0;
     const int fact_len = 1;
 
     char fact_upper = toupper(fact);
 
     /* Pointers for column-major copies if needed */
     slicot_complex_double *z_cm = NULL;
 
     /* Pointers to pass to Fortran */
     const slicot_complex_double *z_ptr;
     int ldz_f;
 
     /* --- Input Parameter Validation (Essential C-level checks) --- */
 
     if (fact_upper != 'F' && fact_upper != 'N') { info = -1; goto cleanup; } // FACT
     if (n < 0) { info = -2; goto cleanup; } // N
     // Z (arg 3)
     if (n > 0 && z == NULL) { info = -3; goto cleanup; }
     
     // LDZ (arg 4)
     int min_ldz_f_val = MAX(1, n);
     if (row_major) {
         if (n > 0 && ldz < n) { info = -4; goto cleanup; } // For row-major C, LDZ is #cols
     } else {
         if (n > 0 && ldz < min_ldz_f_val) { info = -4; goto cleanup; } // For col-major C, LDZ is #rows
     }
 
     if (m < 1) { info = -5; goto cleanup; } // M
     // NBLOCK (arg 6)
     if (m > 0 && nblock == NULL) { info = -6; goto cleanup; }
     // ITYPE (arg 7)
     if (m > 0 && itype == NULL) { info = -7; goto cleanup; }
     
     // X (arg 8) - refined validation
     // X is used if FACT='F' or if NBLOCK[0] != N (and its dimension is > 0).
     // If X is used and is NULL, then it's an error.
     if (m > 0 && nblock != NULL && itype != NULL) { // Need nblock & itype to determine if X is used and its size
         int first_block_size = nblock[0]; // Assuming m > 0, nblock is not NULL
         if (fact_upper == 'F' || (n > 0 && first_block_size != n)) { // Conditions under which X is actively used as input or state
             int num_real_blocks = 0;
             for (int i_mr = 0; i_mr < m; ++i_mr) {
                 if (itype[i_mr] == 1) num_real_blocks++;
             }
             int x_dim_calculated = m + num_real_blocks - 1;
             if (x_dim_calculated > 0 && x == NULL) {
                 info = -8; goto cleanup;
             }
         }
     } else if (m > 0 && (nblock == NULL || itype == NULL)) {
         // This case should have been caught by -6 or -7, but as a safeguard for X logic:
         // If nblock or itype is NULL (and m>0), we can't determine X properties.
         // If FACT='F', X is required.
         if (fact_upper == 'F' && x == NULL) { info = -8; goto cleanup; }
     }
 
 
     if (bound == NULL) { info = -9; goto cleanup; } // BOUND
     // D (arg 10)
     if (n > 0 && d == NULL) { info = -10; goto cleanup; }
     // G (arg 11)
     if (n > 0 && g == NULL) { info = -11; goto cleanup; }
 
     /* --- Workspace Allocation --- */
 
     // IWORK: size MAX(1, MAX(4*M-2, N))
     // Ensure m is not used in MAX if m is 0, but m is already validated >=1
     iwork_size = MAX(1, MAX(4 * m - 2, n));
     iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
     CHECK_ALLOC(iwork);
 
     // DWORK and ZWORK
     if (n == 0) { // M must be >= 1
         ldwork = MAX(1, 9*m*m + 33*m - 11);
         lzwork = MAX(1, 6*m - 3);
     } else { // N > 0 and M >= 1
         // Workspace query
         int info_query = 0;
         double dwork_query_val[1];
         slicot_complex_double zwork_query_val[1];
         int ldwork_for_query = -1;
         int lzwork_for_query = -1;
         
         // For query, pass NULL for Z, X, D, G if they are outputs or their content isn't needed for sizing
         // NBLOCK and ITYPE are needed for structure.
         // LDZ for query must be valid.
         int ldz_for_query = row_major ? MAX(1, n) : (n > 0 ? ldz : 1);
 
 
         F77_FUNC(ab13md, AB13MD)(&fact_upper, &n,
                                  NULL, &ldz_for_query, // Z is NULL for query
                                  &m, nblock, itype,    // NBLOCK, ITYPE are needed
                                  NULL, NULL, NULL, NULL, // X, BOUND, D, G are NULL for query
                                  iwork,                // IWORK might be used by some query paths
                                  dwork_query_val, &ldwork_for_query,
                                  zwork_query_val, &lzwork_for_query, 
                                  &info_query,
                                  fact_len);
         
         // Minimum workspace sizes from documentation formula
         int min_ldwork_formula = MAX(1, 2*n*n*m - n*n + 9*m*m + n*m + 11*n + 33*m - 11);
         int min_lzwork_formula = MAX(1, 6*n*n*m + 12*n*n + 6*m + 6*n - 3);
 
         if (info_query == 0) { // Query successful
             ldwork = MAX((int)dwork_query_val[0], min_ldwork_formula);
             lzwork = MAX((int)SLICOT_COMPLEX_REAL(zwork_query_val[0]), min_lzwork_formula);
         } else { // Query failed or returned unexpected info
             ldwork = min_ldwork_formula;
             lzwork = min_lzwork_formula;
         }
     }
     ldwork = MAX(1, ldwork); // Ensure at least 1
     lzwork = MAX(1, lzwork); // Ensure at least 1
 
     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork);
     zwork = (slicot_complex_double*)malloc((size_t)lzwork * sizeof(slicot_complex_double));
     CHECK_ALLOC(zwork);
     
     // Initialize workspaces (optional, but good practice)
     if (iwork) memset(iwork, 0, (size_t)iwork_size * sizeof(int));
     if (dwork) memset(dwork, 0, (size_t)ldwork * sizeof(double));
     if (zwork) memset(zwork, 0, (size_t)lzwork * sizeof(slicot_complex_double));
 
 
     /* --- Prepare Z Array for Fortran --- */
     size_t z_elem_size = sizeof(slicot_complex_double);
     if (row_major && n > 0) {
         ldz_f = MAX(1, n); // Fortran LD is number of rows
         size_t z_matrix_elements = (size_t)n * n;
         if (z_matrix_elements > 0 && z != NULL) { // z should not be NULL if n > 0 (checked earlier)
             z_cm = (slicot_complex_double*)malloc(z_matrix_elements * z_elem_size);
             CHECK_ALLOC(z_cm);
             // For row_major C, ldz is number of columns.
             slicot_transpose_to_fortran_with_ld(z, z_cm, n, n, ldz, ldz_f, z_elem_size);
             z_ptr = z_cm;
         } else { // Should not happen if n > 0 due to prior checks
             z_ptr = NULL; 
         }
     } else { // Column-major C
         ldz_f = (n > 0) ? ldz : 1; // Use C ldz (rows); if n=0, ldz_f must be >=1
         z_ptr = (n > 0 && z != NULL) ? z : NULL;
     }
     if (n == 0) ldz_f = 1; // Ensure LDZ_F is at least 1 for Fortran call if N=0
 
     /* Call the computational Fortran routine */
     F77_FUNC(ab13md, AB13MD)(&fact_upper, &n,
                              z_ptr, &ldz_f,
                              &m, nblock, itype,
                              x, bound, d, g, iwork,
                              dwork, &ldwork, zwork, &lzwork, &info,
                              fact_len);
 
 cleanup:
     /* --- Cleanup --- */
     free(zwork);
     free(dwork);
     free(iwork);
     free(z_cm); 
 
     return info;
 }
 