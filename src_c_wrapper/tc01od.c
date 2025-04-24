/**
 * @file tc01od.c
 * @brief C wrapper implementation for SLICOT routine TC01OD
 *
 * This file provides a C wrapper implementation for the SLICOT routine TC01OD,
 * which finds the dual right (left) polynomial matrix representation
 * of a given left (right) polynomial matrix representation by transposing
 * the coefficient matrices.
 */

 #include <stdlib.h>
 #include <ctype.h>
 #include <string.h> // For memcpy
 #include <stddef.h> // For size_t

 // Include the header file for this wrapper
 // #include "tc01od.h" // Assuming a header file exists
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, MIN, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * Note the 3D arrays PCOEFF, QCOEFF are passed as flat pointers.
  * Hidden length for CHARACTER argument is added at the end.
  */
 extern void F77_FUNC(tc01od, TC01OD)(
     const char* leri,       // CHARACTER*1 LERI
     const int* m,           // INTEGER M
     const int* p,           // INTEGER P
     const int* indlim,      // INTEGER INDLIM
     double* pcoeff,         // DOUBLE PRECISION PCOEFF(LDPCO1,LDPCO2,*) (in/out)
     const int* ldpco1,      // INTEGER LDPCO1
     const int* ldpco2,      // INTEGER LDPCO2
     double* qcoeff,         // DOUBLE PRECISION QCOEFF(LDQCO1,LDQCO2,*) (in/out) - Needs workspace
     const int* ldqco1,      // INTEGER LDQCO1
     const int* ldqco2,      // INTEGER LDQCO2
     int* info,              // INTEGER INFO (output)
     int leri_len            // Hidden length
 );


 /* C wrapper function definition */
 SLICOT_EXPORT
 int slicot_tc01od(char leri, int m, int p, int indlim,
                   double* pcoeff, int ldpco1, int ldpco2,
                   double* qcoeff, int ldqco1, int ldqco2,
                   int row_major)
 {
     /* Local variables */
     int info = 0;
     // No workspace needed

     const int leri_len = 1;
     char leri_upper = toupper(leri);

     /* Pointers for column-major copies if needed */
     double *pcoeff_cm = NULL, *qcoeff_cm = NULL;
     double *pcoeff_ptr, *qcoeff_ptr;
     int ldpco1_f, ldpco2_f, ldqco1_f, ldqco2_f;

     /* Determine dimensions based on LERI */
     int porm = (leri_upper == 'L') ? p : m; // Input P(s) dimension
     int porp = (leri_upper == 'L') ? m : p; // Input Q(s) other dimension
     int maxmp = MAX(m, p);

     /* --- Input Parameter Validation --- */
     if (leri_upper != 'L' && leri_upper != 'R') { info = -1; goto cleanup; }
     if (m < 0) { info = -2; goto cleanup; }
     if (p < 0) { info = -3; goto cleanup; }
     if (indlim < 1) { info = -4; goto cleanup; }

     // Check leading dimensions based on storage order
     // Note: Fortran workspace requirements may mandate larger LDs than input/output sizes.
     int min_ldpco1_f = MAX(1, porm);
     int min_ldpco2_f = MAX(1, porm);
     int min_ldqco1_f = MAX(1, maxmp); // Needs workspace (size MAXMP x MAXMP)
     int min_ldqco2_f = MAX(1, maxmp); // Needs workspace (size MAXMP x MAXMP)

     if (row_major) {
         // For row-major C, LD is number of columns (second dim for 3D slice)
         int min_ldpco1_rm_rows = porm;
         int min_ldpco2_rm_cols = porm;
         int min_ldqco1_rm_rows = maxmp; // Needs workspace
         int min_ldqco2_rm_cols = maxmp; // Needs workspace

         if (ldpco1 < min_ldpco1_rm_rows) { info = -6; goto cleanup; }
         if (ldpco2 < min_ldpco2_rm_cols) { info = -7; goto cleanup; }
         if (ldqco1 < min_ldqco1_rm_rows) { info = -9; goto cleanup; }
         if (ldqco2 < min_ldqco2_rm_cols) { info = -10; goto cleanup; }
     } else {
         // For column-major C, LD is number of rows (Fortran style)
         if (ldpco1 < min_ldpco1_f) { info = -6; goto cleanup; }
         if (ldpco2 < min_ldpco2_f) { info = -7; goto cleanup; }
         if (ldqco1 < min_ldqco1_f) { info = -9; goto cleanup; }
         if (ldqco2 < min_ldqco2_f) { info = -10; goto cleanup; }
     }

     /* --- Workspace Allocation --- */
     // No workspace needed for TC01OD

     /* --- Prepare Arrays and Call Fortran Routine --- */
     size_t elem_size = sizeof(double);
     if (row_major) {
         /* Calculate total sizes for allocation/transposition */
         size_t pcoeff_slice_size = (size_t)porm * porm;
         size_t qcoeff_slice_size = (size_t)maxmp * maxmp; // Use workspace size for copy
         size_t pcoeff_total_size = pcoeff_slice_size * indlim;
         size_t qcoeff_total_size = qcoeff_slice_size * indlim;

         /* Allocate memory for column-major copies */
         if (pcoeff_total_size > 0) { pcoeff_cm = (double*)malloc(pcoeff_total_size * elem_size); CHECK_ALLOC(pcoeff_cm); }
         if (qcoeff_total_size > 0) { qcoeff_cm = (double*)malloc(qcoeff_total_size * elem_size); CHECK_ALLOC(qcoeff_cm); }

         /* Transpose each 2D slice from C (row-major) to Fortran (column-major) */
         // Determine input dimensions for Q
         int q_in_rows = p;
         int q_in_cols = m;
         size_t q_in_slice_size = (size_t)q_in_rows * q_in_cols;

         for (int k = 0; k < indlim; ++k) {
             if (pcoeff_slice_size > 0 && pcoeff_cm) {
                 slicot_transpose_to_fortran(pcoeff + k * pcoeff_slice_size, pcoeff_cm + k * pcoeff_slice_size, porm, porm, elem_size);
             }
             if (q_in_slice_size > 0 && qcoeff_cm) {
                 // Transpose only the PxM part into the MAXMP x MAXMP buffer
                 slicot_transpose_to_fortran(qcoeff + k * q_in_slice_size, qcoeff_cm + k * qcoeff_slice_size, q_in_rows, q_in_cols, elem_size);
             }
         }

         /* Fortran leading dimensions (use workspace sizes) */
         ldpco1_f = (porm > 0) ? porm : 1;
         ldpco2_f = (porm > 0) ? porm : 1;
         ldqco1_f = (maxmp > 0) ? maxmp : 1;
         ldqco2_f = (maxmp > 0) ? maxmp : 1;

         /* Set pointers */
         pcoeff_ptr = pcoeff_cm;
         qcoeff_ptr = qcoeff_cm;

     } else {
         /* Column-major case - use original arrays */
         ldpco1_f = ldpco1; ldpco2_f = ldpco2;
         ldqco1_f = ldqco1; ldqco2_f = ldqco2;
         pcoeff_ptr = pcoeff;
         qcoeff_ptr = qcoeff;
     }

     /* --- Call the computational routine --- */
     F77_FUNC(tc01od, TC01OD)(&leri_upper, &m, &p, &indlim,
                              pcoeff_ptr, &ldpco1_f, &ldpco2_f,
                              qcoeff_ptr, &ldqco1_f, &ldqco2_f, &info,
                              leri_len);

     /* --- Copy results back to row-major format if needed --- */
     if (row_major && info == 0) {
         // Determine output dimensions
         int porm_out = (leri_upper == 'L') ? m : p; // Output P'(s) dimension
         int q_rows_out = m; // Output Q'(s) is M x P
         int q_cols_out = p;
         size_t pcoeff_out_slice_size = (size_t)porm_out * porm_out;
         size_t qcoeff_out_slice_size = (size_t)q_rows_out * q_cols_out;

         // Determine CM buffer slice sizes used by Fortran
         size_t pcoeff_cm_slice_size = (size_t)porm * porm;
         size_t qcoeff_cm_slice_size = (size_t)maxmp * maxmp;

         for (int k = 0; k < indlim; ++k) {
             if (pcoeff_out_slice_size > 0 && pcoeff_cm) {
                 // Transpose from the CM buffer (which now holds P') to the original RM buffer
                 slicot_transpose_to_c(pcoeff_cm + k * pcoeff_cm_slice_size, pcoeff + k * pcoeff_out_slice_size, porm_out, porm_out, elem_size);
             }
             if (qcoeff_out_slice_size > 0 && qcoeff_cm) {
                 // Transpose from the CM buffer (which now holds Q') to the original RM buffer
                 slicot_transpose_to_c(qcoeff_cm + k * qcoeff_cm_slice_size, qcoeff + k * qcoeff_out_slice_size, q_rows_out, q_cols_out, elem_size);
             }
         }
     }
     // In column-major case, PCOEFF, QCOEFF modified in place.

  cleanup:
     /* --- Cleanup --- */
     free(pcoeff_cm);
     free(qcoeff_cm);

     return info;
 }