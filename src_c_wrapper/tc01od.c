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
 #include "tc01od.h"
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
     double* qcoeff,         // DOUBLE PRECISION QCOEFF(LDQCO1,LDQCO2,*) (in/out)
     const int* ldqco1,      // INTEGER LDQCO1
     const int* ldqco2,      // INTEGER LDQCO2
     int* info,              // INTEGER INFO (output)
     int leri_len            // Hidden length
 );
 
 
 /* C wrapper function definition */
 int slicot_tc01od(char leri, int m, int p, int indlim,
                   double* pcoeff, int ldpco1, int ldpco2,
                   double* qcoeff, int ldqco1, int ldqco2,
                   int row_major)
 {
     /* Local variables */
     int info = 0;
     const int leri_len = 1;
     char leri_upper = toupper(leri);
 
     /* Pointers for column-major copies if needed */
     double *pcoeff_cm = NULL, *qcoeff_cm = NULL;
 
     /* Determine dimensions based on LERI */
     int porm = (leri_upper == 'L') ? p : m; // Input P(s) dimension
     int porp = (leri_upper == 'L') ? m : p; // Input Q(s) other dimension
     int maxmp = MAX(m, p);
 
     /* --- Input Parameter Validation --- */
     if (m < 0) { info = -2; goto cleanup; }
     if (p < 0) { info = -3; goto cleanup; }
     if (indlim < 1) { info = -4; goto cleanup; }
     if (leri_upper != 'L' && leri_upper != 'R') { info = -1; goto cleanup; }
 
     // Check leading dimensions based on storage order
     int min_ldpco1_f = MAX(1, porm);
     int min_ldpco2_f = MAX(1, porm);
     int min_ldqco1_f = MAX(1, maxmp); // Fortran doc implies PxM always? Check example. Example uses MAXMP.
     int min_ldqco2_f = MAX(1, maxmp); // Fortran doc implies PxM always? Check example. Example uses MAXMP.
 
     if (row_major) {
         // For row-major C, LD is number of columns (second dim for 3D slice)
         int min_ldpco1_rm_rows = porm; // 1st dim is rows
         int min_ldpco2_rm_cols = porm; // 2nd dim is cols
         int min_ldqco1_rm_rows = maxmp; // Use maxmp based on example
         int min_ldqco2_rm_cols = maxmp; // Use maxmp based on example
 
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
 
     // Calculate total sizes for allocation/transposition
     size_t pcoeff_slice_size = (size_t)porm * porm;
     size_t qcoeff_slice_size = (size_t)p * m; // Q is always P x M on input? Example uses MAXMP x MAXMP. Let's assume P x M based on docs.
     size_t pcoeff_total_size = pcoeff_slice_size * indlim;
     size_t qcoeff_total_size = qcoeff_slice_size * indlim;
 
     // Dimensions for Q based on Fortran interface (PxM)
     int q_rows_f = p;
     int q_cols_f = m;
 
 
     if (row_major) {
         /* --- Row-Major Case --- */
         /* Allocate memory for column-major copies */
         if (pcoeff_total_size > 0) { pcoeff_cm = (double*)malloc(pcoeff_total_size * elem_size); CHECK_ALLOC(pcoeff_cm); }
         if (qcoeff_total_size > 0) { qcoeff_cm = (double*)malloc(qcoeff_total_size * elem_size); CHECK_ALLOC(qcoeff_cm); }
 
         /* Transpose each 2D slice from C (row-major) to Fortran (column-major) */
         for (int k = 0; k < indlim; ++k) {
             size_t offset = k * pcoeff_slice_size;
             if (pcoeff_slice_size > 0) {
                 slicot_transpose_to_fortran(pcoeff + offset, pcoeff_cm + offset, porm, porm, elem_size);
             }
             offset = k * qcoeff_slice_size;
              if (qcoeff_slice_size > 0) {
                slicot_transpose_to_fortran(qcoeff + offset, qcoeff_cm + offset, q_rows_f, q_cols_f, elem_size);
             }
         }
 
         /* Fortran leading dimensions */
         int ldpco1_f = (porm > 0) ? porm : 1;
         int ldpco2_f = (porm > 0) ? porm : 1;
         int ldqco1_f = (q_rows_f > 0) ? q_rows_f : 1; // Use actual P
         int ldqco2_f = (q_rows_f > 0) ? q_rows_f : 1; // Use actual P for 2nd dim? Example uses MAXMP. Let's stick to P.
 
         /* Call the Fortran routine */
         F77_FUNC(tc01od, TC01OD)(&leri_upper, &m, &p, &indlim,
                                  pcoeff_cm, &ldpco1_f, &ldpco2_f,
                                  qcoeff_cm, &ldqco1_f, &ldqco2_f, &info,
                                  leri_len);
 
         /* Copy back results from column-major temps to original row-major arrays */
         if (info == 0) {
             // Determine output dimensions
             int porm_out = (leri_upper == 'L') ? m : p; // Output P'(s) dimension
             int q_rows_out = m; // Output Q'(s) is M x P
             int q_cols_out = p;
             size_t pcoeff_out_slice_size = (size_t)porm_out * porm_out;
             size_t qcoeff_out_slice_size = (size_t)q_rows_out * q_cols_out;
 
             for (int k = 0; k < indlim; ++k) {
                 size_t offset_p_cm = k * pcoeff_slice_size; // Offset in potentially padded input CM buffer
                 size_t offset_q_cm = k * qcoeff_slice_size;
                 size_t offset_p_rm = k * pcoeff_out_slice_size; // Offset in potentially smaller output RM buffer
                 size_t offset_q_rm = k * qcoeff_out_slice_size;
 
                 if (pcoeff_slice_size > 0 && pcoeff_out_slice_size > 0) {
                     // Transpose from the CM buffer (which now holds P') to the original RM buffer
                     slicot_transpose_to_c(pcoeff_cm + offset_p_cm, pcoeff + offset_p_rm, porm_out, porm_out, elem_size);
                 }
                  if (qcoeff_slice_size > 0 && qcoeff_out_slice_size > 0) {
                     // Transpose from the CM buffer (which now holds Q') to the original RM buffer
                     slicot_transpose_to_c(qcoeff_cm + offset_q_cm, qcoeff + offset_q_rm, q_rows_out, q_cols_out, elem_size);
                 }
             }
         }
         /* Temps freed in cleanup */
 
     } else {
         /* --- Column-Major Case --- */
         /* Call the Fortran routine directly with user-provided arrays */
         F77_FUNC(tc01od, TC01OD)(&leri_upper, &m, &p, &indlim,
                                  pcoeff, &ldpco1, &ldpco2,
                                  qcoeff, &ldqco1, &ldqco2, &info,
                                  leri_len);
         // PCOEFF, QCOEFF modified in place.
     }
 
  cleanup:
     /* --- Cleanup --- */
     free(pcoeff_cm);
     free(qcoeff_cm);
 
     return info;
 }
 