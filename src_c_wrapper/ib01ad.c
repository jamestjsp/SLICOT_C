/**
 * @file ib01ad.c
 * @brief C wrapper implementation for SLICOT routine IB01AD
 *
 * This file provides a C wrapper implementation for the SLICOT routine IB01AD,
 * which estimates the system order and computes the triangular factor
 * of the concatenated block-Hankel matrices for subsequent
 * system identification in the SLICOT IB01 family of routines.
 *
 * @note This version allocates IWORK and DWORK internally.
 */

 #include <stdlib.h>
 #include <ctype.h>
 #include <stddef.h>
 #include <string.h> // For memset, memcpy
 #include <stdbool.h> // For bool type
 #include <math.h>    // For isnan, isinf (optional)
 #include <stdio.h>   // For error messages if needed
 
 #include "ib01ad.h"
 #include "slicot_utils.h" // Provides MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR etc.
 #include "slicot_f77.h"
 
 
 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  */
 extern void F77_FUNC(ib01ad, IB01AD)(
     const char *meth, const char *alg, const char *jobd,
     const char *batch, const char *conct, const char *ctrl,
     const int *nobr, const int *m, const int *l, const int *nsmp,
     double *u, const int *ldu, double *y, const int *ldy,
     int *n, double *r, const int *ldr, double *sv, double *rcond,
     const double *tol, int *iwork, double *dwork, const int *ldwork,
     int *iwarn, int *info,
     const int meth_len, const int alg_len, const int jobd_len,
     const int batch_len, const int conct_len, const int ctrl_len
 );
 
 /**
  * @brief Internal helper to calculate the minimum required workspace size (LDWORK).
  *
  * This function estimates the minimum LDWORK based on the parameters,
  * following the logic described in the SLICOT documentation for IB01AD.
  *
  * @param meth Algorithm method ('M' or 'N').
  * @param alg Algorithm type ('C', 'F', or 'Q').
  * @param jobd Singular value computation ('M' or 'N').
  * @param batch Batch processing ('F', 'I', 'L', 'O').
  * @param conct Concatenation ('C' or 'N').
  * @param m Number of inputs.
  * @param l Number of outputs.
  * @param nobr Number of block rows/columns in Hankel matrices.
  * @param nsmp Number of samples.
  * @param ldr_f Leading dimension of R as seen by Fortran.
  * @return Minimum required value for LDWORK.
  */
 static int calculate_min_ldwork_internal(char meth, char alg, char jobd, char batch, char conct,
                                          int m, int l, int nobr, int nsmp, int ldr_f) {
     int min_ldwork = 1; // Initialize to minimum possible value
     int ml = m + l;
     int ml_nobr = ml * nobr;
     int l_nobr = l * nobr;
     int m_nobr = m * nobr;
     int ns = nsmp - 2 * nobr + 1;
 
     // Convert to uppercase for comparison (already done in main wrapper)
 
     // Specific calculations based on algorithm and parameters (derived from documentation)
     if (alg == 'C') {
         if (batch == 'F' || batch == 'I') {
             min_ldwork = (conct == 'C') ? (4 * nobr - 2) * ml : 1;
         } else { // BATCH == 'L' or 'O'
             if (meth == 'M') {
                 if (jobd == 'M') {
                     min_ldwork = MAX(MAX((2*m-1)*nobr, ml_nobr), 5*l_nobr);
                      if (batch == 'L' && conct == 'C') {
                         min_ldwork = MAX(min_ldwork, (4 * nobr - 2) * ml);
                     }
                 } else { // JOBD == 'N'
                     min_ldwork = 5 * l_nobr;
                      if (batch == 'L' && conct == 'C') {
                         min_ldwork = MAX(min_ldwork, (4 * nobr - 2) * ml);
                     }
                 }
             } else { // METH == 'N'
                 min_ldwork = 5 * ml_nobr + 1;
             }
         }
     } else if (alg == 'F') {
         if (batch == 'O' || (batch == 'L' && conct == 'N')) {
             min_ldwork = ml * 4 * nobr * (ml + 1) + ml * 2 * nobr;
         } else if (batch != 'O' && conct == 'C') { // F/I/L and C
             min_ldwork = ml * 2 * nobr * (ml + 3);
         } else { // F/I and N
             min_ldwork = ml * 2 * nobr * (ml + 1);
         }
     } else if (alg == 'Q') {
         bool ldr_ge_ns = (ldr_f >= ns);
         if (batch == 'I' || batch == 'L') {
             if (conct == 'C') {
                 min_ldwork = 4 * (nobr + 1) * ml_nobr;
             } else { // CONCT == 'N'
                 min_ldwork = 6 * ml_nobr;
             }
         } else { // BATCH == 'F' or 'O'
             if (ldr_ge_ns) {
                 if (batch == 'F') {
                     min_ldwork = 4 * ml_nobr;
                 } else { // BATCH == 'O'
                     if (meth == 'M') {
                         min_ldwork = MAX(4 * ml_nobr, 5 * l_nobr);
                     } else { // METH == 'N'
                         min_ldwork = 5 * ml_nobr + 1;
                     }
                 }
             } else { // LDR < NS
                 min_ldwork = 6 * ml_nobr;
             }
         }
     }
 
     // Ensure at least 1 for workspace allocation
     return MAX(1, min_ldwork);
 }
 
 
 /* C wrapper function definition */
 SLICOT_C_WRAPPER_API
 int slicot_ib01ad(char meth, char alg, char jobd, char batch, char conct, char ctrl,
                   int nobr, int m, int l, int nsmp,
                   double *u, int ldu, double *y, int ldy,
                   int *n, double *r, int ldr, double *sv, double rcond,
                   double tol, int *iwarn, int row_major)
 {
     // 1. Variable declarations
     int info = 0;
     int iwarn_local = 0; // Local variable for Fortran output
     int n_local = 0;     // Local copy for Fortran call, initialized safely
 
     // Pointers for column-major copies
     double *u_cm = NULL;
     double *y_cm = NULL;
 
     // Workspace arrays allocated internally
     int *iwork = NULL;
     double *dwork = NULL;
     int liwork = 0;
     int ldwork = 0;
 
     // Fortran-compatible parameters
     char meth_upper = toupper(meth);
     char alg_upper = toupper(alg);
     char jobd_upper = toupper(jobd);
     char batch_upper = toupper(batch);
     char conct_upper = toupper(conct);
     char ctrl_upper = toupper(ctrl);
     const int meth_len = 1, alg_len = 1, jobd_len = 1;
     const int batch_len = 1, conct_len = 1, ctrl_len = 1;
 
     int ldu_f, ldy_f, ldr_f;
     double *u_ptr, *y_ptr; // Pointers to use in Fortran call
 
     // 2. Input parameter validation
 
     // Safely initialize n_local from input pointer *n
     if (n != NULL) {
         n_local = *n;
     } else {
         info = -15; // Error if n pointer is NULL
         goto cleanup; // Skip workspace allocation if basic params invalid
     }
 
     // Character parameter validation
     if (meth_upper != 'M' && meth_upper != 'N') { info = -1; goto cleanup; }
     if (alg_upper != 'C' && alg_upper != 'F' && alg_upper != 'Q') { info = -2; goto cleanup; }
     if (jobd_upper != 'M' && jobd_upper != 'N') { info = -3; goto cleanup; }
     if (batch_upper != 'F' && batch_upper != 'I' && batch_upper != 'L' && batch_upper != 'O') { info = -4; goto cleanup; }
     if (conct_upper != 'C' && conct_upper != 'N') { info = -5; goto cleanup; }
     if (ctrl_upper != 'C' && ctrl_upper != 'N') { info = -6; goto cleanup; }
 
     // Dimension validation
     if (nobr <= 0) { info = -7; goto cleanup; }
     if (m < 0) { info = -8; goto cleanup; }
     if (l <= 0) { info = -9; goto cleanup; }
     if (nsmp <= 0) { info = -10; goto cleanup; }
 
     // NSMP validation based on BATCH (from HTML docs)
     if (batch_upper == 'O' && nsmp < 2*(m+l+1)*nobr - 1) { info = -10; goto cleanup; }
     if ((batch_upper == 'F' || batch_upper == 'I') && nsmp < 2*nobr) { info = -10; goto cleanup; }
 
     // Calculate minimum required Fortran LDR (rows for Fortran)
     int min_ldr_f = (meth_upper == 'M' && jobd_upper == 'M') ? MAX(2*(m+l)*nobr, 3*m*nobr) : 2*(m+l)*nobr;
     min_ldr_f = MAX(1, min_ldr_f); // Ensure at least 1
 
     // Leading dimension checks (LDU, LDY, LDR)
     if (row_major) {
         if (m > 0 && ldu < m) { info = -12; goto cleanup; }
         if (l > 0 && ldy < l) { info = -14; goto cleanup; }
         if (ldr < 2 * (m + l) * nobr) { info = -17; goto cleanup; }
     } else {
         if (m > 0 && ldu < nsmp) { info = -12; goto cleanup; }
         if (l > 0 && ldy < nsmp) { info = -14; goto cleanup; }
         if (ldr < min_ldr_f) { info = -17; goto cleanup; }
     }
 
     // Set the leading dimension for R to be passed to Fortran.
     ldr_f = ldr;
 
     // Check required pointers (n already checked)
     if (m > 0 && u == NULL) { info = -11; goto cleanup; }
     if (l > 0 && y == NULL) { info = -13; goto cleanup; }
     if (r == NULL) { info = -16; goto cleanup; }
     if (sv == NULL) { info = -18; goto cleanup; }
     // iwarn is output only, can be NULL
 
     // 3. Internal Workspace Allocation
 
     // Calculate LIWORK based on documentation
     liwork = 3; // Minimum
     if (meth_upper == 'N') {
         liwork = MAX(3, (m + l) * nobr);
     } else if (meth_upper == 'M' && alg_upper == 'F') {
         liwork = MAX(3, m + l);
     }
     iwork = (int*)malloc((size_t)liwork * sizeof(int));
     CHECK_ALLOC(iwork); // Macro should handle error and jump to cleanup
 
     // Calculate LDWORK based on documentation
     ldwork = calculate_min_ldwork_internal(meth_upper, alg_upper, jobd_upper, batch_upper, conct_upper,
                                            m, l, nobr, nsmp, ldr_f);
     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Macro should handle error and jump to cleanup
 
 
     // 4. Memory allocation for column-major copies (if needed)
     size_t u_size = (size_t)nsmp * m; if (nsmp == 0 || m == 0) u_size = 0;
     size_t y_size = (size_t)nsmp * l; if (nsmp == 0 || l == 0) y_size = 0;
 
     if (row_major) {
         /* Allocate and convert U and Y from row-major to column-major */
         if (u_size > 0) {
             u_cm = (double*)malloc(u_size * sizeof(double));
             CHECK_ALLOC(u_cm);
             slicot_transpose_to_fortran(u, u_cm, nsmp, m, sizeof(double));
         }
 
         if (y_size > 0) {
             y_cm = (double*)malloc(y_size * sizeof(double));
             CHECK_ALLOC(y_cm);
             slicot_transpose_to_fortran(y, y_cm, nsmp, l, sizeof(double));
         }
 
         /* Set Fortran leading dimensions and pointers for U, Y */
         ldu_f = MAX(1, nsmp);
         ldy_f = MAX(1, nsmp);
         u_ptr = (u_size > 0) ? u_cm : NULL;
         y_ptr = (y_size > 0) ? y_cm : NULL;
     } else {
         /* Column-major case - use original arrays and leading dimensions */
         ldu_f = ldu;
         ldy_f = ldy;
         u_ptr = (u_size > 0) ? u : NULL;
         y_ptr = (y_size > 0) ? y : NULL;
     }
 
     // 5. Call Fortran Computation
     F77_FUNC(ib01ad, IB01AD)(
         &meth_upper, &alg_upper, &jobd_upper,
         &batch_upper, &conct_upper, &ctrl_upper,
         &nobr, &m, &l, &nsmp,
         u_ptr, &ldu_f, y_ptr, &ldy_f,
         &n_local, r, &ldr_f, sv, &rcond,
         &tol, iwork, dwork, &ldwork, // Use internally allocated workspaces
         &iwarn_local, &info,
         meth_len, alg_len, jobd_len,
         batch_len, conct_len, ctrl_len
     );
 
     // 6. Post-processing and Output
     if (info == 0) {
         // Update output parameters only on success
         if (n != NULL) *n = n_local; // Update output order
         if (iwarn != NULL) *iwarn = iwarn_local; // Update output warning
     }
 
 cleanup:
     // 7. Cleanup
     // Free internally allocated workspaces
     free(iwork);
     free(dwork);
     // Free temporary column-major copies if allocated
     free(u_cm);
     free(y_cm);
 
     // Check if info indicates a memory error from CHECK_ALLOC
     if (info == SLICOT_MEMORY_ERROR) {
        // Optionally print an error message or log
        // fprintf(stderr, "Error: Memory allocation failed in slicot_ib01ad.\n");
     }
 
     return info;
 }
