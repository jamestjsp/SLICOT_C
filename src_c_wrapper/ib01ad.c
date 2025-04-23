/**
 * @file ib01ad.c
 * @brief C wrapper implementation for SLICOT routine IB01AD
 *
 * This file provides a C wrapper implementation for the SLICOT routine IB01AD,
 * which estimates the system order and computes the triangular factor 
 * of the concatenated block-Hankel matrices for subsequent 
 * system identification in the SLICOT IB01 family of routines.
 */

#include <stdlib.h>
#include <ctype.h>
#include <stddef.h>
#include <string.h> // For memset, memcpy
#include <stdbool.h> // For bool type

#include "ib01ad.h"
#include "slicot_utils.h"
#include "slicot_f77.h"


/*
 * Declare the external Fortran routine using the F77_FUNC macro.
 * This handles potential name mangling issues between C and Fortran compilers.
 * All arguments are passed by reference (pointers).
 * The hidden string lengths for char arguments are explicitly included at the end.
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

// Helper function to calculate minimum LDWORK based on documentation
static int calculate_min_ldwork(char meth, char alg, char jobd, char batch, char conct,
                                int m, int l, int nobr, int nsmp, int ldr_f) {
    int min_ldwork = 1;
    int ml = m + l;
    int ml_nobr = ml * nobr;
    int l_nobr = l * nobr;
    int m_nobr = m * nobr;
    int ns = nsmp - 2 * nobr + 1;
    
    // Default to a reasonable minimum size for any algorithm
    min_ldwork = MAX(min_ldwork, 6 * ml * nobr);

    // Specific calculations based on algorithm and parameters
    if (alg == 'C') {
        if (batch == 'F' || batch == 'I') {
            min_ldwork = (conct == 'C') ? (4 * nobr - 2) * ml : 1;
        } else { // BATCH == 'L' or 'O'
            if (meth == 'M') {
                if (jobd == 'M') {
                    min_ldwork = MAX(MAX((2 * m - 1) * nobr, ml_nobr), 5 * l_nobr);
                } else { // JOBD == 'N'
                    min_ldwork = 5 * l_nobr;
                }
                if (batch == 'L' && conct == 'C') {
                    min_ldwork = MAX(min_ldwork, (4 * nobr - 2) * ml);
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

    // Use a larger workspace for safety - look at Fortran example's workspace size
    min_ldwork = MAX(min_ldwork, 6 * ml * nobr);
    
    // For workspace query, ensure at least one element
    return MAX(1, min_ldwork);
}

/* C wrapper function definition */
SLICOT_C_WRAPPER_API
int slicot_ib01ad(char meth, char alg, char jobd, char batch, char conct, char ctrl, 
                 int nobr, int m, int l, int nsmp, 
                 double *u, int ldu, double *y, int ldy,
                 int *n, double *r, int ldr, double *sv, double rcond, 
                 double tol, int *iwork, double *dwork, int ldwork, 
                 int *iwarn, int row_major)
{
    // 1. Variable declarations
    int info = 0;
    int iwarn_local = 0; // Local variable for Fortran output
    int n_local = (n != NULL) ? *n : 0; // Local copy for Fortran call

    // Pointers for column-major copies
    double *u_cm = NULL;
    double *y_cm = NULL;
    
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
    if (toupper(meth) != 'M' && toupper(meth) != 'N') { info = -1; goto cleanup; }
    if (toupper(alg) != 'C' && toupper(alg) != 'F' && toupper(alg) != 'Q') { info = -2; goto cleanup; }
    if (toupper(jobd) != 'M' && toupper(jobd) != 'N') { info = -3; goto cleanup; }
    if (toupper(batch) != 'F' && toupper(batch) != 'I' && toupper(batch) != 'L' && toupper(batch) != 'O') { info = -4; goto cleanup; }
    if (toupper(conct) != 'C' && toupper(conct) != 'N') { info = -5; goto cleanup; }
    if (toupper(ctrl) != 'C' && toupper(ctrl) != 'N') { info = -6; goto cleanup; }
    if (nobr <= 0) { info = -7; goto cleanup; }
    if (m < 0) { info = -8; goto cleanup; }
    if (l <= 0) { info = -9; goto cleanup; }
    if (nsmp <= 0) { info = -10; goto cleanup; }

    // NSMP validation based on BATCH (from HTML docs)
    if (batch_upper == 'O' && nsmp < 2*(m+l+1)*nobr - 1) { info = -10; goto cleanup; }
    if ((batch_upper == 'F' || batch_upper == 'I') && nsmp < 2*nobr) { info = -10; goto cleanup; }
    // Note: Total samples check for BATCH='L' is harder to enforce here.

    // Calculate minimum required Fortran LDR
    int min_ldr_f = (meth_upper == 'M' && jobd_upper == 'M') ? MAX(2*(m+l)*nobr, 3*m*nobr) : 2*(m+l)*nobr;
    min_ldr_f = MAX(1, min_ldr_f); // Ensure at least 1

    // Leading dimension checks
    if (row_major) {
        if (m > 0 && ldu < m) { info = -12; goto cleanup; }
        if (l > 0 && ldy < l) { info = -14; goto cleanup; }
        // For row-major, ldr is the number of columns
        if (ldr < 2*(m+l)*nobr) { info = -17; goto cleanup; } 
        ldr_f = min_ldr_f; // Fortran needs the minimum required row count
    } else {
        // For column-major, ldu/ldy/ldr are number of rows
        if (m > 0 && ldu < nsmp) { info = -12; goto cleanup; }
        if (l > 0 && ldy < nsmp) { info = -14; goto cleanup; }
        if (ldr < min_ldr_f) { info = -17; goto cleanup; }
        ldr_f = ldr; // Fortran uses the provided row count
    }

    // Check required pointers
    if (m > 0 && u == NULL) { info = -11; goto cleanup; }
    if (l > 0 && y == NULL) { info = -13; goto cleanup; }
    if (n == NULL) { info = -15; goto cleanup; }
    if (r == NULL) { info = -16; goto cleanup; }
    if (sv == NULL) { info = -18; goto cleanup; }
    // iwork, dwork can be NULL if size is 0
    
    // 3. Memory allocation for column-major copies
    size_t u_size = (size_t)nsmp * m; if (nsmp == 0 || m == 0) u_size = 0;
    size_t y_size = (size_t)nsmp * l; if (nsmp == 0 || l == 0) y_size = 0;
        
    if (row_major) {
        /* Allocate and convert from row-major to column-major */
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
        
        /* Set Fortran leading dimensions and pointers */
        ldu_f = MAX(1, nsmp); // Column-major: rows are the leading dimension
        ldy_f = MAX(1, nsmp);
        u_ptr = (u_size > 0) ? u_cm : NULL; // Use NULL if size is 0
        y_ptr = (y_size > 0) ? y_cm : NULL;
    } else {
        /* Column-major case - use original arrays */
        ldu_f = ldu;
        ldy_f = ldy;
        u_ptr = (u_size > 0) ? u : NULL; // Use NULL if size is 0
        y_ptr = (y_size > 0) ? y : NULL;
    }
      
    // 4. Workspace Validation and Execution
    // Calculate minimum required workspace size
    int min_ldwork = calculate_min_ldwork(meth_upper, alg_upper, jobd_upper, batch_upper, conct_upper,
                                          m, l, nobr, nsmp, ldr_f);
    
    // For workspace query, we need at least 1 element
    if (ldwork == -1) {
        // This is a workspace query, need at least 1 element
        if (dwork == NULL) { info = -22; goto cleanup; }
        
        // Set min_ldwork to a larger value for the query
        min_ldwork = MAX(min_ldwork, 5000);  // Use a reasonable minimum for query
    }
    else if (ldwork < min_ldwork) {
        // Actual computation with insufficient workspace
        info = -23;
        // Optionally store minimum required size in dwork[0] if dwork is not NULL
        if (dwork != NULL && ldwork >= 1) {
            dwork[0] = (double)min_ldwork;
        }
        goto cleanup;
    }
    else if (min_ldwork > 0 && dwork == NULL) {
        // No workspace provided but needed
        info = -22;
        goto cleanup;
    }

    // Basic check for iwork based on documentation
    int min_liwork = 3;
    if (meth_upper == 'N') min_liwork = MAX(3, (m+l)*nobr);
    else if (meth_upper == 'M' && alg_upper == 'F') min_liwork = MAX(3, m+l);
    if (min_liwork > 0 && iwork == NULL) { info = -21; goto cleanup; }

    // Call the actual computation
    F77_FUNC(ib01ad, IB01AD)(
        &meth_upper, &alg_upper, &jobd_upper,
        &batch_upper, &conct_upper, &ctrl_upper,
        &nobr, &m, &l, &nsmp,
        u_ptr, &ldu_f, y_ptr, &ldy_f,
        &n_local, r, &ldr_f, sv, &rcond,
        &tol, iwork, dwork, &ldwork,
        &iwarn_local, &info,
        meth_len, alg_len, jobd_len,
        batch_len, conct_len, ctrl_len
    );
    
    // For workspace query, if successful or info = -23, return required size
    if (ldwork == -1 && (info == 0 || info == -23)) {
        // If we can read the required size from dwork
        if (dwork[0] > 0) {
            // Return success (info = 0) since we got the required size
            info = 0;
        }
    }
    
    // 5. Post-processing and Output
    if (info == 0) {
        if (n != NULL) *n = n_local; // Update output order
        if (iwarn != NULL) *iwarn = iwarn_local; // Update output warning
    }
    
cleanup:
    // 6. Cleanup
    free(u_cm);
    free(y_cm);
    return info;
}
