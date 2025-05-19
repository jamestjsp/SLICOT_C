/**
 * @file mb03vy.c
 * @brief C wrapper implementation for SLICOT routine MB03VY
 *
 * This file provides a C wrapper implementation for the SLICOT routine MB03VY,
 * which generates the orthogonal matrices Q_j from the elementary reflectors
 * computed by MB03VD.
 */

#include <stdlib.h>
#include <stddef.h> // For size_t
#include <string.h> // For memcpy

// Include the header file for this wrapper
#include "mb03vy.h"
// Include necessary SLICOT utility headers
#include "slicot_utils.h" // Assumed to contain MAX, MIN, CHECK_ALLOC, SLICOT_MEMORY_ERROR
#include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

/*
 * Declare the external Fortran routine using the F77_FUNC macro.
 * A is input/output.
 */
extern void F77_FUNC(mb03vy, MB03VY)(
    const int* n,           // INTEGER N
    const int* p,           // INTEGER P
    const int* ilo,         // INTEGER ILO
    const int* ihi,         // INTEGER IHI
    double* a,              // DOUBLE PRECISION A(LDA1,LDA2,P) (in/out)
    const int* lda1,        // INTEGER LDA1
    const int* lda2,        // INTEGER LDA2
    const double* tau,      // DOUBLE PRECISION TAU(LDTAU,P) (input)
    const int* ldtau,       // INTEGER LDTAU
    double* dwork,          // DOUBLE PRECISION DWORK(*)
    const int* ldwork,      // INTEGER LDWORK
    int* info               // INTEGER INFO (output)
);


/* C wrapper function definition */
SLICOT_EXPORT
int slicot_mb03vy(int n, int p, int ilo, int ihi,
                  double* a, int lda1, int lda2,
                  const double* tau, int ldtau,
                  int row_major)
{
    /* Local variables */
    int info = 0;
    double* dwork = NULL; // Workspace
    double* a_cm = NULL;  // Column-major copy if row_major=1
    double* tau_cm = NULL; // Column-major copy of tau if row_major=1
    int ldtau_cm = 0;     // Leading dimension for tau in column-major format
    int lda1_cm = 0;      // Leading dimension for Fortran a (rows)
    int lda2_cm = 0;      // Second dimension for Fortran a (cols)
    int ldwork = 0;       // Size of dwork array

    /* --- Input Parameter Validation --- */

    if (n < 0) { info = -1; goto cleanup; }
    if (p < 1) { info = -2; goto cleanup; }
    if (ilo < 1 || ilo > MAX(1, n)) { info = -3; goto cleanup; }
    if (ihi < MIN(ilo, n) || ihi > n) { info = -4; goto cleanup; }
    
    // For row-major, validate dimensions differently
    if (row_major) {
        if (a == NULL && n > 0) { info = -5; goto cleanup; }
        if (lda1 < n) { info = -6; goto cleanup; }   // Number of rows in 2D slice
        if (lda2 < n) { info = -7; goto cleanup; }   // Number of columns in 2D slice
    } else {
        // For column-major
        if (a == NULL && n > 0) { info = -5; goto cleanup; }
        if (lda1 < MAX(1, n)) { info = -6; goto cleanup; }   // Leading dimension (rows)
        if (lda2 < MAX(1, n)) { info = -7; goto cleanup; }   // Second dimension (cols)
    }
    
    if (tau == NULL && n > 1 && ihi > ilo) { info = -8; goto cleanup; }
    if (ldtau < MAX(1, n - 1)) { info = -9; goto cleanup; }

    if (info != 0) goto cleanup;

    /* --- Workspace Allocation --- */
    
    // Set up Fortran dimensions
    lda1_cm = MAX(1, n);  // Fortran leading dimension (rows)
    lda2_cm = MAX(1, n);  // Fortran second dimension (cols)
    ldtau_cm = MAX(1, n-1); // Fortran leading dimension for tau
    
    // Allocate DWORK (size N*N)
    ldwork = MAX(1, n*n);
    dwork = (double*)malloc((size_t)ldwork * sizeof(double));
    CHECK_ALLOC(dwork);

    /* --- Prepare Arrays for Fortran Routine --- */
    
    double* a_f_ptr = a;     // Pointer to pass to Fortran
    double* tau_f_ptr = (double*)tau; // Pointer to pass to Fortran (removing const)

    if (row_major && n > 0) {
        // Allocate temporary column-major arrays
        size_t a_cm_size = (size_t)lda1_cm * lda2_cm * p;
        a_cm = (double*)malloc(a_cm_size * sizeof(double));
        CHECK_ALLOC(a_cm);
        
        // Copy from row-major to column-major, slice by slice
        for (int k = 0; k < p; ++k) {
            // Calculate starting positions for this slice
            double* a_rm_slice = a + (size_t)k * lda1 * lda2;
            double* a_cm_slice = a_cm + (size_t)k * lda1_cm * lda2_cm;
            
            // Transpose each 2D slice
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    a_cm_slice[j * lda1_cm + i] = a_rm_slice[i * lda2 + j];
                }
            }
        }
        
        a_f_ptr = a_cm; // Pass the column-major copy to Fortran
        
        // For tau, allocate column-major buffer if needed
        if (n > 1 && ihi > ilo) {
            size_t tau_cm_size = (size_t)ldtau_cm * p;
            tau_cm = (double*)malloc(tau_cm_size * sizeof(double));
            CHECK_ALLOC(tau_cm);
            
            // Copy tau from row-major to column-major
            for (int j = 0; j < p; ++j) {
                for (int i = 0; i < n-1; ++i) {
                    tau_cm[j * ldtau_cm + i] = tau[j * ldtau + i];
                }
            }
            
            tau_f_ptr = tau_cm;
        }
    }

    /* Call the Fortran routine */
    F77_FUNC(mb03vy, MB03VY)(&n, &p, &ilo, &ihi,
                             a_f_ptr, &lda1_cm, &lda2_cm,
                             tau_f_ptr, &ldtau_cm,
                             dwork, &ldwork, &info);

    /* Copy results back for row-major case */
    if (row_major && n > 0 && info == 0) {
        // Copy from column-major to row-major, slice by slice
        for (int k = 0; k < p; ++k) {
            // Calculate starting positions for this slice
            double* a_rm_slice = a + (size_t)k * lda1 * lda2;
            double* a_cm_slice = a_cm + (size_t)k * lda1_cm * lda2_cm;
            
            // Transpose each 2D slice
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    a_rm_slice[i * lda2 + j] = a_cm_slice[j * lda1_cm + i];
                }
            }
        }
    }

cleanup:
    /* --- Cleanup --- */
    free(dwork);
    free(a_cm);
    free(tau_cm);

    /* Return the info code from the Fortran routine or memory error */
    return info;
}
