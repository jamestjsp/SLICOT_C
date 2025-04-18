/**
 * @file mb03vd.c
 * @brief C wrapper implementation for SLICOT routine MB03VD
 *
 * This file provides a C wrapper implementation for the SLICOT routine MB03VD,
 * which reduces a product of p real general matrices A = A_1*...*A_p
 * to periodic Hessenberg form H = H_1*...*H_p using orthogonal
 * similarity transformations.
 * NOTE: This wrapper assumes the 3D array 'a' is passed in column-major order.
 */

#include <stdlib.h>
#include <stddef.h> // For size_t

// Include the header file for this wrapper
#include "mb03vd.h"
// Include necessary SLICOT utility headers
#include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR
#include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

/*
 * Declare the external Fortran routine using the F77_FUNC macro.
 * A and TAU are input/output.
 */
extern void F77_FUNC(mb03vd, MB03VD)(
    const int* n,           // INTEGER N
    const int* p,           // INTEGER P
    const int* ilo,         // INTEGER ILO
    const int* ihi,         // INTEGER IHI
    double* a,              // DOUBLE PRECISION A(LDA1,LDA2,P) (in/out)
    const int* lda1,        // INTEGER LDA1
    const int* lda2,        // INTEGER LDA2
    double* tau,            // DOUBLE PRECISION TAU(LDTAU,P) (output)
    const int* ldtau,       // INTEGER LDTAU
    double* dwork,          // DOUBLE PRECISION DWORK(*)
    int* info               // INTEGER INFO (output)
);


/* C wrapper function definition */
int slicot_mb03vd(int n, int p, int ilo, int ihi,
                  double* a, int lda1, int lda2,
                  double* tau, int ldtau,
                  int row_major) // row_major ignored for 'a'
{
    /* Local variables */
    int info = 0;
    double* dwork = NULL; // Workspace
    int dwork_size = 0;

    /* --- Input Parameter Validation --- */

    if (n < 0) { info = -1; goto cleanup; }
    if (p < 1) { info = -2; goto cleanup; }
    if (ilo < 1 || ilo > MAX(1, n)) { info = -3; goto cleanup; }
    if (ihi < MIN(ilo, n) || ihi > n) { info = -4; goto cleanup; }
    if (lda1 < MAX(1, n)) { info = -6; goto cleanup; }
    if (lda2 < MAX(1, n)) { info = -7; goto cleanup; }
    if (ldtau < MAX(1, n - 1)) { info = -9; goto cleanup; }

    /* --- Workspace Allocation --- */

    // Allocate DWORK (size N) - No query needed
    dwork_size = MAX(1, n); // Ensure minimum size 1
    dwork = (double*)malloc((size_t)dwork_size * sizeof(double));
    CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure

    /* --- Prepare Arrays and Call Fortran Routine --- */

    // NOTE: Assuming 'a' is already in column-major (Fortran) layout.
    // No transposition is performed for the 3D array 'a'.
    // The 'row_major' flag is effectively ignored for 'a'.

    /* Call the Fortran routine directly */
    F77_FUNC(mb03vd, MB03VD)(&n, &p, &ilo, &ihi,
                             a, &lda1, &lda2,
                             tau, &ldtau,
                             dwork, &info);
    // A and TAU are modified in place.

cleanup:
    /* --- Cleanup --- */
    free(dwork);

    /* Return the info code from the Fortran routine or SLICOT_MEMORY_ERROR */
    return info;
}
