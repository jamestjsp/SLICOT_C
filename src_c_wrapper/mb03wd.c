/**
 * @file mb03wd.c
 * @brief C wrapper implementation for SLICOT routine MB03WD
 *
 * This file provides a C wrapper implementation for the SLICOT routine MB03WD,
 * which computes the Schur decomposition and eigenvalues of a product
 * of matrices H = H_1*...*H_p in periodic Hessenberg form.
 * NOTE: This wrapper assumes the 3D arrays 'h' and 'z' are passed in column-major order.
 */

#include <stdlib.h>
#include <ctype.h>
#include <stddef.h> // For size_t

// Include the header file for this wrapper
#include "mb03wd.h"
// Include necessary SLICOT utility headers
#include "slicot_utils.h" // Assumed to contain MAX, MIN, CHECK_ALLOC, SLICOT_MEMORY_ERROR
#include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

/*
 * Declare the external Fortran routine using the F77_FUNC macro.
 * H is input/output. Z is input/output depending on COMPZ.
 * WR, WI are output.
 */
extern void F77_FUNC(mb03wd, MB03WD)(
    const char* job,        // CHARACTER*1 JOB
    const char* compz,      // CHARACTER*1 COMPZ
    const int* n,           // INTEGER N
    const int* p,           // INTEGER P
    const int* ilo,         // INTEGER ILO
    const int* ihi,         // INTEGER IHI
    const int* iloz,        // INTEGER ILOZ
    const int* ihiz,        // INTEGER IHIZ
    double* h,              // DOUBLE PRECISION H(LDH1,LDH2,P) (in/out)
    const int* ldh1,        // INTEGER LDH1
    const int* ldh2,        // INTEGER LDH2
    double* z,              // DOUBLE PRECISION Z(LDZ1,LDZ2,P) (in/out)
    const int* ldz1,        // INTEGER LDZ1
    const int* ldz2,        // INTEGER LDZ2
    double* wr,             // DOUBLE PRECISION WR(*) (output)
    double* wi,             // DOUBLE PRECISION WI(*) (output)
    double* dwork,          // DOUBLE PRECISION DWORK(*)
    const int* ldwork,      // INTEGER LDWORK
    int* info,              // INTEGER INFO (output)
    int job_len,            // Hidden length
    int compz_len           // Hidden length
);


/* C wrapper function definition */
SLICOT_C_WRAPPER_API
int slicot_mb03wd(char job, char compz, int n, int p, int ilo, int ihi,
                  int iloz, int ihiz,
                  double* h, int ldh1, int ldh2,
                  double* z, int ldz1, int ldz2,
                  double* wr, double* wi,
                  int row_major) // row_major ignored for 'h', 'z'
{
    /* Local variables */
    int info = 0;
    double* dwork = NULL; // Workspace
    int ldwork = 0;
    const int job_len = 1, compz_len = 1;

    char job_upper = toupper(job);
    char compz_upper = toupper(compz);

    /* --- Input Parameter Validation --- */

    if (n < 0) { info = -3; goto cleanup; }
    if (p < 1) { info = -4; goto cleanup; }
    if (ilo < 1 || ilo > MAX(1, n)) { info = -5; goto cleanup; }
    if (ihi < MIN(ilo, n) || ihi > n) { info = -6; goto cleanup; }
    if (iloz < 1 || iloz > ilo) { info = -7; goto cleanup; }
    if (ihiz < ihi || ihiz > n) { info = -8; goto cleanup; }
    if (job_upper != 'E' && job_upper != 'S') { info = -1; goto cleanup; }
    if (compz_upper != 'N' && compz_upper != 'I' && compz_upper != 'V') { info = -2; goto cleanup; }

    // Check leading dimensions
    if (ldh1 < MAX(1, n)) { info = -10; goto cleanup; }
    if (ldh2 < MAX(1, n)) { info = -11; goto cleanup; }
    if (compz_upper == 'N') {
        if (ldz1 < 1) { info = -13; goto cleanup; }
        if (ldz2 < 1) { info = -14; goto cleanup; }
    } else {
        if (ldz1 < MAX(1, n)) { info = -13; goto cleanup; }
        if (ldz2 < MAX(1, n)) { info = -14; goto cleanup; }
    }


    /* --- Workspace Allocation --- */

    // Allocate DWORK (size IHI-ILO+P-1) - No query needed
    ldwork = MAX(1, ihi - ilo + p - 1); // Ensure minimum size 1
    dwork = (double*)malloc((size_t)ldwork * sizeof(double));
    CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure

    /* --- Prepare Arrays and Call Fortran Routine --- */

    // NOTE: Assuming 'h' and 'z' are already in column-major (Fortran) layout.
    // No transposition is performed for these 3D arrays.
    // The 'row_major' flag is effectively ignored for 'h' and 'z'.

    /* Call the Fortran routine directly */
    F77_FUNC(mb03wd, MB03WD)(&job_upper, &compz_upper, &n, &p, &ilo, &ihi,
                             &iloz, &ihiz,
                             h, &ldh1, &ldh2,
                             (compz_upper == 'N' ? NULL : z), &ldz1, &ldz2,
                             wr, wi,
                             dwork, &ldwork, &info,
                             job_len, compz_len);
    // H, Z (if COMPZ != 'N'), WR, WI are modified in place.

cleanup:
    /* --- Cleanup --- */
    free(dwork);

    /* Return the info code from the Fortran routine or SLICOT_MEMORY_ERROR */
    return info;
}
