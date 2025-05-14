/**
 * @file mb03wd.c
 * @brief C wrapper implementation for SLICOT routine MB03WD
 *
 * This file provides a C wrapper implementation for the SLICOT routine MB03WD,
 * which computes the Schur decomposition and eigenvalues of a product
 * of matrices H = H_1*...*H_p in periodic Hessenberg form.
 * NOTE: This wrapper assumes the 3D arrays 'h' and 'z' are passed in column-major order
 * and the 'row_major' flag is IGNORED for 'h' and 'z'.
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
SLICOT_EXPORT
int slicot_mb03wd(char job, char compz, int n, int p, int ilo, int ihi,
                  int iloz, int ihiz,
                  double* h, int ldh1, int ldh2,
                  double* z, int ldz1, int ldz2,
                  double* wr, double* wi,
                  int row_major) // row_major ignored for 'h', 'z'
{
    /* Local variables */
    int info = 0;
    int ldwork_val = 0; // Renamed to avoid conflict, size is calculated directly
    double* dwork = NULL; // Workspace
    // No iwork needed for this routine.
    const int job_len = 1, compz_len = 1;

    char job_upper = toupper(job);
    char compz_upper = toupper(compz);

    // No _cm pointers needed as 'h', 'z' are assumed column-major and row_major is ignored.

    /* --- Input Parameter Validation --- */
    if (job_upper != 'E' && job_upper != 'S') { info = -1; goto cleanup; } // JOB
    if (compz_upper != 'N' && compz_upper != 'I' && compz_upper != 'V') { info = -2; goto cleanup; } // COMPZ
    if (n < 0) { info = -3; goto cleanup; } // N
    if (p < 1) { info = -4; goto cleanup; } // P

    // ILO (arg 5)
    if (n == 0) { // Special case for N=0
        if (ilo != 1) { info = -5; goto cleanup; }
    } else { // N > 0
        if (ilo < 1 || ilo > n) { info = -5; goto cleanup; }
    }
    // IHI (arg 6)
    if (n == 0) { // Special case for N=0
         if (ihi != 0) { info = -6; goto cleanup; }
    } else { // N > 0
        if (ihi < MIN(ilo, n) || ihi > n) { info = -6; goto cleanup; }
    }
    // ILOZ (arg 7)
     if (n == 0) { // If N=0, ILO=1. ILOZ must be <= ILO.
        if (iloz != 1) { info = -7; goto cleanup; }
    } else { // N > 0
        if (iloz < 1 || iloz > ilo) { info = -7; goto cleanup; }
    }
    // IHIZ (arg 8)
    if (n == 0) { // If N=0, IHI=0. IHIZ must be >= IHI.
        if (ihiz != 0) { info = -8; goto cleanup; }
    } else { // N > 0
         if (ihiz < ihi || ihiz > n) { info = -8; goto cleanup; }
    }
    
    // H (arg 9) - NULL check
    if (n > 0 && h == NULL) { info = -9; goto cleanup; }

    // LDH1 (arg 10)
    if (ldh1 < MAX(1, n)) { info = -10; goto cleanup; }
    // LDH2 (arg 11)
    if (ldh2 < MAX(1, n)) { info = -11; goto cleanup; }

    // Z (arg 12) - NULL check
    if (compz_upper != 'N' && n > 0 && z == NULL) { info = -12; goto cleanup; }

    // LDZ1 (arg 13)
    if (compz_upper == 'N') {
        if (ldz1 < 1) { info = -13; goto cleanup; }
    } else { // COMPZ = 'I' or 'V'
        if (ldz1 < MAX(1, n)) { info = -13; goto cleanup; }
    }
    // LDZ2 (arg 14)
    if (compz_upper == 'N') {
        if (ldz2 < 1) { info = -14; goto cleanup; }
    } else { // COMPZ = 'I' or 'V'
        if (ldz2 < MAX(1, n)) { info = -14; goto cleanup; }
    }

    // WR (arg 15) - NULL check
    if (n > 0 && wr == NULL) { info = -15; goto cleanup; }
    // WI (arg 16) - NULL check
    if (n > 0 && wi == NULL) { info = -16; goto cleanup; }
    
    // The 'row_major' flag is intentionally ignored for validation here due to the 3D arrays.

    /* --- Workspace Allocation --- */

    // Allocate DWORK (size IHI-ILO+P-1) - No query needed for MB03WD
    // Ensure ihi >= ilo for subtraction, which is guaranteed by prior validation if N > 0.
    // If N=0, ILO=1, IHI=0. So IHI-ILO = -1.
    // LDWORK >= IHI-ILO+P-1.
    if (n == 0) { // If N=0, ILO=1, IHI=0. IHI-ILO+P-1 = 0-1+P-1 = P-2
        ldwork_val = MAX(1, p - 2);
    } else { // N > 0
        ldwork_val = MAX(1, ihi - ilo + p - 1);
    }
    
    if (ldwork_val > 0) {
        dwork = (double*)malloc((size_t)ldwork_val * sizeof(double));
        CHECK_ALLOC(dwork); // Sets info to SLICOT_MEMORY_ERROR and jumps to cleanup on failure
    } else { // This case should ideally not be hit if P>=1 due to MAX(1, ...)
        dwork = NULL; // Or handle as an error if ldwork_val becomes non-positive unexpectedly
    }


    /* --- Call the computational routine --- */
    // 'h' and 'z' are used directly as they are assumed column-major.
    // If N=0, relevant pointers (h, z, wr, wi) might be NULL but should not be dereferenced by Fortran.
    // If N>0, they've been checked to be non-NULL.
    double* z_for_fortran = z;
    if (compz_upper == 'N') {
        z_for_fortran = NULL; // Pass actual NULL if Z is not used.
                              // LDZ1, LDZ2 still need to be >= 1 for Fortran interface.
    }
    
    // For N=0 case, if h, wr, wi were NULL, that's fine.
    // If N>0, they are not NULL.
    double* h_for_fortran = (n > 0) ? h : NULL; 
    double* wr_for_fortran = (n > 0) ? wr : NULL;
    double* wi_for_fortran = (n > 0) ? wi : NULL;


    F77_FUNC(mb03wd, MB03WD)(&job_upper, &compz_upper, &n, &p, &ilo, &ihi,
                             &iloz, &ihiz,
                             h_for_fortran, &ldh1, &ldh2,
                             z_for_fortran, &ldz1, &ldz2,
                             wr_for_fortran, wi_for_fortran,
                             dwork, &ldwork_val, &info,
                             job_len, compz_len);

cleanup:
    /* --- Cleanup --- */
    free(dwork);

    return info;
}
