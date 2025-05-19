/**
 * @file mb03wd.c
 * @brief C wrapper implementation for SLICOT routine MB03WD
 *
 * This file provides a C wrapper implementation for the SLICOT routine MB03WD,
 * which computes the Schur decomposition and eigenvalues of a product
 * of matrices H = H_1*...*H_p in periodic Hessenberg form.
 */

#include <stdlib.h>
#include <stddef.h> // For size_t
#include <string.h> // For memcpy
#include <ctype.h> // For toupper

// Include the header file for this wrapper
#include "mb03wd.h"
// Include necessary SLICOT utility headers
#include "slicot_utils.h" // For MAX, MIN, CHECK_ALLOC, SLICOT_MEMORY_ERROR
#include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

/*
 * Declare the external Fortran routine using the F77_FUNC macro.
 * H and Z are input/output.
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
    double* h,              // DOUBLE PRECISION H(LDH1,LDH2,P)
    const int* ldh1,        // INTEGER LDH1
    const int* ldh2,        // INTEGER LDH2
    double* z,              // DOUBLE PRECISION Z(LDZ1,LDZ2,P)
    const int* ldz1,        // INTEGER LDZ1
    const int* ldz2,        // INTEGER LDZ2
    double* wr,             // DOUBLE PRECISION WR(N)
    double* wi,             // DOUBLE PRECISION WI(N)
    double* dwork,          // DOUBLE PRECISION DWORK(*)
    const int* ldwork,      // INTEGER LDWORK
    int* info,              // INTEGER INFO
    int job_len,            // JOB length (added by f2c/C wrapper)
    int compz_len           // COMPZ length (added by f2c/C wrapper)
);

/* C wrapper function definition */
SLICOT_EXPORT
int slicot_mb03wd(char job, char compz, int n, int p, int ilo, int ihi,
                  int iloz, int ihiz,
                  double* h, int ldh1, int ldh2,
                  double* z, int ldz1, int ldz2,
                  double* wr, double* wi,
                  int row_major)
{
    /* Local variables */
    int info = 0;
    char job_f = toupper(job);
    char compz_f = toupper(compz);
    double* dwork = NULL; // Workspace
    double* h_cm = NULL;  // Column-major copy if row_major=1
    double* z_cm = NULL;  // Column-major copy if row_major=1
    int ldh1_cm = 0;      // Leading dimension for Fortran h (rows)
    int ldh2_cm = 0;      // Second dimension for Fortran h (cols)
    int ldz1_cm = 0;      // Leading dimension for Fortran z (rows)
    int ldz2_cm = 0;      // Second dimension for Fortran z (cols)
    int ldwork = 0;       // Size of dwork array
    int need_z = 0;       // Flag to indicate if Z is needed

    /* --- Input Parameter Validation --- */
    if (job_f != 'E' && job_f != 'S') { info = -1; goto cleanup; }
    if (compz_f != 'N' && compz_f != 'I' && compz_f != 'V') { info = -2; goto cleanup; }
    if (n < 0) { info = -3; goto cleanup; }
    if (p < 1) { info = -4; goto cleanup; }
    if (ilo < 1 || ilo > MAX(1, n)) { info = -5; goto cleanup; }
    if (ihi < MIN(ilo, n) || ihi > n) { info = -6; goto cleanup; }
    if (iloz < 1 || iloz > ilo) { info = -7; goto cleanup; }
    if (ihiz < ihi || ihiz > n) { info = -8; goto cleanup; }
    
    need_z = (compz_f == 'I' || compz_f == 'V');
    
    // For row-major, validate dimensions differently
    if (row_major) {
        if (h == NULL && n > 0) { info = -9; goto cleanup; }
        if (ldh1 < n) { info = -10; goto cleanup; }  // Number of rows in 2D slice
        if (ldh2 < n) { info = -11; goto cleanup; }  // Number of columns in 2D slice
        
        if (need_z) {
            if (z == NULL && n > 0) { info = -12; goto cleanup; }
            if (ldz1 < n) { info = -13; goto cleanup; }  // Number of rows in 2D slice
            if (ldz2 < n) { info = -14; goto cleanup; }  // Number of columns in 2D slice
        } else {
            if (ldz1 < 1) { info = -13; goto cleanup; }  // Minimal dimension for unused Z
            if (ldz2 < 1) { info = -14; goto cleanup; }  // Minimal dimension for unused Z
        }
    } else {
        // For column-major
        if (h == NULL && n > 0) { info = -9; goto cleanup; }
        if (ldh1 < MAX(1, n)) { info = -10; goto cleanup; }  // Leading dimension (rows)
        if (ldh2 < MAX(1, n)) { info = -11; goto cleanup; }  // Second dimension (cols)
        
        if (need_z) {
            if (z == NULL && n > 0) { info = -12; goto cleanup; }
            if (ldz1 < MAX(1, n)) { info = -13; goto cleanup; }  // Leading dimension (rows)
            if (ldz2 < MAX(1, n)) { info = -14; goto cleanup; }  // Second dimension (cols)
        } else {
            if (ldz1 < 1) { info = -13; goto cleanup; }  // Minimal dimension for unused Z
            if (ldz2 < 1) { info = -14; goto cleanup; }  // Minimal dimension for unused Z
        }
    }
    
    if (wr == NULL && n > 0) { info = -15; goto cleanup; }
    if (wi == NULL && n > 0) { info = -16; goto cleanup; }

    if (info != 0) goto cleanup;

    /* --- Workspace Allocation --- */
    
    // Set up Fortran dimensions
    ldh1_cm = MAX(1, n);  // Fortran leading dimension (rows)
    ldh2_cm = MAX(1, n);  // Fortran second dimension (cols)
    ldz1_cm = need_z ? MAX(1, n) : 1;  // Fortran leading dimension (rows)
    ldz2_cm = need_z ? MAX(1, n) : 1;  // Fortran second dimension (cols)
    
    // Calculate workspace size required by MB03WD
    ldwork = MAX(1, 2 * p * n * n + n * n + 8 * n);
    dwork = (double*)malloc((size_t)ldwork * sizeof(double));
    CHECK_ALLOC(dwork);

    /* --- Prepare Arrays for Fortran Routine --- */
    
    double* h_f_ptr = h;     // Pointer to pass to Fortran
    double* z_f_ptr = z;     // Pointer to pass to Fortran

    if (row_major && n > 0) {
        // Allocate temporary column-major arrays for H
        size_t h_cm_size = (size_t)ldh1_cm * ldh2_cm * p;
        h_cm = (double*)malloc(h_cm_size * sizeof(double));
        CHECK_ALLOC(h_cm);
        
        // Copy H from row-major to column-major, slice by slice
        for (int k = 0; k < p; ++k) {
            // Calculate starting positions for this slice
            double* h_rm_slice = h + (size_t)k * ldh1 * ldh2;
            double* h_cm_slice = h_cm + (size_t)k * ldh1_cm * ldh2_cm;
            
            // Transpose each 2D slice
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    h_cm_slice[j * ldh1_cm + i] = h_rm_slice[i * ldh2 + j];
                }
            }
        }
        
        h_f_ptr = h_cm; // Pass the column-major copy to Fortran
        
        // Handle Z similarly if needed
        if (need_z) {
            size_t z_cm_size = (size_t)ldz1_cm * ldz2_cm * p;
            z_cm = (double*)malloc(z_cm_size * sizeof(double));
            CHECK_ALLOC(z_cm);
            
            // If compz='I', Fortran will initialize Z to identity
            // If compz='V', we need to copy user's Z to column-major
            if (compz_f == 'V' && z != NULL) {
                for (int k = 0; k < p; ++k) {
                    // Calculate starting positions for this slice
                    double* z_rm_slice = z + (size_t)k * ldz1 * ldz2;
                    double* z_cm_slice = z_cm + (size_t)k * ldz1_cm * ldz2_cm;
                    
                    // Transpose each 2D slice
                    for (int i = 0; i < n; ++i) {
                        for (int j = 0; j < n; ++j) {
                            z_cm_slice[j * ldz1_cm + i] = z_rm_slice[i * ldz2 + j];
                        }
                    }
                }
            }
            
            z_f_ptr = z_cm; // Pass the column-major copy to Fortran
        } else {
            z_f_ptr = NULL; // Not needed if compz='N'
        }
    }

    /* Call the Fortran routine */
    F77_FUNC(mb03wd, MB03WD)(
        &job_f, &compz_f, &n, &p, &ilo, &ihi, &iloz, &ihiz,
        h_f_ptr, &ldh1_cm, &ldh2_cm,
        z_f_ptr, &ldz1_cm, &ldz2_cm,
        wr, wi, dwork, &ldwork, &info, 1, 1);

    /* Copy results back for row-major case */
    if (row_major && n > 0 && info == 0) {
        // Copy H from column-major back to row-major if job='S'
        if (job_f == 'S') {
            for (int k = 0; k < p; ++k) {
                // Calculate starting positions for this slice
                double* h_rm_slice = h + (size_t)k * ldh1 * ldh2;
                double* h_cm_slice = h_cm + (size_t)k * ldh1_cm * ldh2_cm;
                
                // Transpose each 2D slice
                for (int i = 0; i < n; ++i) {
                    for (int j = 0; j < n; ++j) {
                        h_rm_slice[i * ldh2 + j] = h_cm_slice[j * ldh1_cm + i];
                    }
                }
            }
        }
        
        // Copy Z from column-major back to row-major if needed
        if (need_z && z != NULL && z_cm != NULL) {
            for (int k = 0; k < p; ++k) {
                // Calculate starting positions for this slice
                double* z_rm_slice = z + (size_t)k * ldz1 * ldz2;
                double* z_cm_slice = z_cm + (size_t)k * ldz1_cm * ldz2_cm;
                
                // Transpose each 2D slice
                for (int i = 0; i < n; ++i) {
                    for (int j = 0; j < n; ++j) {
                        z_rm_slice[i * ldz2 + j] = z_cm_slice[j * ldz1_cm + i];
                    }
                }
            }
        }
    }

cleanup:
    /* --- Cleanup --- */
    free(dwork);
    free(h_cm);
    free(z_cm);

    /* Return the info code from the Fortran routine or memory error */
    return info;
}
