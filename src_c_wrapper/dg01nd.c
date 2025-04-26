/**  
 * @file dg01nd.c  
 * @brief C wrapper for SLICOT routine DG01ND.  
 * @details Computes the discrete Fourier transform, or inverse Fourier transform, of a real signal.
 * Workspace is allocated internally by this wrapper.
 * Input/output format is handled via the row_major parameter.  
 */

#include <stdlib.h> // For malloc, free  
#include <ctype.h>  // For toupper  
#include <stddef.h> // For size_t  
#include <stdio.h>  // For error logging

#include "dg01nd.h" // Public header for this wrapper  
#include "slicot_utils.h"  // Provides CHECK_ALLOC, SLICOT_MEMORY_ERROR, etc.  
#include "slicot_f77.h"    // Provides F77_FUNC macro for Fortran name mangling

/* External Fortran routine declaration */  
extern void F77_FUNC(dg01nd, DG01ND)(  
    const char* indi,
    const int* n,
    double* xr,
    double* xi,
    int* info,
    int indi_len  // Hidden string length
);

/* C wrapper function definition */  
SLICOT_EXPORT
int slicot_dg01nd(char indi, int n, double* xr, double* xi, int row_major) {  
    // 1. Variable declarations  
    int info = 0;
    
    // Pointers for column-major copies (if row_major is used)  
    double* xr_cm = NULL;  
    double* xi_cm = NULL;

    // 2. Input parameter validation (Check BEFORE allocating memory)  
    if (n < 2) { info = -2; goto cleanup; } // Check minimum value for n

    // Check if n is a power of 2
    int n_check = n;
    while (n_check > 1) {
        if (n_check % 2 != 0) {
            info = -2; // n is not a power of 2
            goto cleanup;
        }
        n_check /= 2;
    }

    char indi_upper = toupper(indi); // Convert character option once  
    if (indi_upper != 'D' && indi_upper != 'I') { info = -1; goto cleanup; }

    // Check pointers for required arrays
    if (xr == NULL) { info = -3; goto cleanup; }
    if (xi == NULL) { info = -4; goto cleanup; }

    // Exit if any validation failed before allocating memory  
    if (info != 0) { goto cleanup; }

    // 3. Memory allocation for column-major copies (if row_major)
    // Calculate the sizes needed for input and output
    size_t input_size = (indi_upper == 'D') ? n : (n + 1); // Size required for input
    size_t output_size = (indi_upper == 'D') ? (n + 1) : n; // Size required for output
    size_t max_size = (input_size > output_size) ? input_size : output_size;
    
    if (row_major) {
        xr_cm = (double*)malloc(max_size * sizeof(double));
        CHECK_ALLOC(xr_cm);
        
        xi_cm = (double*)malloc(max_size * sizeof(double));
        CHECK_ALLOC(xi_cm);
        
        // Copy input data to column-major arrays
        for (size_t i = 0; i < input_size; i++) {
            xr_cm[i] = xr[i];
            xi_cm[i] = xi[i];
        }
    }

    // 4. Prepare Fortran parameters
    double* xr_ptr = row_major ? xr_cm : xr;
    double* xi_ptr = row_major ? xi_cm : xi;
    const int indi_len = 1;

    // 5. Call Fortran function  
    F77_FUNC(dg01nd, DG01ND)(
        &indi_upper,
        &n,
        xr_ptr,
        xi_ptr,
        &info,
        indi_len
    );

    // 6. Copy results back to row-major arrays if needed
    if (row_major && info == 0) {
        // Copy results back to the original arrays
        for (size_t i = 0; i < output_size; i++) {
            xr[i] = xr_ptr[i];
            xi[i] = xi_ptr[i];
        }
    }

cleanup:  
    /* --- Cleanup --- */  
    // Free temporary arrays if allocated  
    free(xr_cm);  
    free(xi_cm);

    // Check if info was set by CHECK_ALLOC during allocation  
    if (info == SLICOT_MEMORY_ERROR) {  
       fprintf(stderr, "Error: Memory allocation failed in slicot_dg01nd.\n");  
    }  
    
    return info;  
}
