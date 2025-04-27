/**  
 * @file dk01md.c  
 * @brief C wrapper for SLICOT routine DK01MD.  
 * @details Apply an anti-aliasing window to a real signal.
 * Workspace is allocated internally by this wrapper.
 * Input/output format is handled via the row_major parameter.  
 */

#include <stdlib.h> // For malloc, free  
#include <ctype.h>  // For toupper  
#include <stddef.h> // For size_t  
#include <stdio.h>  // For error logging

#include "dk01md.h" // Public header for this wrapper  
#include "slicot_utils.h"  // Provides CHECK_ALLOC, SLICOT_MEMORY_ERROR, etc.  
#include "slicot_f77.h"    // Provides F77_FUNC macro for Fortran name mangling

/* External Fortran routine declaration */  
extern void F77_FUNC(dk01md, DK01MD)(  
    const char* type,
    const int* n,
    double* a,
    int* info,
    int type_len  // Hidden string length
);

/* C wrapper function definition */  
SLICOT_EXPORT
int slicot_dk01md(char type, int n, double* a, int row_major)  
{  
    // 1. Variable declarations  
    int info = 0;
    
    // Pointers for column-major copies (if row_major is used)  
    double* a_cm = NULL;

    // 2. Input parameter validation (Check BEFORE allocating memory)
    if (n < 1) { info = -2; goto cleanup; } // Check dimension n
    if (a == NULL) { info = -3; goto cleanup; } // Check array pointer
    
    char type_upper = toupper(type);
    if (type_upper != 'M' && type_upper != 'N' && type_upper != 'Q') {
        info = -1; goto cleanup;
    }

    // Exit if any validation failed before allocating memory  
    if (info != 0) { goto cleanup; }

    // 3. Memory allocation for column-major copies (if row_major)
    if (row_major) {
        a_cm = (double*)malloc((size_t)n * sizeof(double));
        CHECK_ALLOC(a_cm);
        
        // Copy input data to column-major array
        for (size_t i = 0; i < (size_t)n; i++) {
            a_cm[i] = a[i];
        }
    }

    // 4. Prepare Fortran parameters
    double* a_ptr = row_major ? a_cm : a;
    const int type_len = 1;

    // 5. Call Fortran function  
    F77_FUNC(dk01md, DK01MD)(
        &type_upper,
        &n,
        a_ptr,
        &info,
        type_len
    );

    // 6. Copy results back to row-major arrays if needed
    if (row_major && info == 0) {
        // Copy results back to the original array
        for (size_t i = 0; i < (size_t)n; i++) {
            a[i] = a_ptr[i];
        }
    }

cleanup:  
    /* --- Cleanup --- */  
    // Free temporary array if allocated  
    free(a_cm);

    // Check if info was set by CHECK_ALLOC during allocation  
    if (info == SLICOT_MEMORY_ERROR) {  
       fprintf(stderr, "Error: Memory allocation failed in slicot_dk01md.\n");  
    }  
    
    return info;  
}
