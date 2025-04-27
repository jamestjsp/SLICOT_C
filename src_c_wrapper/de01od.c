/**
 * @file de01od.c
 * @brief C wrapper for SLICOT routine DE01OD.
 * @details This routine computes the convolution or deconvolution of two real signals.
 * The algorithm uses an FFT approach. No internal workspace is needed as
 * the Fortran routine does not require additional workspace arrays.
 * Input/output matrix format is handled via the row_major parameter.
 */

#include <stdlib.h>
#include <ctype.h>
#include <stddef.h>
#include <stdio.h>

#include "de01od.h"
#include "slicot_utils.h"
#include "slicot_f77.h"

/* External Fortran routine declaration */
extern void F77_FUNC(de01od, DE01OD)(
    const char* conv,
    const int* n,
    double* a,
    double* b,
    int* info,
    int conv_len
);

/* C wrapper function definition */
SLICOT_EXPORT
int slicot_de01od(char conv, int n, double* a, double* b, int row_major)
{
    int info = 0;
    
    // Pointers for column-major copies (if row_major is used)
    double* a_cm = NULL;
    double* b_cm = NULL;
    
    // 1. Input parameter validation
    // Check scalar parameters first
    if (n < 2) { info = -2; goto cleanup; } // N must be >= 2
    
    // Check if N is a power of 2
    int power_of_2 = 1;
    while (power_of_2 < n) {
        power_of_2 <<= 1;
    }
    if (power_of_2 != n) { info = -2; goto cleanup; } // N must be a power of 2
    
    // Check option parameter
    char conv_upper = toupper(conv);
    if (conv_upper != 'C' && conv_upper != 'D') { info = -1; goto cleanup; } // CONV must be 'C' or 'D'
    
    // Check pointers for required arrays
    if (a == NULL) { info = -3; goto cleanup; } // A cannot be NULL
    if (b == NULL) { info = -4; goto cleanup; } // B cannot be NULL
    
    // 2. Memory allocation for column-major copies (if row_major)
    if (row_major) {
        a_cm = (double*)malloc((size_t)n * sizeof(double));
        CHECK_ALLOC(a_cm);
        
        b_cm = (double*)malloc((size_t)n * sizeof(double));
        CHECK_ALLOC(b_cm);
        
        // Copy data to column-major format (in this case just copying the arrays)
        for (int i = 0; i < n; i++) {
            a_cm[i] = a[i];
            b_cm[i] = b[i];
        }
    }
    
    // 3. Prepare Fortran parameters
    double* a_ptr = row_major ? a_cm : a;
    double* b_ptr = row_major ? b_cm : b;
    
    // 4. Call Fortran function
    F77_FUNC(de01od, DE01OD)(
        &conv_upper,
        &n,
        a_ptr,
        b_ptr,
        &info,
        1 /* conv_len = 1 character */
    );
    
    // 5. Copy results back to row-major format if needed
    if (row_major && info == 0) {
        for (int i = 0; i < n; i++) {
            a[i] = a_cm[i];
        }
    }
    
cleanup:
    free(a_cm);
    free(b_cm);
    
    if (info == SLICOT_MEMORY_ERROR) {
        fprintf(stderr, "Error: Memory allocation failed in slicot_de01od.\n");
    }
    
    return info;
}
