/**
 * @file mc01td.c
 * @brief C wrapper implementation for SLICOT routine MC01TD
 *
 * This file provides a C wrapper implementation for the SLICOT routine MC01TD,
 * which checks the stability of a given real polynomial.
 */

 #include <stdlib.h>
 #include <ctype.h>
 #include <stddef.h> // For size_t
 
 // Include the header file for this wrapper
 #include "mc01td.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions
 
 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * Note: Fortran LOGICAL maps to int* in C.
  */
 extern void F77_FUNC(mc01td, MC01TD)(
     const char* dico,       // CHARACTER*1 DICO
     int* dp,                // INTEGER DP (in/out)
     const double* p,        // DOUBLE PRECISION P(*)
     int* stable,            // LOGICAL STABLE (output) -> int*
     int* nz,                // INTEGER NZ (output)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     int* iwarn,             // INTEGER IWARN (output)
     int* info,              // INTEGER INFO (output)
     int dico_len            // Hidden length
 );
 
 
 /* C wrapper function definition */
 SLICOT_C_WRAPPER_API
 int slicot_mc01td(char dico, int* dp, const double* p,
                   int* stable, int* nz, int* iwarn)
 {
     /* Local variables */
     int info = 0;
     double* dwork = NULL; // Workspace
     int dwork_size = 0;
     const int dico_len = 1;
     int dp_val = (dp != NULL) ? *dp : -1; // Get degree value for validation/workspace
 
     char dico_upper = toupper(dico);
 
     /* --- Input Parameter Validation --- */
 
     if (dp == NULL || p == NULL || stable == NULL || nz == NULL || iwarn == NULL) {
         // Check for NULL pointers for essential arguments
         info = -99; // Assign a custom error code for NULL pointers
         goto cleanup;
     }
     if (dp_val < 0) { info = -2; goto cleanup; }
     if (dico_upper != 'C' && dico_upper != 'D') { info = -1; goto cleanup; }
 
 
     /* --- Workspace Allocation --- */
 
     // Allocate DWORK (size 2*DP+2) - No query needed
     // Use dp_val (degree on entry) for allocation size
     dwork_size = MAX(1, 2 * dp_val + 2); // Ensure minimum size 1
     dwork = (double*)malloc((size_t)dwork_size * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure
 
     /* --- Prepare Arrays and Call Fortran Routine --- */
 
     // Arrays P is 1D input, STABLE, NZ, IWARN are scalar outputs. No row-major needed.
 
     /* Call the Fortran routine directly */
     F77_FUNC(mc01td, MC01TD)(&dico_upper, dp, p, stable, nz,
                              dwork, iwarn, &info, dico_len);
     // DP, STABLE, NZ, IWARN are modified in place.
 
 cleanup:
     /* --- Cleanup --- */
     free(dwork);
 
     /* Return the info code from the Fortran routine or SLICOT_MEMORY_ERROR */
     return info;
 }
 