/**
 * @file mc01td.c
 * @brief C wrapper implementation for SLICOT routine MC01TD
 *
 * This file provides a C wrapper implementation for the SLICOT routine MC01TD,
 * which checks the stability of a given real polynomial.
 * NOTE: This routine does not involve 2D arrays requiring row/column major handling.
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
     const double* p,        // DOUBLE PRECISION P(*) (in)
     int* stable,            // LOGICAL STABLE (output) -> int*
     int* nz,                // INTEGER NZ (output)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     int* iwarn,             // INTEGER IWARN (output)
     int* info,              // INTEGER INFO (output)
     int dico_len            // Hidden length
 );


 /* C wrapper function definition */
 SLICOT_EXPORT
 int slicot_mc01td(char dico, int* dp, const double* p,
                   int* stable, int* nz, int* iwarn)
                   // No row_major parameter needed
 {
     /* Local variables */
     int info = 0;
     int ldwork = 0; // Size calculated directly
     double* dwork = NULL; // Workspace
     // No iwork needed for this routine
     const int dico_len = 1;
     int dp_val = (dp != NULL) ? *dp : -1; // Get degree value for validation/workspace

     char dico_upper = toupper(dico);
     // No _cm pointers needed.

     /* --- Input Parameter Validation --- */

     if (dico_upper != 'C' && dico_upper != 'D') { info = -1; goto cleanup; }
     // Check for NULL pointers for essential arguments AFTER checking DICO
     if (dp == NULL || p == NULL || stable == NULL || nz == NULL || iwarn == NULL) {
         // Assign a custom error code for NULL pointers if needed, or let info=-2 handle dp==NULL
         info = -99; // Example custom error code
         goto cleanup;
     }
     if (dp_val < 0) { info = -2; goto cleanup; }


     /* --- Workspace Allocation --- */

     // Allocate DWORK (size 2*DP+2) - No query needed
     // Use dp_val (degree on entry) for allocation size
     ldwork = MAX(1, 2 * dp_val + 2); // Ensure minimum size 1
     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure

     /* --- Call the computational routine --- */

     // Array P is 1D input. DP, STABLE, NZ, IWARN are scalar outputs.
     // No row-major conversion needed.
     F77_FUNC(mc01td, MC01TD)(&dico_upper, dp, p, stable, nz,
                              dwork, iwarn, &info, dico_len);
     // DP, STABLE, NZ, IWARN are modified in place.

     /* --- Copy results back to row-major format if needed --- */
     // Not applicable for this routine.

 cleanup:
     /* --- Cleanup --- */
     free(dwork);

     /* Return the info code from the Fortran routine or SLICOT_MEMORY_ERROR */
     return info;
 }