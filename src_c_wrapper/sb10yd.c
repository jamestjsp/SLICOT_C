/**
 * @file sb10yd.c
 * @brief C wrapper for SLICOT routine SB10YD.
 * @details Fits supplied frequency response data with a stable, minimum phase SISO system.
 * Workspace (IWORK, DWORK, ZWORK) is allocated internally.
 * Output matrices A, B, C, D and output scalar N are modified.
 * Input/output matrix format for A is handled via the row_major parameter.
 * B, C are vectors, D is scalar-like. RFRDAT, IFRDAT, OMEGA are 1D inputs.
 */

#include <stdlib.h> // For malloc, free
#include <string.h> // For memcpy
#include <math.h>   // For fmax, fmin (from C standard library)
#include "sb10yd.h" // Public header for this wrapper
#include "slicot_utils.h"  // Provides CHECK_ALLOC, SLICOT_MEMORY_ERROR, MAX, MIN, transpose functions etc.
#include "slicot_f77.h"    // Provides F77_FUNC macro for Fortran name mangling

/* External Fortran routine declaration */
extern void F77_FUNC(sb10yd, SB10YD)(
    const int* discfl, const int* flag, const int* lendat,
    const double* rfrdat, const double* ifrdat, const double* omega,
    int* n_fortran, /* input/output */
    double* a_fortran, const int* lda_fortran,
    double* b_fortran, double* c_fortran, double* d_fortran,
    const double* tol,
    int* iwork, double* dwork, const int* ldwork,
    slicot_complex_double* zwork, const int* lzwork, /* COMPLEX*16 */
    int* info
);

/* C wrapper function definition */
SLICOT_EXPORT
int slicot_sb10yd(int discfl, int flag, int lendat,
                  const double* rfrdat, const double* ifrdat, const double* omega,
                  int* n_io, /* input/output */
                  double* a, int lda,
                  double* b, double* c, double* d,
                  double tol_param,
                  int row_major)
{
    // 1. Variable declarations
    int info = 0;
    int *iwork = NULL;
    double *dwork = NULL;
    slicot_complex_double *zwork = NULL; // For Fortran COMPLEX*16
    
    int liwork = 0;
    int ldwork = 0;
    int lzwork = 0;

    double *a_cm = NULL; // For column-major copy of A if row_major is true

    // Fortran-style leading dimension for A
    int lda_f;
    double tol_f = tol_param; // Pass tol by reference

    // Initial value of N for workspace calculation and Fortran call
    int n_val_in = (n_io != NULL) ? *n_io : 0;


    // 2. Input parameter validation
    if (discfl != 0 && discfl != 1) { info = -1; goto cleanup; }
    if (flag != 0 && flag != 1) { info = -2; goto cleanup; }
    if (lendat < 2) { info = -3; goto cleanup; }
    if (rfrdat == NULL && lendat > 0) { info = -4; goto cleanup; }
    if (ifrdat == NULL && lendat > 0) { info = -5; goto cleanup; }
    if (omega == NULL && lendat > 0) { info = -6; goto cleanup; }
    if (n_io == NULL) { info = -7; goto cleanup; }
    if (n_val_in < 0 || n_val_in > lendat -1) {info = -7; goto cleanup;}


    // A, B, C, D are outputs. Pointers must be valid if n_val_in (on entry) > 0 for A, B, C.
    // D is always output (scalar).
    if (a == NULL && n_val_in > 0) { info = -8; goto cleanup; }
    if (b == NULL && n_val_in > 0) { info = -10; goto cleanup; } // B is N-vector
    if (c == NULL && n_val_in > 0) { info = -11; goto cleanup; } // C is N-vector
    if (d == NULL) { info = -12; goto cleanup; } // D is scalar

    // Check leading dimension for A if N_in > 0
    int min_lda_f = MAX(1, n_val_in);
    if (n_val_in > 0) {
        if (row_major) { // C LDA is number of columns
            if (lda < n_val_in) { info = -9; goto cleanup; }
        } else { // Column-major C (Fortran-style LDA)
            if (lda < min_lda_f) { info = -9; goto cleanup; }
        }
    } else { // N_in == 0
        if (lda < 1 && a != NULL) {info = -9; goto cleanup;} // LDA must be >=1 even if N=0 for non-NULL A
    }
    if (info != 0) { goto cleanup; }


    // 3. Internal Workspace Allocation
    // LIWORK = max(2, 2*N_in+1)
    liwork = MAX(2, 2 * n_val_in + 1);
    iwork = (int*)malloc((size_t)liwork * sizeof(int));
    CHECK_ALLOC(iwork);

    // LDWORK calculation
    int hpts = 2048; // Constant from documentation
    int lw1 = 2 * lendat + 4 * hpts;
    int lw2 = lendat + 6 * hpts;
    int mn_ws = MIN(2 * lendat, 2 * n_val_in + 1);
    int lw3;
    if (n_val_in > 0) {
        lw3 = 2 * lendat * (2 * n_val_in + 1) + MAX(2 * lendat, 2 * n_val_in + 1) +
              MAX(mn_ws + 6 * n_val_in + 4, 2 * mn_ws + 1);
    } else { // N_in == 0
        lw3 = 4 * lendat + 5;
    }
    int lw4 = 0;
    if (flag == 1) {
        lw4 = MAX(n_val_in * n_val_in + 5 * n_val_in, 6 * n_val_in + 1 + MIN(1, n_val_in));
    }
    ldwork = MAX(2, MAX(lw1, MAX(lw2, MAX(lw3, lw4))));
    
    dwork = (double*)malloc((size_t)ldwork * sizeof(double));
    CHECK_ALLOC(dwork);

    // LZWORK calculation
    if (n_val_in > 0) {
        lzwork = lendat * (2 * n_val_in + 3);
    } else { // N_in == 0
        lzwork = lendat;
    }
    lzwork = MAX(1, lzwork); // Ensure at least 1 if lendat is small and N=0
                             // (though lendat >= 2 is a precondition)
    
    zwork = (slicot_complex_double*)malloc((size_t)lzwork * sizeof(slicot_complex_double));
    CHECK_ALLOC(zwork);


    // 4. Memory allocation for column-major copy of A (output)
    // Size based on initial N, as N can change.
    // The Fortran routine will write to a_ptr. If row_major, a_ptr is a_cm.
    // The final N (N_out) determines how much to copy back.
    size_t a_size_elem = (size_t)n_val_in * n_val_in;
    if (n_val_in == 0) a_size_elem = 0; // If N_in is 0, A is not dimensioned by N_in.

    if (row_major && n_val_in > 0) { // Only needed if N_in > 0 and row_major
        a_cm = (double*)malloc(a_size_elem * sizeof(double));
        CHECK_ALLOC(a_cm);
    }

    // 5. Prepare Fortran parameters
    double* a_ptr = a;
    lda_f = lda;

    if (row_major) {
        if (n_val_in > 0) { // If N_in is 0, a_ptr remains 'a', lda_f remains 'lda' (must be >=1)
            lda_f = MAX(1, n_val_in);
            a_ptr = a_cm; // Fortran writes to a_cm
        } else { // N_in == 0
            lda_f = MAX(1, lda); // Use user's LDA, must be >=1
            a_ptr = a; // Fortran might not touch A if N=0, or expect a dummy placeholder
        }
    } else { // Column-major C
        if (n_val_in == 0) {
             a_ptr = a; // Pass original 'a', LDA must be >=1
             lda_f = MAX(1, lda);
        } else {
            a_ptr = a;
            lda_f = lda; // Already Fortran-style
        }
    }
    
    // N is input/output, pass directly to Fortran
    int n_for_fortran = n_val_in;
    int ldwork_f_call = ldwork;
    int lzwork_f_call = lzwork;


    // 7. Call Fortran function
    F77_FUNC(sb10yd, SB10YD)(&discfl, &flag, &lendat,
                             rfrdat, ifrdat, omega,
                             &n_for_fortran, /* n_io is passed directly */
                             a_ptr, &lda_f,
                             b, c, d, /* B, C are N-vectors, D is scalar, passed directly */
                             &tol_f,
                             iwork, dwork, &ldwork_f_call,
                             zwork, &lzwork_f_call,
                             &info);

    // Update C n_io from Fortran output N
    if (n_io != NULL) {
        *n_io = n_for_fortran;
    }
    int n_val_out = n_for_fortran; // The actual order of the system identified.

    // 8. Convert result A back to row-major if needed
    if (row_major && info == 0) {
        if (n_val_out > 0 && a_cm != NULL && a != NULL) {
            // Transpose from a_cm (where Fortran wrote) to user's 'a'.
            // The dimensions for transpose are N_out x N_out.
            // User's 'a' array must be large enough for N_out x N_out with C_LDA 'lda'.
            // Fortran LDA for a_cm was MAX(1, n_val_in).
            // C LDA for 'a' is 'lda'.
            slicot_transpose_to_c_with_ld(a_cm, a, n_val_out, n_val_out, MAX(1, n_val_in) /*lda_f used by fortran*/, lda /*user's lda*/, sizeof(double));
        } else if (n_val_out == 0 && n_val_in > 0 && a != NULL) {
            // If system becomes 0-order, and original N was >0, what to do with 'a'?
            // Conventionally, A would be empty. The caller should check N_out.
            // No transpose needed.
        }
    }
    
    // DWORK(1) and DWORK(2) might contain optimal LDWORK and LZWORK if INFO=0
    // This wrapper uses calculated values, but this info could be returned if desired.

cleanup:
    free(iwork);
    free(dwork);
    free(zwork);
    if (row_major) {
        free(a_cm);
    }
    
    if (info == SLICOT_MEMORY_ERROR) {
       // Memory allocation error was caught by CHECK_ALLOC
    }
    return info;
}
