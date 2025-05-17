/**
 * @file sg03ad.c
 * @brief C wrapper for SLICOT routine SG03AD.
 * @details Solves continuous- or discrete-time generalized Lyapunov equations
 * and estimates separation.
 * Matrices A, E, Q, Z, X are input/output.
 * Workspace (IWORK, DWORK) is allocated internally.
 * Input/output matrix format is handled via the row_major parameter.
 */

#include <stdlib.h> // For malloc, free
#include <string.h> // For memcpy
#include <ctype.h>  // For toupper
#include <math.h>   // For fmax (from C standard library)

#include "sg03ad.h"       // Public header for this wrapper
#include "slicot_utils.h" // Provides CHECK_ALLOC, SLICOT_MEMORY_ERROR, MAX, MIN, transpose functions etc.
#include "slicot_f77.h"   // Provides F77_FUNC macro

/* External Fortran routine declaration */
extern void F77_FUNC(sg03ad, SG03AD)(
    const char* dico, const char* job, const char* fact, const char* trans, const char* uplo,
    const int* n,
    double* a, const int* lda, /* input/output */
    double* e, const int* lde, /* input/output */
    double* q, const int* ldq, /* input/output */
    double* z, const int* ldz, /* input/output */
    double* x, const int* ldx, /* input: Y, output: X */
    double* scale, double* sep, double* ferr, /* output */
    double* alphar, double* alfai, double* beta, /* output */
    int* iwork, double* dwork, const int* ldwork, /* workspace */
    int* info, /* output */
    size_t dico_len, size_t job_len, size_t fact_len, size_t trans_len, size_t uplo_len
);

SLICOT_EXPORT
int slicot_sg03ad(
    char dico_param, char job_param, char fact_param, char trans_param, char uplo_param,
    int n_param,
    double* a_io, int lda,
    double* e_io, int lde,
    double* q_io, int ldq, /* If FACT='N', output. If FACT='F', input. */
    double* z_io, int ldz, /* If FACT='N', output. If FACT='F', input. */
    double* y_rhs_in_x_out, int ldx, /* Input: Y, Output: X solution */
    double* scale_out, double* sep_out, double* ferr_out,
    double* alphar_out, double* alfai_out, double* beta_out,
    int row_major)
{
    // 1. Variable declarations
    int info = 0;
    int *iwork = NULL;
    double *dwork = NULL;
    int liwork = 0;
    int ldwork = 0;

    // Column-major copies for input/output matrices
    double *a_cm = NULL, *e_cm = NULL, *q_cm = NULL, *z_cm = NULL, *x_cm = NULL;

    // Fortran-style leading dimensions
    int lda_f, lde_f, ldq_f, ldz_f, ldx_f;
    
    // Character parameters for Fortran
    char dico_f = toupper(dico_param);
    char job_f  = toupper(job_param);
    char fact_f = toupper(fact_param);
    char trans_f= toupper(trans_param);
    char uplo_f = toupper(uplo_param);

    // 2. Input parameter validation
    if (strchr("CD", dico_f) == NULL) { info = -1; goto cleanup; }
    if (strchr("XSB", job_f) == NULL) { info = -2; goto cleanup; }
    if (strchr("NF", fact_f) == NULL) { info = -3; goto cleanup; }
    if (strchr("NT", trans_f) == NULL) { info = -4; goto cleanup; }
    if (strchr("LU", uplo_f) == NULL) { info = -5; goto cleanup; }
    if (n_param < 0) { info = -6; goto cleanup; }

    // Validate pointers for required arrays
    // A, E are always N x N
    if (a_io == NULL && n_param > 0) { info = -7; goto cleanup; }
    if (e_io == NULL && n_param > 0) { info = -9; goto cleanup; }
    // Q, Z are N x N. If FACT='F', they are input. If FACT='N', they are output.
    if (q_io == NULL && n_param > 0) { info = -11; goto cleanup; }
    if (z_io == NULL && n_param > 0) { info = -13; goto cleanup; }
    // X (input Y) is N x N, required if JOB is 'X' or 'B'
    if ((job_f == 'X' || job_f == 'B') && y_rhs_in_x_out == NULL && n_param > 0) { info = -15; goto cleanup; }
    
    // Output scalars/arrays
    if (scale_out == NULL) { info = -17; goto cleanup; } // SCALE is always output
    if ((job_f == 'S' || job_f == 'B') && sep_out == NULL) { info = -18; goto cleanup; } // SEP output if JOB='S' or 'B'
    if (job_f == 'B' && ferr_out == NULL) { info = -19; goto cleanup; } // FERR output if JOB='B'
    // ALPHAR, ALPHAI, BETA are output if FACT='N'
    if (fact_f == 'N' && n_param > 0) {
        if (alphar_out == NULL) { info = -20; goto cleanup; }
        if (alfai_out == NULL) { info = -21; goto cleanup; }
        if (beta_out == NULL) { info = -22; goto cleanup; }
    }
    
    // Validate leading dimensions
    int min_ld_f = MAX(1, n_param); // For A, E, Q, Z, X

    if (row_major) { // C LDA is number of columns
        if (n_param > 0 && lda < n_param) { info = -8; goto cleanup; }
        if (n_param > 0 && lde < n_param) { info = -10; goto cleanup; }
        if (n_param > 0 && ldq < n_param) { info = -12; goto cleanup; }
        if (n_param > 0 && ldz < n_param) { info = -14; goto cleanup; }
        if (n_param > 0 && (job_f == 'X' || job_f == 'B') && ldx < n_param) { info = -16; goto cleanup; }
    } else { // Column-major C (Fortran-style LDs)
        if (n_param > 0 && lda < min_ld_f) { info = -8; goto cleanup; }
        if (n_param > 0 && lde < min_ld_f) { info = -10; goto cleanup; }
        if (n_param > 0 && ldq < min_ld_f) { info = -12; goto cleanup; }
        if (n_param > 0 && ldz < min_ld_f) { info = -14; goto cleanup; }
        if (n_param > 0 && (job_f == 'X' || job_f == 'B') && ldx < min_ld_f) { info = -16; goto cleanup; }
    }
    if (info != 0) { goto cleanup; }

    // 3. Internal Workspace Allocation
    // IWORK(N*N)
    liwork = (n_param > 0) ? (size_t)n_param * n_param : 1;
    liwork = MAX(1, liwork); // Ensure at least 1
    iwork = (int*)malloc(liwork * sizeof(int));
    CHECK_ALLOC(iwork);

    // LDWORK query
    int ldwork_query_val = -1;
    double dwork_query_result[1];
    F77_FUNC(sg03ad, SG03AD)(&dico_f, &job_f, &fact_f, &trans_f, &uplo_f, &n_param,
                             NULL, &min_ld_f, NULL, &min_ld_f, NULL, &min_ld_f, NULL, &min_ld_f, NULL, &min_ld_f,
                             scale_out, sep_out, ferr_out, NULL, NULL, NULL,
                             iwork, dwork_query_result, &ldwork_query_val, &info,
                             1,1,1,1,1);
    
    if (info == 0 && ldwork_query_val == -1) { // Query successful
        ldwork = (int)dwork_query_result[0];
    } else { // Fallback to formula if query fails or not supported as expected
        info = 0; // Reset info if it was from query with -1
        if (job_f == 'X') {
            ldwork = (fact_f == 'F') ? MAX(1, n_param) : MAX(1, 4 * n_param);
        } else { // JOB = 'B' or 'S'
            ldwork = (fact_f == 'F') ? MAX(1, 2 * n_param * n_param) : MAX(1, MAX(2 * n_param * n_param, 4 * n_param));
        }
    }
    ldwork = MAX(1, ldwork);
    dwork = (double*)malloc((size_t)ldwork * sizeof(double));
    CHECK_ALLOC(dwork);

    // 4. Memory allocation for column-major copies
    size_t n_sq = (size_t)n_param * n_param;
    if (n_param == 0) n_sq = 0;

    if (row_major && n_param > 0) {
        a_cm = (double*)malloc(n_sq * sizeof(double)); CHECK_ALLOC(a_cm);
        e_cm = (double*)malloc(n_sq * sizeof(double)); CHECK_ALLOC(e_cm);
        q_cm = (double*)malloc(n_sq * sizeof(double)); CHECK_ALLOC(q_cm);
        z_cm = (double*)malloc(n_sq * sizeof(double)); CHECK_ALLOC(z_cm);
        if (job_f == 'X' || job_f == 'B') {
            x_cm = (double*)malloc(n_sq * sizeof(double)); CHECK_ALLOC(x_cm);
        }
    }

    // 5. Prepare Fortran parameters and perform conversions
    double* a_ptr = a_io; double* e_ptr = e_io; double* q_ptr = q_io;
    double* z_ptr = z_io; double* x_ptr = y_rhs_in_x_out;

    lda_f = lda; lde_f = lde; ldq_f = ldq; ldz_f = ldz; ldx_f = ldx;

    if (row_major && n_param > 0) {
        lda_f = min_ld_f; lde_f = min_ld_f; ldq_f = min_ld_f;
        ldz_f = min_ld_f; ldx_f = min_ld_f;

        slicot_transpose_to_fortran_with_ld(a_io, a_cm, n_param, n_param, lda, lda_f, sizeof(double)); a_ptr = a_cm;
        slicot_transpose_to_fortran_with_ld(e_io, e_cm, n_param, n_param, lde, lde_f, sizeof(double)); e_ptr = e_cm;
        
        if (fact_f == 'F') { // Q and Z are inputs if FACT='F'
            slicot_transpose_to_fortran_with_ld(q_io, q_cm, n_param, n_param, ldq, ldq_f, sizeof(double)); q_ptr = q_cm;
            slicot_transpose_to_fortran_with_ld(z_io, z_cm, n_param, n_param, ldz, ldz_f, sizeof(double)); z_ptr = z_cm;
        } else { // Q and Z are outputs if FACT='N', Fortran writes to _cm buffers
            q_ptr = q_cm;
            z_ptr = z_cm;
        }
        
        if (job_f == 'X' || job_f == 'B') {
            slicot_transpose_to_fortran_with_ld(y_rhs_in_x_out, x_cm, n_param, n_param, ldx, ldx_f, sizeof(double)); x_ptr = x_cm;
        } else { // JOB = 'S', X is not referenced as input
            x_ptr = x_cm; // Still point to cm buffer for potential output, though not strictly defined as output if JOB='S'
        }
    } else { // Column-major C or N=0
        if (n_param == 0) { // For N=0, pass NULL if original was NULL, or original if non-NULL but with LD=1
            a_ptr = (a_io == NULL) ? NULL : a_io; lda_f = (a_io == NULL) ? 1: lda;
            e_ptr = (e_io == NULL) ? NULL : e_io; lde_f = (e_io == NULL) ? 1: lde;
            q_ptr = (q_io == NULL) ? NULL : q_io; ldq_f = (q_io == NULL) ? 1: ldq;
            z_ptr = (z_io == NULL) ? NULL : z_io; ldz_f = (z_io == NULL) ? 1: ldz;
            x_ptr = (y_rhs_in_x_out == NULL && (job_f == 'X' || job_f == 'B')) ? NULL : y_rhs_in_x_out;
            ldx_f = (y_rhs_in_x_out == NULL && (job_f == 'X' || job_f == 'B')) ? 1 : ldx;
        }
    }
    
    int n_f_call = n_param;
    int ldwork_f_call = ldwork;

    F77_FUNC(sg03ad, SG03AD)(&dico_f, &job_f, &fact_f, &trans_f, &uplo_f, &n_f_call,
                             a_ptr, &lda_f, e_ptr, &lde_f, q_ptr, &ldq_f, z_ptr, &ldz_f, x_ptr, &ldx_f,
                             scale_out, sep_out, ferr_out, alphar_out, alfai_out, beta_out,
                             iwork, dwork, &ldwork_f_call, &info,
                             1,1,1,1,1); // Lengths of char params

    // 8. Convert results back to row-major
    if (row_major && info == 0 && n_param > 0) { // Only if N > 0
        // A, E, Q, Z, X are modified
        if (a_cm != NULL && a_io != NULL) slicot_transpose_to_c_with_ld(a_cm, a_io, n_param, n_param, lda_f, lda, sizeof(double));
        if (e_cm != NULL && e_io != NULL) slicot_transpose_to_c_with_ld(e_cm, e_io, n_param, n_param, lde_f, lde, sizeof(double));
        if (q_cm != NULL && q_io != NULL) slicot_transpose_to_c_with_ld(q_cm, q_io, n_param, n_param, ldq_f, ldq, sizeof(double)); // Q is output if FACT='N'
        if (z_cm != NULL && z_io != NULL) slicot_transpose_to_c_with_ld(z_cm, z_io, n_param, n_param, ldz_f, ldz, sizeof(double)); // Z is output if FACT='N'
        
        if ((job_f == 'X' || job_f == 'B') && x_cm != NULL && y_rhs_in_x_out != NULL) {
             slicot_transpose_to_c_with_ld(x_cm, y_rhs_in_x_out, n_param, n_param, ldx_f, ldx, sizeof(double));
        }
    }

cleanup:
    free(iwork);
    free(dwork);
    // bwork_c was not used in this wrapper as SG03AD does not have BWORK argument.
    // free(bwork_c); 

    if (row_major) {
        free(a_cm); free(e_cm); free(q_cm); free(z_cm); free(x_cm);
    }
    
    if (info == SLICOT_MEMORY_ERROR) {
       // Memory allocation error
    }
    return info;
}
