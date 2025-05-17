/**
 * @file sb02mt.c
 * @brief C wrapper implementation for SLICOT routine SB02MT
 *
 * This file provides a C wrapper implementation for the SLICOT routine SB02MT,
 * which converts optimal problems with coupling weighting terms to standard problems.
 * Workspace (IWORK, DWORK) is allocated internally by this wrapper.
 * Input/output matrix format is handled via the row_major parameter.
 */

#include <stdlib.h>
#include <ctype.h>  // For toupper
#include <stddef.h> // For size_t
#include <math.h>   // For MAX/MIN if not in slicot_utils.h

// Include the header file for this wrapper
#include "sb02mt.h"
// Include necessary SLICOT utility headers
#include "slicot_utils.h" // Assumed to contain MAX, MIN, CHECK_ALLOC, transpose etc.
#include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

/**
 * Declare the external Fortran routine using the F77_FUNC macro.
 */
extern void F77_FUNC(sb02mt, SB02MT)(
    const char* jobg,       // CHARACTER*1 JOBG
    const char* jobl,       // CHARACTER*1 JOBL
    const char* fact,       // CHARACTER*1 FACT
    const char* uplo,       // CHARACTER*1 UPLO
    const int* n,           // INTEGER N
    const int* m,           // INTEGER M
    double* a,              // DOUBLE PRECISION A(LDA,*)
    const int* lda,         // INTEGER LDA
    double* b,              // DOUBLE PRECISION B(LDB,*)
    const int* ldb,         // INTEGER LDB
    double* q,              // DOUBLE PRECISION Q(LDQ,*)
    const int* ldq,         // INTEGER LDQ
    double* r,              // DOUBLE PRECISION R(LDR,*)
    const int* ldr,         // INTEGER LDR
    double* l,              // DOUBLE PRECISION L(LDL,*)
    const int* ldl,         // INTEGER LDL
    int* ipiv,              // INTEGER IPIV(*)
    int* oufact,            // INTEGER OUFACT
    double* g,              // DOUBLE PRECISION G(LDG,*)
    const int* ldg,         // INTEGER LDG
    int* iwork,             // INTEGER IWORK(*)
    double* dwork,          // DOUBLE PRECISION DWORK(*)
    const int* ldwork,      // INTEGER LDWORK
    int* info,              // INTEGER INFO (output)
    int jobg_len,           // Hidden length for jobg
    int jobl_len,           // Hidden length for jobl
    int fact_len,           // Hidden length for fact
    int uplo_len            // Hidden length for uplo
);


 /* C wrapper function definition */
 SLICOT_EXPORT
 int slicot_sb02mt(char jobg, char jobl, char fact, char uplo,
                   int n, int m,
                   double* a_io, int lda,
                   double* b_io, int ldb,
                   double* q_io, int ldq,
                   double* r_io, int ldr,
                   double* l_io, int ldl,
                   int* ipiv_io,
                   double* g_out, int ldg,
                   int* oufact_out,
                   int row_major)
 {
     /* Local variables */
     int info = 0; // Main info for the routine
     int ldwork_val = 0;
     double* dwork_ptr = NULL;
     int* iwork_ptr = NULL;
     int iwork_size_val = 0; // For IWORK array in Fortran
     const int jobg_len = 1, jobl_len = 1, fact_len = 1, uplo_len = 1;

     char jobg_upper = toupper(jobg);
     char jobl_upper = toupper(jobl);
     char fact_upper = toupper(fact);
     char uplo_upper = toupper(uplo);

     /* Arrays for column-major versions if needed */
     double *a_cm = NULL, *b_cm = NULL, *q_cm = NULL, *r_cm = NULL, *l_cm = NULL, *g_cm = NULL;

     /* Pointers to pass to Fortran */
     double *a_f_ptr, *b_f_ptr, *q_f_ptr, *r_f_ptr, *l_f_ptr, *g_f_ptr;
     int *ipiv_f_ptr;

     int lda_f, ldb_f, ldq_f, ldr_f, ldl_f, ldg_f;

    /* --- Input Parameter Validation (BEFORE allocating memory) --- */
    if (jobg_upper != 'G' && jobg_upper != 'N') { info = -1; goto cleanup; }
    if (jobl_upper != 'Z' && jobl_upper != 'N') { info = -2; goto cleanup; }
    if (fact_upper != 'N' && fact_upper != 'C' && fact_upper != 'U') { info = -3; goto cleanup; }
    if (uplo_upper != 'U' && uplo_upper != 'L') { info = -4; goto cleanup; }
    if (n < 0) { info = -5; goto cleanup; }
    if (m < 0) { info = -6; goto cleanup; }

    // A
    if (jobl_upper == 'N') {
        if (n > 0 && a_io == NULL) { info = -7; goto cleanup; }
        if (row_major) { if (n > 0 && lda < n) { info = -8; goto cleanup; } }
        else { if (n > 0 && lda < n) { info = -8; goto cleanup; } }
    } else { if (lda < 1) { info = -8; goto cleanup; } }
    // B
    if (n > 0 && m > 0 && b_io == NULL) { info = -9; goto cleanup; }
    if (row_major) { if (n > 0 && m > 0 && ldb < m) { info = -10; goto cleanup; } }
    else { if (n > 0 && m > 0 && ldb < n) { info = -10; goto cleanup; } }
    // Q
    if (jobl_upper == 'N') {
        if (n > 0 && q_io == NULL) { info = -11; goto cleanup; }
        if (row_major) { if (n > 0 && ldq < n) { info = -12; goto cleanup; } }
        else { if (n > 0 && ldq < n) { info = -12; goto cleanup; } }
    } else { if (ldq < 1) { info = -12; goto cleanup; } }
    // R
    if (m > 0 && r_io == NULL) { info = -13; goto cleanup; }
    if (row_major) { if (m > 0 && ldr < m) { info = -14; goto cleanup; } }
    else { if (m > 0 && ldr < m) { info = -14; goto cleanup; } }
    // L
    if (jobl_upper == 'N') {
        if (n > 0 && m > 0 && l_io == NULL) { info = -15; goto cleanup; }
        if (row_major) { if (n > 0 && m > 0 && ldl < m) { info = -16; goto cleanup; } }
        else { if (n > 0 && m > 0 && ldl < n) { info = -16; goto cleanup; } }
    } else { if (ldl < 1) { info = -16; goto cleanup; } }
    // IPIV
    if (fact_upper == 'U' && m > 0 && ipiv_io == NULL) { info = -17; goto cleanup; }
    // OUFACT
    if (oufact_out == NULL) { info = -18; goto cleanup; }
    // G
    if (jobg_upper == 'G') {
        if (n > 0 && g_out == NULL) { info = -19; goto cleanup; }
        if (row_major) { if (n > 0 && ldg < n) { info = -20; goto cleanup; } }
        else { if (n > 0 && ldg < n) { info = -20; goto cleanup; } }
    } else { if (ldg < 1) { info = -20; goto cleanup; } }

    if (info != 0) goto cleanup;

    /* --- Workspace Query & Allocation --- */
    // IWORK
    iwork_size_val = (fact_upper == 'N' && m > 0) ? m : 1;
    if (fact_upper == 'N' && m > 0) {
        iwork_ptr = (int*)malloc((size_t)iwork_size_val * sizeof(int));
        CHECK_ALLOC(iwork_ptr);
    } else {
        iwork_ptr = (int*)malloc(1 * sizeof(int)); // Dummy allocation
        CHECK_ALLOC(iwork_ptr);
        iwork_ptr[0] = 0; // Initialize
    }

    // DWORK - Workspace Query
    double dwork_query_val[1]; // Changed name to avoid conflict
    int query_ldwork = -1;
    int info_query; // Use a separate info for the query

    // Define Fortran-style dummy LDs for the query call
    int lda_q = MAX(1, n);
    int ldb_q = MAX(1, n);
    int ldq_q = MAX(1, n);
    int ldr_q = MAX(1, m);
    int ldl_q = MAX(1, n);
    int ldg_q = MAX(1, n);
    int temp_oufact_q; // Dummy for oufact output in query

    F77_FUNC(sb02mt, SB02MT)(&jobg_upper, &jobl_upper, &fact_upper, &uplo_upper,
                             &n, &m,
                             NULL, &lda_q, NULL, &ldb_q, NULL, &ldq_q,
                             NULL, &ldr_q, NULL, &ldl_q,
                             NULL /*ipiv*/, &temp_oufact_q /*oufact*/, NULL, &ldg_q,
                             NULL /*iwork, Fortran will use DWORK for this if needed*/,
                             dwork_query_val, &query_ldwork, &info_query,
                             jobg_len, jobl_len, fact_len, uplo_len);

    if (info_query == 0) { // Query successful
        ldwork_val = (int)dwork_query_val[0];
    } else { // Query failed or not supported as expected, use documented minimums
        // info_query is not the main info, so don't assign it to main info here
        if (fact_upper == 'C' || (fact_upper == 'U' && jobg_upper == 'N' && jobl_upper == 'Z')) {
            ldwork_val = 1;
        } else if (fact_upper == 'N' && jobg_upper == 'N' && jobl_upper == 'Z') {
            ldwork_val = MAX(2, 3 * m);
        } else if (fact_upper == 'N' && (jobg_upper == 'G' || jobl_upper == 'N')) {
            ldwork_val = MAX(MAX(2, 3 * m), n * m); // N*M can be 0 if N or M is 0
             if (n * m == 0 && ldwork_val < MAX(2, 3*m)) ldwork_val = MAX(2,3*m); // ensure this part
             else if (n*m > 0) ldwork_val = MAX(ldwork_val, n*m);


        } else if (fact_upper == 'U' && (jobg_upper == 'G' || jobl_upper == 'N')) {
            ldwork_val = MAX(1, n * m);
        } else {
            ldwork_val = 1; // Default minimum
        }
    }
    // Ensure ldwork_val is at least 1, especially if n*m calculation results in 0.
    if (n == 0 && m == 0 && fact_upper == 'N' && jobg_upper == 'N' && jobl_upper == 'Z') {
         // Specific case from doc: LDWORK >= MAX(2,3*M) -> MAX(2,0) = 2
         ldwork_val = MAX(ldwork_val, 2);
    }
    ldwork_val = MAX(1, ldwork_val);


    dwork_ptr = (double*)malloc((size_t)ldwork_val * sizeof(double));
    CHECK_ALLOC(dwork_ptr);


    /* --- Prepare Arrays and Call Fortran Routine --- */
    size_t elem_size = sizeof(double);
    static double slicot_static_dummy_double = 0.0;
    static int slicot_static_dummy_int = 0;

    size_t a_elements = (jobl_upper == 'N' && n > 0) ? (size_t)n * n : 0;
    size_t b_elements = (n > 0 && m > 0) ? (size_t)n * m : 0;
    size_t q_elements = (jobl_upper == 'N' && n > 0) ? (size_t)n * n : 0;
    size_t r_elements = (m > 0) ? (size_t)m * m : 0;
    size_t l_elements = (jobl_upper == 'N' && n > 0 && m > 0) ? (size_t)n * m : 0;
    size_t g_elements = (jobg_upper == 'G' && n > 0) ? (size_t)n * n : 0;

    if (row_major) {
        lda_f = MAX(1, n); ldb_f = MAX(1, n); ldq_f = MAX(1, n);
        ldr_f = MAX(1, m); ldl_f = MAX(1, n); ldg_f = MAX(1, n);

        if (a_elements > 0) { a_cm = (double*)malloc(a_elements * elem_size); CHECK_ALLOC(a_cm); slicot_transpose_to_fortran_with_ld(a_io, a_cm, n, n, lda, lda_f, elem_size); a_f_ptr = a_cm; }
        else { a_f_ptr = (jobl_upper == 'N' ? a_io : &slicot_static_dummy_double); }

        if (b_elements > 0) { b_cm = (double*)malloc(b_elements * elem_size); CHECK_ALLOC(b_cm); slicot_transpose_to_fortran_with_ld(b_io, b_cm, n, m, ldb, ldb_f, elem_size); b_f_ptr = b_cm; }
        else { b_f_ptr = (n > 0 && m > 0 ? b_io : &slicot_static_dummy_double); }

        if (q_elements > 0) { q_cm = (double*)malloc(q_elements * elem_size); CHECK_ALLOC(q_cm); slicot_transpose_to_fortran_with_ld(q_io, q_cm, n, n, ldq, ldq_f, elem_size); q_f_ptr = q_cm; }
        else { q_f_ptr = (jobl_upper == 'N' ? q_io : &slicot_static_dummy_double); }

        if (r_elements > 0) { r_cm = (double*)malloc(r_elements * elem_size); CHECK_ALLOC(r_cm); slicot_transpose_to_fortran_with_ld(r_io, r_cm, m, m, ldr, ldr_f, elem_size); r_f_ptr = r_cm; }
        else { r_f_ptr = (m > 0 ? r_io : &slicot_static_dummy_double); }

        if (l_elements > 0) { l_cm = (double*)malloc(l_elements * elem_size); CHECK_ALLOC(l_cm); slicot_transpose_to_fortran_with_ld(l_io, l_cm, n, m, ldl, ldl_f, elem_size); l_f_ptr = l_cm; }
        else { l_f_ptr = (jobl_upper == 'N' && n > 0 && m > 0 ? l_io : &slicot_static_dummy_double); }

        if (g_elements > 0) { g_cm = (double*)malloc(g_elements * elem_size); CHECK_ALLOC(g_cm); g_f_ptr = g_cm; }
        else { g_f_ptr = (jobg_upper == 'G' ? g_out : &slicot_static_dummy_double); }

    } else { /* Column-Major Case */
        lda_f = lda; ldb_f = ldb; ldq_f = ldq; ldr_f = ldr; ldl_f = ldl; ldg_f = ldg;

        a_f_ptr = (jobl_upper == 'N' && n > 0) ? a_io : &slicot_static_dummy_double;
        b_f_ptr = (n > 0 && m > 0) ? b_io : &slicot_static_dummy_double;
        q_f_ptr = (jobl_upper == 'N' && n > 0) ? q_io : &slicot_static_dummy_double;
        r_f_ptr = (m > 0) ? r_io : &slicot_static_dummy_double;
        l_f_ptr = (jobl_upper == 'N' && n > 0 && m > 0) ? l_io : &slicot_static_dummy_double;
        g_f_ptr = (jobg_upper == 'G' && n > 0) ? g_out : &slicot_static_dummy_double;
    }

    ipiv_f_ptr = (fact_upper == 'U' && m > 0) ? ipiv_io : &slicot_static_dummy_int;
    if (fact_upper == 'N' && m == 0 && ipiv_io != NULL) { // If M=0 but IPIV was from test
        ipiv_f_ptr = &slicot_static_dummy_int; // Use dummy if M=0 for FACT='N'
    }


    /* Ensure Fortran LDs are at least 1 for the actual call */
    lda_f = MAX(1, lda_f); ldb_f = MAX(1, ldb_f); ldq_f = MAX(1, ldq_f);
    ldr_f = MAX(1, ldr_f); ldl_f = MAX(1, ldl_f); ldg_f = MAX(1, ldg_f);

    if (jobl_upper == 'Z') {
        if(n == 0 && a_f_ptr == &slicot_static_dummy_double) lda_f = 1;
        if(n == 0 && q_f_ptr == &slicot_static_dummy_double) ldq_f = 1;
        if(m == 0 && l_f_ptr == &slicot_static_dummy_double && n == 0) ldl_f = 1;
    }
     if (jobg_upper == 'N') {
        if(n == 0 && g_f_ptr == &slicot_static_dummy_double) ldg_f = 1;
    }
    if (m == 0 && fact_upper == 'U' && ipiv_f_ptr == &slicot_static_dummy_int) {
        // If M=0, IPIV not referenced even if FACT='U'
    }


    /* --- Call Fortran Routine --- */
    F77_FUNC(sb02mt, SB02MT)(&jobg_upper, &jobl_upper, &fact_upper, &uplo_upper,
                             &n, &m,
                             a_f_ptr, &lda_f, b_f_ptr, &ldb_f, q_f_ptr, &ldq_f,
                             r_f_ptr, &ldr_f, l_f_ptr, &ldl_f,
                             ipiv_f_ptr, oufact_out, g_f_ptr, &ldg_f,
                             iwork_ptr, dwork_ptr, &ldwork_val, &info, // Main info
                             jobg_len, jobl_len, fact_len, uplo_len);

    /* --- Copy Results Back to Row-Major Arrays if Needed --- */
    if (row_major && info == 0) { // Check main info
        if (jobl_upper == 'N' && a_cm && a_io && a_elements > 0) slicot_transpose_to_c_with_ld(a_cm, a_io, n, n, lda_f, lda, elem_size);
        // According to SB02MT.html, B is modified if oufact_out=1.
        if (*oufact_out == 1 && b_cm && b_io && b_elements > 0) slicot_transpose_to_c_with_ld(b_cm, b_io, n, m, ldb_f, ldb, elem_size);
        if (jobl_upper == 'N' && q_cm && q_io && q_elements > 0) slicot_transpose_to_c_with_ld(q_cm, q_io, n, n, ldq_f, ldq, elem_size);
        // R is modified if oufact_out=1 or 2
        if ((*oufact_out == 1 || *oufact_out == 2) && r_cm && r_io && r_elements > 0) slicot_transpose_to_c_with_ld(r_cm, r_io, m, m, ldr_f, ldr, elem_size);
        // L is modified if jobl='N' and oufact_out=1
        if (jobl_upper == 'N' && *oufact_out == 1 && l_cm && l_io && l_elements > 0) slicot_transpose_to_c_with_ld(l_cm, l_io, n, m, ldl_f, ldl, elem_size);

        if (jobg_upper == 'G' && g_cm && g_out && g_elements > 0) slicot_transpose_to_c_with_ld(g_cm, g_out, n, n, ldg_f, ldg, elem_size);
    }
    // IPIV is directly modified in ipiv_io if fact_upper='U' and oufact_out=2, no transpose needed for 1D int array.

cleanup:
    free(dwork_ptr);
    free(iwork_ptr);
    if (row_major) {
        free(a_cm); free(b_cm); free(q_cm);
        free(r_cm); free(l_cm); free(g_cm);
    }
    return info; // Return main info
 }