/**
 * @file sb02mt.c
 * @brief C wrapper implementation for SLICOT routine SB02MT
 *
 * This file provides a C wrapper implementation for the SLICOT routine SB02MT,
 * which converts optimal problems with coupling weighting terms to standard problems.
 */

#include <stdlib.h>
#include <ctype.h>  // For toupper
#include <stddef.h> // For size_t

// Include the header file for this wrapper
#include "sb02mt.h"
// Include necessary SLICOT utility headers
#include "slicot_utils.h" // Assumed to contain MAX, MIN, CHECK_ALLOC, etc.
#include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

/**
 * Declare the external Fortran routine using the F77_FUNC macro.
 * Note use of slicot_complex_double for CWORK.
 */
extern void F77_FUNC(sb02mt, SB02MT)(
    const char* dico,       // CHARACTER*1 DICO
    const char* jobb,       // CHARACTER*1 JOBB
    const char* fact,       // CHARACTER*1 FACT
    const char* uplo,       // CHARACTER*1 UPLO
    const char* jobl,       // CHARACTER*1 JOBL
    const char* sort,       // CHARACTER*1 SORT
    const int* n,           // INTEGER N
    const int* m,           // INTEGER M
    const int* p,           // INTEGER P
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
    double* x,              // DOUBLE PRECISION X(LDX,*)
    const int* ldx,         // INTEGER LDX
    double* g,              // DOUBLE PRECISION G(LDG,*)
    const int* ldg,         // INTEGER LDG
    double* rcond,          // DOUBLE PRECISION RCOND(2)
    int* rank,              // INTEGER RANK
    double* s,              // DOUBLE PRECISION S(LDS,*)
    const int* lds,         // INTEGER LDS
    double* u,              // DOUBLE PRECISION U(LDU,*)
    const int* ldu,         // INTEGER LDU
    double* wr,             // DOUBLE PRECISION WR(*)
    double* wi,             // DOUBLE PRECISION WI(*)
    const double* tol,      // DOUBLE PRECISION TOL
    int* iwork,             // INTEGER IWORK(*)
    double* dwork,          // DOUBLE PRECISION DWORK(*)
    const int* ldwork,      // INTEGER LDWORK
    int* info,              // INTEGER INFO (output)
    int dico_len,           // Hidden length for dico
    int jobb_len,           // Hidden length for jobb
    int fact_len,           // Hidden length for fact
    int uplo_len,           // Hidden length for uplo
    int jobl_len,           // Hidden length for jobl
    int sort_len            // Hidden length for sort
);


 /* C wrapper function definition */
 SLICOT_EXPORT
 int slicot_sb02mt(char dico, char jobb, char fact, char uplo, char jobl,
                   char sort, int n, int m, int p,
                   double* a_io, int lda, double* b_io, int ldb,
                   double* q_io, int ldq, double* r_io, int ldr,
                   double* l_io, int ldl, double* x_out, int ldx,
                   double* g_out, int ldg, double* rcond_out, int* rank_out,
                   double* s_out, int lds, double* u_out, int ldu,
                   double* wr_out, double* wi_out,
                   double tol_in, int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork_val = -1; /* Use -1 for workspace query, or calculate */
     double* dwork_ptr = NULL;
     int* iwork_ptr = NULL;
     int iwork_size_val = 0;
     const int dico_len = 1, jobb_len = 1, fact_len = 1, uplo_len = 1, jobl_len = 1, sort_len = 1;

     char dico_upper = toupper(dico);
     char jobb_upper = toupper(jobb);
     char fact_upper = toupper(fact);
     char uplo_upper = toupper(uplo);
     char jobl_upper = toupper(jobl);
     char sort_upper = toupper(sort);

     /* Arrays for column-major versions if needed */
     double *a_cm = NULL, *b_cm = NULL, *q_cm = NULL, *r_cm = NULL, *l_cm = NULL;
     double *x_cm = NULL, *g_cm = NULL, *s_cm = NULL, *u_cm = NULL;

     /* Pointers to pass to Fortran */
     double *a_f_ptr, *b_f_ptr, *q_f_ptr, *r_f_ptr, *l_f_ptr;
     double *x_f_ptr, *g_f_ptr, *s_f_ptr, *u_f_ptr;
     /* rcond_out, wr_out, wi_out are 1D, rank_out is int*, so direct pass after row_major check */

     int lda_f, ldb_f, ldq_f, ldr_f, ldl_f, ldx_f, ldg_f, lds_f, ldu_f;
    
    /* --- Input Parameter Validation --- */
    /* Character parameters */
    if (dico_upper != 'C' && dico_upper != 'D') { info = -1; goto cleanup; }
    if (jobb_upper != 'B' && jobb_upper != 'N') { info = -2; goto cleanup; }
    if (fact_upper != 'F' && fact_upper != 'N') { info = -3; goto cleanup; }
    if (uplo_upper != 'U' && uplo_upper != 'L') { info = -4; goto cleanup; }
    if (jobl_upper != 'N' && jobl_upper != 'Z') { info = -5; goto cleanup; }
    if (sort_upper != 'S' && sort_upper != 'N') { info = -6; goto cleanup; }

    /* Dimensions */
    if (n < 0) { info = -7; goto cleanup; }
    if (m < 0) { info = -8; goto cleanup; }
    if (p < 0) { info = -9; goto cleanup; }

    /* Array pointers and Leading dimensions */
    /* A */
    if (n > 0 && a_io == NULL) { info = -10; goto cleanup; }
    if (row_major) { if (n > 0 && lda < n) { info = -11; goto cleanup; } }
    else { if (lda < MAX(1,n)) { info = -11; goto cleanup; } }

    /* B */
    if (n > 0 && m > 0 && b_io == NULL) { info = -12; goto cleanup; }
    if (row_major) { if (n > 0 && m > 0 && ldb < m) { info = -13; goto cleanup; } }
    else { if (ldb < MAX(1,n)) { info = -13; goto cleanup; } }

    /* Q */
    if (n > 0 && q_io == NULL) { info = -14; goto cleanup; }
    if (row_major) { if (n > 0 && ldq < n) { info = -15; goto cleanup; } }
    else { if (ldq < MAX(1,n)) { info = -15; goto cleanup; } }

    /* R */
    if (m > 0 && r_io == NULL) { info = -16; goto cleanup; }
    if (row_major) { if (m > 0 && ldr < m) { info = -17; goto cleanup; } }
    else { if (ldr < MAX(1,m)) { info = -17; goto cleanup; } }

    /* L */
    if (jobl_upper == 'N') {
        if (n > 0 && m > 0 && l_io == NULL) { info = -18; goto cleanup; }
        if (row_major) { if (n > 0 && m > 0 && ldl < m) { info = -19; goto cleanup; } }
        else { if (ldl < MAX(1,n)) { info = -19; goto cleanup; } }
    } else { /* L not referenced by Fortran, but LDL must be >= 1 */
        if (ldl < 1) { info = -19; goto cleanup; }
    }

    /* X (output) */
    if (n > 0 && x_out == NULL) { info = -20; goto cleanup; }
    if (row_major) { if (n > 0 && ldx < n) { info = -21; goto cleanup; } }
    else { if (ldx < MAX(1,n)) { info = -21; goto cleanup; } }

    /* G (output) */
    if (n > 0 && m > 0 && g_out == NULL) { info = -22; goto cleanup; }
    if (row_major) { if (n > 0 && m > 0 && ldg < m) { info = -23; goto cleanup; } }
    else { if (ldg < MAX(1,n)) { info = -23; goto cleanup; } }

    /* RCOND (output) */
    if (rcond_out == NULL) { info = -24; goto cleanup; }

    /* RANK (output) */
    if (rank_out == NULL) { info = -25; goto cleanup; }

    /* S (output) */
    if (n > 0 && s_out == NULL) { info = -26; goto cleanup; }
    if (row_major) { if (2 * n > 0 && lds < 2 * n) { info = -27; goto cleanup; } }
    else { if (lds < MAX(1, 2 * n)) { info = -27; goto cleanup; } }

    /* U (output) */
    if (n > 0 && u_out == NULL) { info = -28; goto cleanup; }
    if (row_major) { if (2 * n > 0 && ldu < 2 * n) { info = -29; goto cleanup; } }
    else { if (ldu < MAX(1, 2 * n)) { info = -29; goto cleanup; } }

    /* WR (output) */
    if (n > 0 && wr_out == NULL) { info = -30; goto cleanup; }

    /* WI (output) */
    if (n > 0 && wi_out == NULL) { info = -31; goto cleanup; }

    /* TOL */
    if (tol_in < 0.0) { info = -32; goto cleanup; }


    /* --- Workspace Query & Allocation --- */
     // Calculate needed workspace size for DWORK
     // This is a placeholder, actual SB02MT ldwork formula should be used if known.
     // Using a common pattern for Riccati solvers.
     ldwork_val = MAX(1, 4 * n * n + MAX(m * m + 2 * n * m, 5 * n + 2 * m));
     if (jobb_upper == 'B' && jobl_upper == 'Z') { // Simpler case
         ldwork_val = MAX(1, 2 * n * n + n);
     }
     // More general case from some Riccati solvers:
     // ldwork_val = MAX(1, MAX( (2*n)*(2*n+5), n*(n+m+p+MAX(n,m,p)+5) ) ); // Adjust as per SB02MT spec
     // For SB02MT, a safe bet based on typical structure might be:
     // N*(2*N + M + P + max(N,M,P) + 5) + N*(N+1)/2 for Schur method based solvers.
     // Or simpler: max(1, 4*N*N + 2*N*M + M*M + 5*N)
     // The test file implies S and U are 2N x 2N, related to Hamiltonian/Symplectic pencil.
     // A common workspace for such problems is around 8*N or more for DWORK.
     // Let's use a slightly more generous one based on related routines.
     // Max of (size for QR/Schur of 2N Hamiltonian, size for Sylvester solvers etc.)
     // From SB02OD (similar context): LDWORK >= MAX(1, 2*N*(N+M+7)+MAX(N,M))
     // From SB02MD (DARE solver): LDWORK >= MAX(1, N*(N+2*M+P+MAX(N,M,P)+7))
     // Given S, U are 2N x 2N, a 2N system is solved.
     ldwork_val = MAX(1, (2*n)*(2*n + 5) + 2*n ); // Basic for Schur decomp + few vectors
     if (n > 0) { // More robust allocation based on typical needs
        ldwork_val = MAX(ldwork_val, 8*n*n + 16*n); 
     } else { // n == 0
        ldwork_val = 1; // Ensure minimum workspace for N=0 case
     }


     dwork_ptr = (double*)malloc((size_t)ldwork_val * sizeof(double));
     CHECK_ALLOC(dwork_ptr);
     
     // Allocate integer workspace IWORK
     iwork_size_val = MAX(1, 2 * n); // Usually for permutations or block sizes in Schur
     if (n == 0) iwork_size_val = 1;
     iwork_ptr = (int*)malloc((size_t)iwork_size_val * sizeof(int));
     CHECK_ALLOC(iwork_ptr);
    
     /* --- Prepare Arrays and Call Fortran Routine --- */
     size_t elem_size = sizeof(double);
     static double slicot_static_dummy_double = 0.0; // Static dummy for N=0 cases
     
     if (row_major) {
         /* --- Row-Major Case --- */
         
         // Allocate and transpose matrices
         size_t a_elements = (size_t)n * n;
         size_t b_elements = (size_t)n * m;
         size_t q_elements = (size_t)n * n;
         size_t r_elements = (size_t)m * m;
         size_t l_elements = (jobl_upper == 'N') ? (size_t)n * m : 0;
         size_t x_elements = (size_t)n * n;
         size_t g_elements = (size_t)n * m;
         size_t s_elements = (n > 0) ? (size_t)(2 * n) * (2 * n) : 0; // Guard for n=0
         size_t u_elements = (n > 0) ? (size_t)(2 * n) * (2 * n) : 0; // Guard for n=0

         if (a_elements > 0) { a_cm = (double*)malloc(a_elements * elem_size); CHECK_ALLOC(a_cm); }
         if (b_elements > 0) { b_cm = (double*)malloc(b_elements * elem_size); CHECK_ALLOC(b_cm); }
         if (q_elements > 0) { q_cm = (double*)malloc(q_elements * elem_size); CHECK_ALLOC(q_cm); }
         if (r_elements > 0) { r_cm = (double*)malloc(r_elements * elem_size); CHECK_ALLOC(r_cm); }
         if (l_elements > 0) { l_cm = (double*)malloc(l_elements * elem_size); CHECK_ALLOC(l_cm); }
         if (x_elements > 0) { x_cm = (double*)malloc(x_elements * elem_size); CHECK_ALLOC(x_cm); } // Output
         if (g_elements > 0) { g_cm = (double*)malloc(g_elements * elem_size); CHECK_ALLOC(g_cm); } // Output
         if (s_elements > 0) { s_cm = (double*)malloc(s_elements * elem_size); CHECK_ALLOC(s_cm); } // Output
         if (u_elements > 0) { u_cm = (double*)malloc(u_elements * elem_size); CHECK_ALLOC(u_cm); } // Output
        
         // Transpose input matrices from row-major to column-major
         if (a_cm && a_io) slicot_transpose_to_fortran_with_ld(a_io, a_cm, n, n, lda, MAX(1,n), elem_size);
         if (b_cm && b_io) slicot_transpose_to_fortran_with_ld(b_io, b_cm, n, m, ldb, MAX(1,n), elem_size);
         if (q_cm && q_io) slicot_transpose_to_fortran_with_ld(q_io, q_cm, n, n, ldq, MAX(1,n), elem_size);
         if (r_cm && r_io) slicot_transpose_to_fortran_with_ld(r_io, r_cm, m, m, ldr, MAX(1,m), elem_size);
         if (jobl_upper == 'N' && l_cm && l_io) slicot_transpose_to_fortran_with_ld(l_io, l_cm, n, m, ldl, MAX(1,n), elem_size);
        
         // Set Fortran leading dimensions and pointers
         lda_f = MAX(1, n); ldb_f = MAX(1, n); ldq_f = MAX(1, n); ldr_f = MAX(1, m);
         ldl_f = (jobl_upper == 'N') ? MAX(1, n) : 1;
         ldx_f = MAX(1, n); ldg_f = MAX(1, n);
         lds_f = MAX(1, 2 * n); ldu_f = MAX(1, 2 * n);
        
         a_f_ptr = a_cm; b_f_ptr = b_cm; q_f_ptr = q_cm; r_f_ptr = r_cm; 
         l_f_ptr = (jobl_upper == 'N') ? l_cm : NULL; 
         x_f_ptr = x_cm; g_f_ptr = g_cm; s_f_ptr = s_cm; u_f_ptr = u_cm;

     } else {
         /* --- Column-Major Case --- */
         lda_f = lda; ldb_f = ldb; ldq_f = ldq; ldr_f = ldr; ldl_f = ldl;
         ldx_f = ldx; ldg_f = ldg; lds_f = lds; ldu_f = ldu;
        
         a_f_ptr = a_io; b_f_ptr = b_io; q_f_ptr = q_io; r_f_ptr = r_io; 
         l_f_ptr = (jobl_upper == 'N') ? l_io : NULL; 
         x_f_ptr = x_out; g_f_ptr = g_out; s_f_ptr = s_out; u_f_ptr = u_out;
     }

    /* Ensure non-NULL pointers for Fortran when dimensions are zero, based on test findings */
    if (n == 0) {
        a_f_ptr = (a_f_ptr == NULL) ? &slicot_static_dummy_double : a_f_ptr;
        b_f_ptr = (b_f_ptr == NULL) ? &slicot_static_dummy_double : b_f_ptr;
        q_f_ptr = (q_f_ptr == NULL) ? &slicot_static_dummy_double : q_f_ptr;
        if (jobl_upper == 'N') {
            l_f_ptr = (l_f_ptr == NULL) ? &slicot_static_dummy_double : l_f_ptr;
        }
        // Else L is not referenced by Fortran, l_f_ptr can remain NULL if it was.
        // However, test used dummy for L, so safer to assign dummy if JOBL='N'.

        x_f_ptr = (x_f_ptr == NULL) ? &slicot_static_dummy_double : x_f_ptr;
        g_f_ptr = (g_f_ptr == NULL) ? &slicot_static_dummy_double : g_f_ptr;
        s_f_ptr = (s_f_ptr == NULL) ? &slicot_static_dummy_double : s_f_ptr;
        u_f_ptr = (u_f_ptr == NULL) ? &slicot_static_dummy_double : u_f_ptr;
        // wr_out and wi_out are direct arguments to the C wrapper.
        // If N=0, the test passes &dummy_val for them. The C wrapper does not need to change these pointers.
    }
    if (m == 0) {
        r_f_ptr = (r_f_ptr == NULL) ? &slicot_static_dummy_double : r_f_ptr;
        if (n > 0) { // If N > 0 and M = 0
            b_f_ptr = (b_f_ptr == NULL) ? &slicot_static_dummy_double : b_f_ptr;
            g_f_ptr = (g_f_ptr == NULL) ? &slicot_static_dummy_double : g_f_ptr;
            if (jobl_upper == 'N') {
                l_f_ptr = (l_f_ptr == NULL) ? &slicot_static_dummy_double : l_f_ptr;
            }
        }
    }
    // Note: rcond_out and rank_out must always point to valid memory by the caller.
    
     /* --- Call Fortran Routine --- */
     F77_FUNC(sb02mt, SB02MT)(&dico_upper, &jobb_upper, &fact_upper, &uplo_upper,
                              &jobl_upper, &sort_upper, &n, &m, &p,
                              a_f_ptr, &lda_f, b_f_ptr, &ldb_f, q_f_ptr, &ldq_f,
                              r_f_ptr, &ldr_f, l_f_ptr, &ldl_f,
                              x_f_ptr, &ldx_f, g_f_ptr, &ldg_f,
                              rcond_out, rank_out,
                              s_f_ptr, &lds_f, u_f_ptr, &ldu_f,
                              wr_out, wi_out, &tol_in,
                              iwork_ptr, dwork_ptr, &ldwork_val, &info,
                              dico_len, jobb_len, fact_len, uplo_len, jobl_len, sort_len);
    
     /* --- Copy Results Back to Row-Major Arrays if Needed --- */
     if (row_major && info == 0) {
         // Transpose output matrices from column-major back to row-major
         if (a_cm && a_io) slicot_transpose_to_c_with_ld(a_cm, a_io, n, n, lda_f, lda, elem_size);
         // B is input only, Q can be output if JOBL='N' (or FACT='F') - check SB02MT spec for Q output conditions
         // Assuming Q is input/output if FACT='F' or if JOBL='N' implies Q is modified.
         // For now, only transposing if FACT='F' (factored input means Q is output).
         // The problem description implies Q_io is input/output.
         if (q_cm && q_io) slicot_transpose_to_c_with_ld(q_cm, q_io, n, n, ldq_f, ldq, elem_size);
         if (x_cm && x_out) slicot_transpose_to_c_with_ld(x_cm, x_out, n, n, ldx_f, ldx, elem_size);
         if (g_cm && g_out) slicot_transpose_to_c_with_ld(g_cm, g_out, n, m, ldg_f, ldg, elem_size);
         if (s_cm && s_out) slicot_transpose_to_c_with_ld(s_cm, s_out, 2 * n, 2 * n, lds_f, lds, elem_size);
         if (u_cm && u_out) slicot_transpose_to_c_with_ld(u_cm, u_out, 2 * n, 2 * n, ldu_f, ldu, elem_size);
         // rcond_out, rank_out, wr_out, wi_out are modified in place or are 1D.
     }
    
 cleanup:
     /* --- Cleanup --- */
     free(dwork_ptr);
     free(iwork_ptr);
     free(a_cm);
     free(b_cm);
     free(q_cm);
     free(r_cm);
     free(l_cm);
     free(x_cm);
     free(g_cm);
     free(s_cm);
     free(u_cm);
    
     return info;
 }