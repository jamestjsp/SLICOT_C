/**
 * @file ib03ad.c
 * @brief C wrapper for SLICOT routine IB03AD.
 * @details Estimates parameters of a Wiener system using Levenberg-Marquardt.
 * Workspace (IWORK, DWORK) is allocated internally by this wrapper.
 * Input/output matrix format for U and Y is handled via the row_major parameter.
 * C-level validation for zero dimensions aims to align with Fortran routine behavior.
 */

#include <stdlib.h> // For malloc, free
#include <ctype.h>  // For toupper
#include <stddef.h> // For size_t
#include <math.h>   // For MAX/MIN if needed (often provided by slicot_utils.h)
#include <stdio.h>  // For error logging (optional)

#include "ib03ad.h"       // Public header for this wrapper
#include "slicot_utils.h" // Provides CHECK_ALLOC, SLICOT_MEMORY_ERROR, MAX/MIN, transpose functions etc.
#include "slicot_f77.h"   // Provides F77_FUNC macro for Fortran name mangling

/* External Fortran routine declaration */
extern void F77_FUNC(ib03ad, IB03AD)(
    const char* init, const char* alg, const char* stor,
    const int* nobr, const int* m, const int* l, const int* nsmp,
    int* n, const int* nn,
    const int* itmax1, const int* itmax2, const int* nprint,
    const double* u, const int* ldu,
    const double* y, const int* ldy,
    double* x, int* lx,
    const double* tol1, const double* tol2,
    int* iwork, double* dwork, const int* ldwork,
    int* iwarn, int* info,
    int init_len, int alg_len, int stor_len);

/* C wrapper function definition */
SLICOT_EXPORT
int slicot_ib03ad(
    char init_char, char alg_char, char stor_char,
    int nobr, int m, int l, int nsmp,
    int* n_ptr, /* N is in/out */
    int nn,
    int itmax1, int itmax2, int nprint,
    const double* u, int ldu,
    const double* y, int ldy,
    double* x, int* lx_ptr, /* LX is in/out */
    double tol1, double tol2,
    int* out_iwork_summary, /* Optional: IWORK(1-3) */
    double* out_dwork_summary, /* Optional: DWORK(1-10) or DWORK(1-(10+IWORK(3))) */
    int* iwarn_ptr, /* Output warning indicator */
    int row_major)
{
    // 1. Variable declarations
    int info = 0;
    int local_iwarn = 0;
    char init_upper, alg_upper, stor_upper;

    int *iwork_alloc = NULL;
    double *dwork_alloc = NULL;
    int liwork_calc = 0;
    int ldwork_calc = 0;

    double* u_cm = NULL;
    double* y_cm = NULL;

    const double* u_f_ptr = u;
    const double* y_f_ptr = y;

    int ldu_f, ldy_f;
    size_t u_size = 0, y_size = 0;

    const int init_len = 1;
    const int alg_len = 1;
    const int stor_len = 1;

    int n_val = *n_ptr; // Dereference for local use and validation
    int lx_val = *lx_ptr;
    const int original_n = n_val;
    const int original_lx = lx_val;


    // 2. Input parameter validation
    init_upper = toupper(init_char);
    alg_upper = toupper(alg_char);
    stor_upper = toupper(stor_char);

    if (init_upper != 'L' && init_upper != 'S' && init_upper != 'B' && init_upper != 'N') { info = -1; goto cleanup; }
    if (alg_upper != 'D' && alg_upper != 'I') { info = -2; goto cleanup; }
    if (alg_upper == 'D' && stor_upper != 'F' && stor_upper != 'P') { info = -3; goto cleanup; }

    if (init_upper == 'L' || init_upper == 'B') {
        if (nobr <= 0) { info = -4; goto cleanup; }
    }
    if (m < 0) { info = -5; goto cleanup; }
    if (l < 0) { info = -6; goto cleanup; } // L >= 0
    if (init_upper == 'L' || init_upper == 'B') {
        if (l <= 0) { info = -6; goto cleanup; } // L > 0 if INIT = L/B
    }

    if (nsmp < 0) { info = -7; goto cleanup; }
    if (init_upper == 'L' || init_upper == 'B') {
        if (nsmp < 2 * (m + l + 1) * nobr - 1) { info = -7; goto cleanup; }
    }

    // N validation
    if (init_upper == 'S' || init_upper == 'N') {
        if (n_val < 0) { info = -8; goto cleanup; }
    } else { // INIT = 'L' or 'B'
        if (n_val >= 0) { // N is provided
            if (n_val >= nobr || n_val == 0) { info = -8; goto cleanup; } // N=0 not acceptable if INIT=L/B
        }
        // If n_val < 0, it's an output, no check here.
    }

    if (nn < 0) { info = -9; goto cleanup; }

    if (init_upper == 'S' || init_upper == 'B') {
        if (itmax1 < 0) { info = -10; goto cleanup; }
    }
    if (itmax2 < 0) { info = -11; goto cleanup; }
    // NPRINT can be any integer.

    // U matrix and LDU
    if (m > 0) { // U is referenced if M > 0
        if (u == NULL && nsmp > 0) { info = -13; goto cleanup; }
        if (u != NULL || nsmp > 0) { // Check LDU if U is provided or should be
             if (row_major) { // C LDU is cols (M)
                if (ldu < m) { info = -14; goto cleanup; }
            } else { // Fortran LDU is rows (NSMP)
                if (ldu < MAX(1,nsmp)) { info = -14; goto cleanup; }
            }
        }
    } else { // M == 0
        if (u != NULL && ldu < 1) { info = -14; goto cleanup; } // LDU must be >=1 if u not NULL
    }


    // Y matrix and LDY
    if (l > 0) { // Y is referenced if L > 0
        if (y == NULL && nsmp > 0) { info = -15; goto cleanup; }
         if (y != NULL || nsmp > 0) { // Check LDY if Y is provided or should be
            if (row_major) { // C LDY is cols (L)
                if (ldy < l) { info = -16; goto cleanup; }
            } else { // Fortran LDY is rows (NSMP)
                if (ldy < MAX(1,nsmp)) { info = -16; goto cleanup; }
            }
        }
    } else { // L == 0 (only possible if INIT = 'S' or 'N')
         if (y != NULL && ldy < 1) { info = -16; goto cleanup; }
    }


    // X array and LX
    // Calculate required NX
    int n_for_lx_calc = (n_val < 0 && (init_upper == 'L' || init_upper == 'B')) ? MAX(0, nobr - 1) : MAX(0, n_val);
    int bsn_lx = nn * (MAX(0,l) + 2) + 1; // Use MAX(0,l) in case L=0 for INIT=S/N
    int lths_lx = n_for_lx_calc * (MAX(0,l) + MAX(0,m) + 1) + MAX(0,l) * MAX(0,m);
    int required_nx = bsn_lx * MAX(0,l) + lths_lx;
    if (required_nx < 0) required_nx = 0; // Should not happen with MAX(0,...)

    if (x == NULL && required_nx > 0) { info = -17; goto cleanup; } // X required if params exist
    if (n_val >= 0) { // If N is known, LX must be sufficient
        if (lx_val < required_nx) { info = -18; goto cleanup; }
    }
    // If N < 0 (estimated), LX check is tricky. Fortran routine handles if LX is too small.
    // We ensure LX is at least 1 if X is not NULL.
    if (x != NULL && lx_val < 1 && required_nx > 0) { info = -18; goto cleanup; }
    if (required_nx == 0 && x != NULL && lx_val < 1) { info = -18; goto cleanup; }


    // TOL1, TOL2 can be any double.

    if (info != 0) { goto cleanup; }

    // 3. Internal Workspace Allocation

    // IWORK
    int liw1 = 0, liw2 = 0;
    if (init_upper == 'S' || init_upper == 'N') {
        liw1 = 0;
        liw2 = 0;
    } else { // INIT = 'L' or 'B'
        liw1 = m + l;
        liw2 = MAX(m * nobr + n_for_lx_calc, m * (n_for_lx_calc + l));
    }
    liwork_calc = MAX(3, MAX(liw1, liw2));
    if (liwork_calc > 0) {
        iwork_alloc = (int*)malloc((size_t)liwork_calc * sizeof(int));
        CHECK_ALLOC(iwork_alloc);
    } else {
        iwork_alloc = NULL; // Should not happen due to MAX(3,...)
    }

    // DWORK
    // N for LDWORK formulas: "N should be taken not larger than NOBR - 1, if N < 0 on entry."
    int n_eff_ldwork = (n_val < 0 && (init_upper == 'L' || init_upper == 'B')) ? MAX(0, nobr - 1) : MAX(0, n_val);

    long long m_ll = m;
    long long l_ll = l;
    long long nsmp_ll = nsmp;
    long long nobr_ll = nobr;
    long long n_ll_eff = n_eff_ldwork; // Effective N for workspace calculation
    long long nn_ll = nn;

    long long bsn  = nn_ll * (l_ll + 2) + 1;
    long long lths = n_ll_eff * (l_ll + m_ll + 1) + l_ll * m_ll;
    long long nx_calc = bsn * l_ll + lths;
    if (nx_calc < 0) nx_calc = 0;


    long long ldw1 = 0, ldw2 = 0, ldw3_val = 0, ldw4_val = 0; // Use _val to avoid conflict
    long long ldw5 = 0, ldw6 = 0, ldw7 = 0, ldw8 = 0;
    long long lw1 = 0, lw2_val = 0, lw3 = 0, lw4 = 0; // Use _val for lw2 to avoid conflict

    if (init_upper == 'S' || init_upper == 'N') {
        lw1 = 0;
    } else { // INIT = 'L' or 'B'
        ldw1 = MAX( 2LL*(l_ll*nobr_ll-l_ll)*n_ll_eff+2LL*n_ll_eff, (l_ll*nobr_ll-l_ll)*n_ll_eff+n_ll_eff*n_ll_eff+7LL*n_ll_eff );
        ldw1 = MAX( ldw1, l_ll*nobr_ll*n_ll_eff +
                        MAX( (l_ll*nobr_ll-l_ll)*n_ll_eff+2LL*n_ll_eff + (2LL*m_ll+l_ll)*nobr_ll+l_ll,
                             2LL*(l_ll*nobr_ll-l_ll)*n_ll_eff+n_ll_eff*n_ll_eff+8LL*n_ll_eff ) );
        ldw1 = MAX( ldw1, MAX( n_ll_eff+4LL*(m_ll*nobr_ll+n_ll_eff)+1LL, m_ll*nobr_ll+3LL*n_ll_eff+l_ll ) );

        if (m_ll == 0) ldw2 = 0;
        else ldw2 = l_ll*nobr_ll*n_ll_eff + m_ll*nobr_ll*(n_ll_eff+l_ll)*(m_ll*(n_ll_eff+l_ll)+1LL) +
                    MAX( (n_ll_eff+l_ll)*(n_ll_eff+l_ll), 4LL*m_ll*(n_ll_eff+l_ll)+1LL );

        ldw3_val = nsmp_ll*l_ll*(n_ll_eff+1LL) + 2LL*n_ll_eff + MAX( 2LL*n_ll_eff*n_ll_eff, 4LL*n_ll_eff );
        ldw4_val = n_ll_eff*(n_ll_eff+1LL) + 2LL*n_ll_eff +
                   MAX( n_ll_eff*l_ll*(n_ll_eff+1LL) + 2LL*n_ll_eff*n_ll_eff + l_ll*n_ll_eff, 4LL*n_ll_eff );
        ldw5 = nsmp_ll*l_ll + (n_ll_eff+l_ll)*(n_ll_eff+m_ll) + 3LL*n_ll_eff+m_ll+l_ll;
        ldw6 = nsmp_ll*l_ll + (n_ll_eff+l_ll)*(n_ll_eff+m_ll) + n_ll_eff +
               MAX(1LL, n_ll_eff*n_ll_eff*l_ll + n_ll_eff*l_ll + n_ll_eff ); // Simplified MAX part
        ldw6 = MAX(ldw6, nsmp_ll*l_ll + (n_ll_eff+l_ll)*(n_ll_eff+m_ll) + n_ll_eff + n_ll_eff*n_ll_eff +
                         MAX(n_ll_eff*n_ll_eff + n_ll_eff*MAX(n_ll_eff,l_ll) + 6LL*n_ll_eff + MIN(n_ll_eff,l_ll),
                             n_ll_eff*m_ll) );


        lw1 = MAX( 2LL*(m_ll+l_ll)*nobr_ll*(2LL*(m_ll+l_ll)*(nobr_ll+1LL)+3LL) + l_ll*nobr_ll,
                   4LL*(m_ll+l_ll)*nobr_ll*(m_ll+l_ll)*nobr_ll + (n_ll_eff+l_ll)*(n_ll_eff+m_ll) + MAX( ldw1, ldw2 ) );
        lw1 = MAX( lw1, (n_ll_eff+l_ll)*(n_ll_eff+m_ll) + n_ll_eff + n_ll_eff*n_ll_eff + 2LL + n_ll_eff*(n_ll_eff+m_ll+l_ll) +
                        MAX( MAX( 5LL*n_ll_eff, 2LL ), MAX( MIN( ldw3_val, ldw4_val ), MAX(ldw5, ldw6) ) ) );
    }

    if (init_upper == 'L' || init_upper == 'N') {
        lw2_val = 0;
        ldw8 = 0; // Not used if lw2_val is 0
    } else { // INIT = 'S' or 'B'
        if (alg_upper == 'D' && stor_upper == 'F') ldw7 = bsn*bsn;
        else if (alg_upper == 'D' && stor_upper == 'P') ldw7 = bsn*(bsn+1LL)/2LL;
        else /* ALG = 'I' */ ldw7 = 3LL*bsn + nsmp_ll; // Assuming NSMP*L in doc is NSMP for CG part
        lw2_val = nsmp_ll*l_ll + MAX( 5LL, nsmp_ll + 2LL*bsn + nsmp_ll*bsn + MAX( 2LL*nn_ll + bsn, ldw7 ) );

        if (m_ll > 0) ldw8 = nsmp_ll*l_ll + (n_ll_eff+l_ll)*(n_ll_eff+m_ll) + 3LL*n_ll_eff+m_ll+l_ll;
        else ldw8 = nsmp_ll*l_ll + (n_ll_eff+l_ll)*n_ll_eff + 2LL*n_ll_eff+l_ll;
    }
    // LW3 is always calculated, but its value might be overridden by LW2 if INIT=S/B
    // The formula for LW3 in doc seems to be for the case when INIT=L or N,
    // because it doesn't depend on BSN or NN.
    // The doc says "LW2 = LW3 = 0, if INIT = 'L' or 'N'; otherwise,"
    // This implies LW3 is specific to INIT = S or B.
    // Let's assume LW3 is for the "whole optimization" part if non-linear init is done.
    // If INIT = 'S' or 'B', LW3 is calculated. Otherwise it's 0.
    if (init_upper == 'S' || init_upper == 'B') {
         lw3 = MAX( ldw8, nsmp_ll*l_ll + (n_ll_eff+l_ll)*(2LL*n_ll_eff+m_ll) + 2LL*n_ll_eff );
    } else {
        lw3 = 0;
    }


    long long l0, l1_val, l2_val; // Use _val for l1, l2
    if (m_ll > 0) l0 = MAX( n_ll_eff*(n_ll_eff+l_ll), n_ll_eff+m_ll+l_ll );
    else l0 = MAX( n_ll_eff*(n_ll_eff+l_ll), l_ll );
    l1_val = nsmp_ll*l_ll + MAX( 2LL*nn_ll, (n_ll_eff+l_ll)*(n_ll_eff+m_ll) + 2LL*n_ll_eff + l0);

    if (alg_upper == 'D' && stor_upper == 'F') l2_val = nx_calc*nx_calc;
    else if (alg_upper == 'D' && stor_upper == 'P') l2_val = nx_calc*(nx_calc+1LL)/2LL;
    else /* ALG = 'I' */ l2_val = 3LL*nx_calc + nsmp_ll*l_ll;

    lw4 = MAX( 30LL, nsmp_ll*l_ll + 2LL*nx_calc + nsmp_ll*l_ll*( bsn + lths ) +
                  MAX( MAX( l1_val + nx_calc, nsmp_ll*l_ll + l1_val), l2_val ) );

    ldwork_calc = (int)MAX(lw1, MAX(lw2_val, MAX(lw3, lw4)));

    // Ensure minimum workspace when dimensions are zero
    // Fortran routine requires minimum workspace even for edge cases
    if (ldwork_calc < 1) {
        if (init_upper == 'S' || init_upper == 'N') {
            // For INIT='S' or 'N', minimum workspace depends on optimization parameters
            ldwork_calc = MAX(30, nsmp_ll * l_ll + 2);
        } else {
            // For INIT='L' or 'B', minimum workspace for linear estimation
            ldwork_calc = MAX(1, 2 * (m_ll + l_ll) * nobr_ll);
        }
    }
    ldwork_calc = MAX(1, ldwork_calc);

    if (ldwork_calc > 0) {
        dwork_alloc = (double*)malloc((size_t)ldwork_calc * sizeof(double));
        CHECK_ALLOC(dwork_alloc);
    } else {
        dwork_alloc = NULL; // Should not happen
    }


    // 4. Memory allocation for column-major copies (if row_major)
    if (nsmp > 0 && m > 0) u_size = (size_t)nsmp * m; else u_size = 0;
    if (nsmp > 0 && l > 0) y_size = (size_t)nsmp * l; else y_size = 0;

    if (row_major) {
        if (u_size > 0) { u_cm = (double*)malloc(u_size * sizeof(double)); CHECK_ALLOC(u_cm); }
        if (y_size > 0) { y_cm = (double*)malloc(y_size * sizeof(double)); CHECK_ALLOC(y_cm); }
    }

    // 5. Prepare Fortran parameters and perform conversions
    ldu_f = (m == 0) ? 1 : MAX(1, nsmp);
    ldy_f = (l == 0) ? 1 : MAX(1, nsmp); // L can be 0 if INIT=S/N

    if (row_major) {
        if (u_size > 0) { slicot_transpose_to_fortran_with_ld(u, u_cm, nsmp, m, ldu, ldu_f, sizeof(double)); u_f_ptr = u_cm; }
        else { u_f_ptr = NULL; }
        if (y_size > 0) { slicot_transpose_to_fortran_with_ld(y, y_cm, nsmp, l, ldy, ldy_f, sizeof(double)); y_f_ptr = y_cm; }
        else { y_f_ptr = NULL; }
    } else { // Column-major C
        if (u_size == 0) u_f_ptr = NULL; ldu_f = ldu;
        if (y_size == 0) y_f_ptr = NULL; ldy_f = ldy;
    }
    if (u_f_ptr != NULL && ldu_f < 1) ldu_f = 1;
    if (y_f_ptr != NULL && ldy_f < 1) ldy_f = 1;


    // 7. Call Fortran function
    for (int attempt = 0; attempt < 2; ++attempt) {
        F77_FUNC(ib03ad, IB03AD)(&init_upper, &alg_upper, &stor_upper,
                                 &nobr, &m, &l, &nsmp,
                                 n_ptr, /* Pass pointer for N */
                                 &nn, &itmax1, &itmax2, &nprint,
                                 u_f_ptr, &ldu_f,
                                 y_f_ptr, &ldy_f,
                                 x, lx_ptr, /* Pass pointer for LX */
                                 &tol1, &tol2,
                                 iwork_alloc, dwork_alloc, &ldwork_calc,
                                 &local_iwarn, &info,
                                 init_len, alg_len, stor_len);

        if (info != -23 || dwork_alloc == NULL) {
            break;
        }

        double required_ldwork_raw = dwork_alloc[0];
        int required_ldwork = (int)ceil(required_ldwork_raw);
        if (!isfinite(required_ldwork_raw) || required_ldwork <= ldwork_calc || required_ldwork <= 0) {
            required_ldwork = ldwork_calc + 64;
        }

        free(dwork_alloc);
        ldwork_calc = required_ldwork;
        dwork_alloc = (double*)malloc((size_t)ldwork_calc * sizeof(double));
        CHECK_ALLOC(dwork_alloc);

        info = 0;
        local_iwarn = 0;
        *n_ptr = original_n;
        *lx_ptr = original_lx;
    }

    if (iwarn_ptr != NULL) {
        *iwarn_ptr = local_iwarn;
    }

    // Populate output summaries if requested
    if (info == 0 || info == -23 || local_iwarn != 0) { // Also populate on warning or if min_ldwork returned
        if (out_iwork_summary != NULL && iwork_alloc != NULL && liwork_calc >= 3) {
            out_iwork_summary[0] = iwork_alloc[0]; // func evals
            out_iwork_summary[1] = iwork_alloc[1]; // jac evals
            out_iwork_summary[2] = iwork_alloc[2]; // rcond count
        }
        if (out_dwork_summary != NULL && dwork_alloc != NULL) {
            if (ldwork_calc >= 1) out_dwork_summary[0] = dwork_alloc[0]; // opt/min LDWORK
            if (ldwork_calc >= 2) out_dwork_summary[1] = dwork_alloc[1]; // residual norm
            if (ldwork_calc >= 3) out_dwork_summary[2] = dwork_alloc[2]; // iterations
            if (ldwork_calc >= 4) out_dwork_summary[3] = dwork_alloc[3]; // cg iterations
            if (ldwork_calc >= 5) out_dwork_summary[4] = dwork_alloc[4]; // Levenberg factor
            if (init_upper == 'S' || init_upper == 'B') {
                if (ldwork_calc >= 6)  out_dwork_summary[5] = dwork_alloc[5];
                if (ldwork_calc >= 7)  out_dwork_summary[6] = dwork_alloc[6];
                if (ldwork_calc >= 8)  out_dwork_summary[7] = dwork_alloc[7];
                if (ldwork_calc >= 9)  out_dwork_summary[8] = dwork_alloc[8];
                if (ldwork_calc >= 10) out_dwork_summary[9] = dwork_alloc[9];
            }
            // For DWORK(11:10+IWORK(3)), the caller needs to manage a potentially larger array
            // For simplicity, this wrapper only returns the first 5 or 10 elements.
            // A more complex wrapper could return more based on IWORK(3).
        }
    }

    // 8. Convert results back to row-major: Not applicable for X (1D) or U, Y (inputs)

cleanup:
    free(iwork_alloc);
    free(dwork_alloc);
    if (row_major) {
        free(u_cm);
        free(y_cm);
    }
    return info;
}
