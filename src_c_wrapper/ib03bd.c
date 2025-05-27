/**
 * @file ib03bd.c
 * @brief C wrapper for SLICOT routine IB03BD.
 * @details Wiener system identification using a MINPACK-like Levenberg-Marquardt algorithm.
 * Workspace (IWORK, DWORK) is allocated internally by this wrapper.
 * Input/output matrix format for U and Y is handled via the row_major parameter.
 * C-level validation for zero dimensions aims to align with Fortran routine behavior.
 */

#include <stdlib.h> // For malloc, free
#include <ctype.h>  // For toupper
#include <stddef.h> // For size_t
#include <math.h>   // For MAX/MIN if needed (often provided by slicot_utils.h)
#include <stdio.h>  // For error logging (optional)

#include "ib03bd.h"       // Public header for this wrapper (defines MAX_EXPECTED_RCONDS_IN_DWORK_HDR_IB03BD)
#include "slicot_utils.h" // Provides CHECK_ALLOC, SLICOT_MEMORY_ERROR, MAX/MIN, transpose functions etc.
#include "slicot_f77.h"   // Provides F77_FUNC macro for Fortran name mangling

/* External Fortran routine declaration */
extern void F77_FUNC(ib03bd, IB03BD)(
    const char* init, const int* nobr, const int* m, const int* l, const int* nsmp,
    int* n, /* N is INOUT */
    const int* nn, const int* itmax1, const int* itmax2, const int* nprint,
    const double* u, const int* ldu,
    const double* y, const int* ldy,
    double* x, int* lx, /* LX is INOUT */
    const double* tol1, const double* tol2,
    int* iwork, double* dwork, const int* ldwork,
    int* iwarn, int* info,
    int init_len);

/* C wrapper function definition */
SLICOT_EXPORT
int slicot_ib03bd(
    char init_char, int nobr_in, int m_in, int l_in, int nsmp_in,
    int* n_ptr, 
    int nn_in, 
    int itmax1_in, int itmax2_in, int nprint_in,
    const double* u, int ldu,
    const double* y, int ldy,
    double* x, int* lx_ptr, 
    double tol1, double tol2,
    int* out_iwork_summary, 
    double* out_dwork_summary, 
    int* iwarn_ptr,
    int row_major)
{
    // 1. Variable declarations
    int info = 0;
    int local_iwarn = 0;
    char init_upper;

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

    int n_val = *n_ptr; 
    int lx_val = *lx_ptr;

    // 2. Input parameter validation
    init_upper = toupper(init_char);

    if (init_upper != 'L' && init_upper != 'S' && init_upper != 'B' && init_upper != 'N') { info = -1; goto cleanup; }

    if (init_upper == 'L' || init_upper == 'B') {
        if (nobr_in <= 0) { info = -2; goto cleanup; } // NOBR is 2nd arg in Fortran, 2nd in C after char
    }
    if (m_in < 0) { info = -3; goto cleanup; }
    if (l_in < 0) { info = -4; goto cleanup; } 
    if (init_upper == 'L' || init_upper == 'B') {
        if (l_in <= 0) { info = -4; goto cleanup; } 
    }

    if (nsmp_in < 0) { info = -5; goto cleanup; }
    if (init_upper == 'L' || init_upper == 'B') {
        if (nsmp_in < 2 * (m_in + l_in + 1) * nobr_in - 1) { info = -5; goto cleanup; }
    }
    
    // N validation (arg 6)
    if (init_upper == 'S' || init_upper == 'N') {
        if (n_val < 0) { info = -6; goto cleanup; }
    } else { // INIT = 'L' or 'B'
        if (n_val >= 0) { 
            if (n_val >= nobr_in && nobr_in > 0) { info = -6; goto cleanup; } 
            if (n_val == 0 && nobr_in > 0 ) { info = -6; goto cleanup; } 
        }
    }

    if (nn_in < 0) { info = -7; goto cleanup; } // NN is arg 7

    if (init_upper == 'S' || init_upper == 'B') { // ITMAX1 is arg 8
        if (itmax1_in < 0) { info = -8; goto cleanup; }
    }
    if (itmax2_in < 0) { info = -9; goto cleanup; } // ITMAX2 is arg 9
    // NPRINT (arg 10) can be any integer.

    // U (arg 11), LDU (arg 12)
    if (m_in > 0) {
        if (u == NULL && nsmp_in > 0) { info = -11; goto cleanup; }
        if (u != NULL || nsmp_in > 0) {
             if (row_major) { 
                if (ldu < m_in) { info = -12; goto cleanup; }
            } else { 
                if (ldu < MAX(1,nsmp_in)) { info = -12; goto cleanup; }
            }
        }
    } else { 
        if (u != NULL && ldu < 1) { info = -12; goto cleanup; }
    }

    // Y (arg 13), LDY (arg 14)
    if (l_in > 0) {
        if (y == NULL && nsmp_in > 0) { info = -13; goto cleanup; }
         if (y != NULL || nsmp_in > 0) {
            if (row_major) { 
                if (ldy < l_in) { info = -14; goto cleanup; }
            } else { 
                if (ldy < MAX(1,nsmp_in)) { info = -14; goto cleanup; }
            }
        }
    } else { 
         if (y != NULL && ldy < 1) { info = -14; goto cleanup; }
    }
    
    // X (arg 15), LX (arg 16)
    int n_for_lx_calc = (n_val < 0 && (init_upper == 'L' || init_upper == 'B')) ? MAX(0, nobr_in - 1) : MAX(0, n_val);
    int bsn_lx = nn_in * (MAX(0,l_in) + 2) + 1;
    int lths_lx = n_for_lx_calc * (MAX(0,l_in) + MAX(0,m_in) + 1) + MAX(0,l_in) * MAX(0,m_in);
    int required_nx = bsn_lx * MAX(0,l_in) + lths_lx;
    if (required_nx < 0) required_nx = 0; 

    if (x == NULL && required_nx > 0) { info = -15; goto cleanup; } 
    if (n_val >= 0) { 
        if (lx_val < required_nx) { info = -16; goto cleanup; }
    }
    if (x != NULL && lx_val < 1 && required_nx > 0) { info = -16; goto cleanup; }
    if (required_nx == 0 && x != NULL && lx_val < 1) { info = -16; goto cleanup; }
    // TOL1 (arg 17), TOL2 (arg 18) can be any double.

    if (info != 0) { goto cleanup; }

    // 3. Internal Workspace Allocation
    int n_eff_ws = (n_val < 0 && (init_upper == 'L' || init_upper == 'B')) ? MAX(0, nobr_in - 1) : MAX(0, n_val);

    long long m_ll = m_in;
    long long l_ll = l_in;
    long long nsmp_ll = nsmp_in;
    long long nobr_ll = nobr_in;
    long long n_ll_eff_ws = n_eff_ws;
    long long nn_ll = nn_in;
    
    long long bsn_ws  = nn_ll * (l_ll + 2) + 1;
    long long lths_ws = n_ll_eff_ws * (l_ll + m_ll + 1) + l_ll * m_ll;
    long long nx_ws = bsn_ws * l_ll + lths_ws;
    if (nx_ws < 0) nx_ws = 0;

    long long liw1_f = 0, liw2_f = 0, liw3_f = 0;
    if (init_upper == 'S' || init_upper == 'N') {
        liw1_f = 0; liw2_f = 0;
    } else {
        liw1_f = m_ll + l_ll;
        liw2_f = MAX(m_ll * nobr_ll + n_ll_eff_ws, m_ll * (n_ll_eff_ws + l_ll));
    }
    if (init_upper == 'S' || init_upper == 'B') {
        liw3_f = 3LL + MAX(nn_ll * (l_ll + 2) + 2LL, nx_ws + l_ll);
    } else { // INIT = 'L' or 'N'
        liw3_f = 3LL + nx_ws + l_ll;
    }
    liwork_calc = (int)MAX(MAX(liw1_f, liw2_f), liw3_f);
    liwork_calc = MAX(1, liwork_calc); // IWORK must be at least 1, though formulas likely give >0.

    if (liwork_calc > 0) {
        iwork_alloc = (int*)malloc((size_t)liwork_calc * sizeof(int));
        CHECK_ALLOC(iwork_alloc);
    } else { 
        iwork_alloc = NULL;
    }
    
    long long ldw1_d = 0, ldw2_d = 0, ldw3_d = 0, ldw4_d = 0;
    long long ldw5_d = 0, ldw6_d = 0, ldw7_d = 0;
    long long lw1_d = 0, lw2_d = 0, lw3_d = 0, lw4_d = 0;

    if (init_upper == 'S' || init_upper == 'N') {
        lw1_d = 0;
    } else {
        ldw1_d = MAX( 2LL*(l_ll*nobr_ll-l_ll)*n_ll_eff_ws+2LL*n_ll_eff_ws, (l_ll*nobr_ll-l_ll)*n_ll_eff_ws+n_ll_eff_ws*n_ll_eff_ws+7LL*n_ll_eff_ws );
        ldw1_d = MAX( ldw1_d, l_ll*nobr_ll*n_ll_eff_ws +
                        MAX( (l_ll*nobr_ll-l_ll)*n_ll_eff_ws+2LL*n_ll_eff_ws + (2LL*m_ll+l_ll)*nobr_ll+l_ll,
                             2LL*(l_ll*nobr_ll-l_ll)*n_ll_eff_ws+n_ll_eff_ws*n_ll_eff_ws+8LL*n_ll_eff_ws ) );
        ldw1_d = MAX( ldw1_d, MAX( n_ll_eff_ws+4LL*(m_ll*nobr_ll+n_ll_eff_ws)+1LL, m_ll*nobr_ll+3LL*n_ll_eff_ws+l_ll ) );

        if (m_ll == 0) ldw2_d = 0;
        else ldw2_d = l_ll*nobr_ll*n_ll_eff_ws + m_ll*nobr_ll*(n_ll_eff_ws+l_ll)*(m_ll*(n_ll_eff_ws+l_ll)+1LL) +
                    MAX( (n_ll_eff_ws+l_ll)*(n_ll_eff_ws+l_ll), 4LL*m_ll*(n_ll_eff_ws+l_ll)+1LL );
        
        ldw3_d = nsmp_ll*l_ll*(n_ll_eff_ws+1LL) + 2LL*n_ll_eff_ws + MAX( 2LL*n_ll_eff_ws*n_ll_eff_ws, 4LL*n_ll_eff_ws );
        ldw4_d = n_ll_eff_ws*(n_ll_eff_ws+1LL) + 2LL*n_ll_eff_ws +
                   MAX( n_ll_eff_ws*l_ll*(n_ll_eff_ws+1LL) + 2LL*n_ll_eff_ws*n_ll_eff_ws + l_ll*n_ll_eff_ws, 4LL*n_ll_eff_ws );
        ldw5_d = nsmp_ll*l_ll + (n_ll_eff_ws+l_ll)*(n_ll_eff_ws+m_ll) + 3LL*n_ll_eff_ws+m_ll+l_ll;
        ldw6_d = nsmp_ll*l_ll + (n_ll_eff_ws+l_ll)*(n_ll_eff_ws+m_ll) + n_ll_eff_ws +
               MAX(1LL, n_ll_eff_ws*n_ll_eff_ws*l_ll + n_ll_eff_ws*l_ll + n_ll_eff_ws );
        ldw6_d = MAX(ldw6_d, nsmp_ll*l_ll + (n_ll_eff_ws+l_ll)*(n_ll_eff_ws+m_ll) + n_ll_eff_ws + n_ll_eff_ws*n_ll_eff_ws +
                         MAX(n_ll_eff_ws*n_ll_eff_ws + n_ll_eff_ws*MAX(n_ll_eff_ws,l_ll) + 6LL*n_ll_eff_ws + MIN(n_ll_eff_ws,l_ll),
                             n_ll_eff_ws*m_ll) );

        lw1_d = MAX( 2LL*(m_ll+l_ll)*nobr_ll*(2LL*(m_ll+l_ll)*(nobr_ll+1LL)+3LL) + l_ll*nobr_ll,
                   4LL*(m_ll+l_ll)*nobr_ll*(m_ll+l_ll)*nobr_ll + (n_ll_eff_ws+l_ll)*(n_ll_eff_ws+m_ll) + MAX( ldw1_d, ldw2_d ) );
        lw1_d = MAX( lw1_d, (n_ll_eff_ws+l_ll)*(n_ll_eff_ws+m_ll) + n_ll_eff_ws + n_ll_eff_ws*n_ll_eff_ws + 2LL + n_ll_eff_ws*(n_ll_eff_ws+m_ll+l_ll) +
                        MAX( MAX( 5LL*n_ll_eff_ws, 2LL ), MAX( MIN( ldw3_d, ldw4_d ), MAX(ldw5_d, ldw6_d) ) ) );
    }

    if (init_upper == 'L' || init_upper == 'N') {
        lw2_d = 0; lw3_d = 0; // LW2 and LW3 from doc are for INIT=S/B
    } else { // INIT = 'S' or 'B'
        lw2_d = nsmp_ll*l_ll + bsn_ws +
                MAX( 4LL, nsmp_ll +
                        MAX( nsmp_ll*bsn_ws + MAX( 2LL*nn_ll, 5LL*bsn_ws + 1LL ),
                             bsn_ws*bsn_ws + bsn_ws +
                             MAX( nsmp_ll + 2LL*nn_ll, 5LL*bsn_ws ) ) );
        
        if (m_ll > 0) ldw7_d = nsmp_ll*l_ll + (n_ll_eff_ws+l_ll)*(n_ll_eff_ws+m_ll) + 3LL*n_ll_eff_ws+m_ll+l_ll;
        else ldw7_d = nsmp_ll*l_ll + (n_ll_eff_ws+l_ll)*n_ll_eff_ws + 2LL*n_ll_eff_ws+l_ll;
        lw3_d = MAX( ldw7_d, nsmp_ll*l_ll + (n_ll_eff_ws+l_ll)*(2LL*n_ll_eff_ws+m_ll) + 2LL*n_ll_eff_ws );
    }
    
    long long l0_d, l1_d_val, l2_d_val, l3_d_val; 
    if (m_ll > 0) l0_d = MAX( n_ll_eff_ws*(n_ll_eff_ws+l_ll), n_ll_eff_ws+m_ll+l_ll );
    else l0_d = MAX( n_ll_eff_ws*(n_ll_eff_ws+l_ll), l_ll );
    l1_d_val = nsmp_ll*l_ll + MAX( 2LL*nn_ll, (n_ll_eff_ws+l_ll)*(n_ll_eff_ws+m_ll) + 2LL*n_ll_eff_ws + l0_d);

    if (l_ll <= 1 || bsn_ws == 0) l2_d_val = 4LL*nx_ws + 1LL;
    else l2_d_val = bsn_ws + MAX(3LL*bsn_ws+1LL,lths_ws);
    
    if (nsmp_ll > bsn_ws) l2_d_val = MAX(l2_d_val,4LL*lths_ws+1LL);
    if (bsn_ws < nsmp_ll && nsmp_ll < 2LL*bsn_ws) l2_d_val = MAX(l2_d_val,(nsmp_ll-bsn_ws)*(l_ll-1LL));
    
    if (l_ll <= 1 || bsn_ws == 0) l3_d_val = 4LL*nx_ws;
    else l3_d_val = lths_ws*bsn_ws + 2LL*nx_ws + 2LL*MAX(bsn_ws,lths_ws);


    lw4_d = nsmp_ll*l_ll + nx_ws +
            MAX( 4LL, nsmp_ll*l_ll +
                    MAX( nsmp_ll*l_ll*( bsn_ws + lths_ws ) +
                         MAX( nsmp_ll*l_ll + l1_d_val, l2_d_val + nx_ws ),
                              nx_ws*( bsn_ws + lths_ws ) +
                              nx_ws +
                              MAX( nsmp_ll*l_ll + l1_d_val, nx_ws + l3_d_val ) ) );
    
    ldwork_calc = (int)MAX(lw1_d, MAX(lw2_d, MAX(lw3_d, lw4_d)));
    ldwork_calc = MAX(1, ldwork_calc); 

    if (ldwork_calc > 0) {
        dwork_alloc = (double*)malloc((size_t)ldwork_calc * sizeof(double));
        CHECK_ALLOC(dwork_alloc);
    } else {
        dwork_alloc = NULL; 
    }

    if ((init_upper == 'S' || init_upper == 'B') && dwork_alloc != NULL && ldwork_calc >=4) {
        if(out_dwork_summary != NULL) { 
            if (dwork_alloc != NULL) { 
                dwork_alloc[0] = out_dwork_summary[0];
                dwork_alloc[1] = out_dwork_summary[1];
                dwork_alloc[2] = out_dwork_summary[2];
                dwork_alloc[3] = out_dwork_summary[3];
            }
        } else { 
            if (dwork_alloc != NULL) {
                dwork_alloc[0] = 1998.0; dwork_alloc[1] = 1999.0; dwork_alloc[2] = 2000.0; dwork_alloc[3] = 2001.0;
            }
        }
    }

    // 4. Memory allocation for column-major copies (if row_major)
    if (nsmp_in > 0 && m_in > 0) u_size = (size_t)nsmp_in * m_in; else u_size = 0;
    if (nsmp_in > 0 && l_in > 0) y_size = (size_t)nsmp_in * l_in; else y_size = 0;

    if (row_major) {
        if (u_size > 0) { u_cm = (double*)malloc(u_size * sizeof(double)); CHECK_ALLOC(u_cm); }
        if (y_size > 0) { y_cm = (double*)malloc(y_size * sizeof(double)); CHECK_ALLOC(y_cm); }
    }

    // 5. Prepare Fortran parameters and perform conversions
    ldu_f = (m_in == 0) ? 1 : MAX(1, nsmp_in);
    ldy_f = (l_in == 0) ? 1 : MAX(1, nsmp_in);

    if (row_major) {
        if (u_size > 0 && u != NULL && u_cm != NULL) { 
            slicot_transpose_to_fortran_with_ld(u, u_cm, nsmp_in, m_in, ldu, ldu_f, sizeof(double)); 
            u_f_ptr = u_cm; 
        } else if (u_size == 0) {
             u_f_ptr = NULL;
        } 
        
        if (y_size > 0 && y != NULL && y_cm != NULL) { 
            slicot_transpose_to_fortran_with_ld(y, y_cm, nsmp_in, l_in, ldy, ldy_f, sizeof(double)); 
            y_f_ptr = y_cm; 
        } else if (y_size == 0) {
            y_f_ptr = NULL;
        }
    } else { 
        if (u_size == 0) u_f_ptr = NULL; ldu_f = ldu; 
        if (y_size == 0) y_f_ptr = NULL; ldy_f = ldy; 
    }
    
    if (u_f_ptr != NULL && ldu_f < 1) ldu_f = 1;
    if (y_f_ptr != NULL && ldy_f < 1) ldy_f = 1;

    // 7. Call Fortran function
    F77_FUNC(ib03bd, IB03BD)(&init_upper, &nobr_in, &m_in, &l_in, &nsmp_in,
                             n_ptr, 
                             &nn_in, &itmax1_in, &itmax2_in, &nprint_in,
                             u_f_ptr, &ldu_f,
                             y_f_ptr, &ldy_f,
                             x, lx_ptr, 
                             &tol1, &tol2,
                             iwork_alloc, dwork_alloc, &ldwork_calc,
                             &local_iwarn, &info,
                             init_len);

    if (iwarn_ptr != NULL) {
        *iwarn_ptr = local_iwarn;
    }

    if (info == 0 || info == -21 || local_iwarn != 0) { 
        if (out_iwork_summary != NULL && iwork_alloc != NULL && liwork_calc >= 3) {
            out_iwork_summary[0] = iwork_alloc[0]; 
            out_iwork_summary[1] = iwork_alloc[1]; 
            out_iwork_summary[2] = iwork_alloc[2]; 
        } else if (out_iwork_summary != NULL) { 
             if (liwork_calc >=1) out_iwork_summary[0] = (iwork_alloc != NULL) ? iwork_alloc[0] : -1; else out_iwork_summary[0] = -1;
             if (liwork_calc >=2 && out_iwork_summary != NULL) out_iwork_summary[1] = (iwork_alloc != NULL) ? iwork_alloc[1] : -1; else if (out_iwork_summary != NULL && liwork_calc >=1) out_iwork_summary[1] = -1; else if (out_iwork_summary != NULL) out_iwork_summary[1] = -1;
             if (liwork_calc >=3 && out_iwork_summary != NULL) out_iwork_summary[2] = (iwork_alloc != NULL) ? iwork_alloc[2] : -1; else if (out_iwork_summary != NULL && liwork_calc >=1) out_iwork_summary[2] = -1; else if (out_iwork_summary != NULL) out_iwork_summary[2] = -1;
        }

        if (out_dwork_summary != NULL && dwork_alloc != NULL) {
            int num_rconds = 0;
            if (iwork_alloc != NULL && liwork_calc >=3 && (init_upper == 'L' || init_upper == 'B')) {
                num_rconds = iwork_alloc[2];
            }
            // DWORK(1-4) main results, DWORK(5-8) init results, DWORK(9:8+IWORK(3)) RCONDS
            int elements_to_copy = 8 + ((num_rconds > 0) ? num_rconds : 0);
            elements_to_copy = MIN(elements_to_copy, ldwork_calc); 
            elements_to_copy = MIN(elements_to_copy, 8 + MAX_EXPECTED_RCONDS_IN_DWORK_HDR_IB03BD);


            for(int i=0; i < elements_to_copy; ++i) {
                 out_dwork_summary[i] = dwork_alloc[i];
            }
            for(int i=elements_to_copy; i < 8 + MAX_EXPECTED_RCONDS_IN_DWORK_HDR_IB03BD; ++i) {
                 if (i < (int)(sizeof(out_dwork_summary)/sizeof(out_dwork_summary[0]))) { // Check bounds
                    out_dwork_summary[i] = -999.0; 
                 }
            }
        }
    }

cleanup:
    free(iwork_alloc);
    free(dwork_alloc);
    if (row_major) {
        free(u_cm);
        free(y_cm);
    }
    return info;
}
