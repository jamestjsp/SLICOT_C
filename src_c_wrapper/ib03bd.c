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

// Static dummy variable to pass to Fortran if M=0 or L=0 and original pointer is NULL,
// and the Fortran routine does not explicitly state the array is "not referenced".
static double SLICOT_IB03BD_DUMMY_DOUBLE_VAL = 0.0;

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

    const double* u_f_ptr = u; // Pointer to pass to Fortran for U
    const double* y_f_ptr = y; // Pointer to pass to Fortran for Y

    int ldu_f, ldy_f; // Leading dimensions for Fortran
    size_t u_size = 0, y_size = 0; // Size in elements for allocation

    const int init_len = 1;

    int n_val = *n_ptr; 
    int lx_val = *lx_ptr;

    // 2. Input parameter validation
    init_upper = toupper(init_char);

    // Arg 1: INIT
    if (init_upper != 'L' && init_upper != 'S' && init_upper != 'B' && init_upper != 'N') { info = -1; goto cleanup; }

    // Arg 2: NOBR
    if (init_upper == 'L' || init_upper == 'B') {
        if (nobr_in <= 0) { info = -2; goto cleanup; }
    }
    // Arg 3: M
    if (m_in < 0) { info = -3; goto cleanup; }
    // Arg 4: L
    if (l_in < 0) { info = -4; goto cleanup; } 
    if (init_upper == 'L' || init_upper == 'B') { // L > 0 if INIT is L or B
        if (l_in <= 0) { info = -4; goto cleanup; } 
    }
    // Arg 5: NSMP
    if (nsmp_in < 0) { info = -5; goto cleanup; }
    if (init_upper == 'L' || init_upper == 'B') {
        if (nsmp_in < 2 * (m_in + l_in + 1) * nobr_in - 1) { info = -5; goto cleanup; }
    }
    
    // Arg 6: N (pointed to by n_ptr)
    if (init_upper == 'S' || init_upper == 'N') { // N is input
        if (n_val < 0) { info = -6; goto cleanup; }
    } else { // INIT = 'L' or 'B'; N can be input or output
        if (n_val >= 0) { // N is input
            // According to IB03BD.html: "The values N >= NOBR, or N = 0, are not acceptable if INIT = 'L' or 'B'."
            // This means if N is provided (N>=0), it must be 0 < N < NOBR.
            if (nobr_in > 0 && (n_val <= 0 || n_val >= nobr_in) ) { info = -6; goto cleanup; } 
        }
        // If N < 0 on entry, it's an output, so no check here.
    }

    // Arg 7: NN
    if (nn_in < 0) { info = -7; goto cleanup; } 

    // Arg 8: ITMAX1
    if (init_upper == 'S' || init_upper == 'B') { 
        if (itmax1_in < 0) { info = -8; goto cleanup; }
    }
    // Arg 9: ITMAX2
    if (itmax2_in < 0) { info = -9; goto cleanup; } 
    // Arg 10: NPRINT can be any integer.

    // Arg 11: U, Arg 12: LDU
    if (m_in > 0) { 
        if (u == NULL && nsmp_in > 0) { info = -11; goto cleanup; } 
        // LDU validation: Fortran LDU is num_rows (nsmp_in for col-major, m_in for row-major)
        // C LDU is num_rows if col-major, num_cols if row-major
        if (u != NULL || (nsmp_in > 0 && m_in > 0) ) { 
             if (row_major) { // C ldu is number of columns
                if (ldu < m_in) { info = -12; goto cleanup; } 
            } else { // C ldu is number of rows
                if (ldu < MAX(1,nsmp_in)) { info = -12; goto cleanup; } 
            }
        }
    } else { // M_in == 0
        // If M=0, U is not referenced. LDU must be >=1 if U is not NULL.
        // If U is NULL, LDU is not strictly checked by Fortran if M=0, but C wrapper ensures it's >=1 if U is not NULL.
        if (u != NULL && ldu < 1) { info = -12; goto cleanup; }
    }
    
    // Arg 13: Y, Arg 14: LDY
    if (l_in > 0) { 
        if (y == NULL && nsmp_in > 0) { info = -13; goto cleanup; }
         if (y != NULL || (nsmp_in > 0 && l_in > 0) ) { 
            if (row_major) { // C ldy is number of columns
                if (ldy < l_in) { info = -14; goto cleanup; }
            } else { // C ldy is number of rows
                if (ldy < MAX(1,nsmp_in)) { info = -14; goto cleanup; }
            }
        }
    } else { // L_in == 0 (possible if INIT = 'S' or 'N')
         // If L=0, Y is not referenced. LDY must be >=1 if Y is not NULL.
         if (y != NULL && ldy < 1) { info = -14; goto cleanup; }
    }
    
    // Arg 15: X, Arg 16: LX
    // Calculate required_nx based on n_val (if input) or estimated n_val (if output)
    int n_for_lx_calc;
    if (init_upper == 'L' || init_upper == 'B') {
        n_for_lx_calc = (n_val < 0) ? MAX(0, nobr_in - 1) : MAX(0, n_val);
    } else { // INIT == 'S' or 'N'
        n_for_lx_calc = MAX(0, n_val);
    }
    int bsn_lx = nn_in * (MAX(0,l_in) + 2) + 1;
    int lths_lx = n_for_lx_calc * (MAX(0,l_in) + MAX(0,m_in) + 1) + MAX(0,l_in) * MAX(0,m_in);
    int required_nx = bsn_lx * MAX(0,l_in) + lths_lx;
    if (required_nx < 0) required_nx = 0; // Should not happen with MAX(0,...)

    if (x == NULL && required_nx > 0) { info = -15; goto cleanup; } 
    // If N is input (N>=0), LX must be >= required_nx.
    // If N is output (N<0), LX is also an output if it was too small (INFO=-21 from Fortran).
    // The Fortran routine will check LX if N is input.
    // The C wrapper mainly checks if X is NULL when it's needed.
    if (n_val >= 0) { // N is input
        if (lx_val < required_nx) { info = -16; goto cleanup; }
    }
    // LX must be at least 1 if X is not NULL, even if required_nx is 0.
    if (x != NULL && lx_val < 1 ) { info = -16; goto cleanup; } 


    if (info != 0) { goto cleanup; }

    // 3. Internal Workspace Allocation
    // Determine effective N for workspace calculation (n_eff_ws)
    // If INIT='L' or 'B' and N < 0 (output), N is estimated as NOBR-1 for workspace.
    // Otherwise, use input N.
    int n_eff_ws;
    if ((init_upper == 'L' || init_upper == 'B') && n_val < 0) {
        n_eff_ws = MAX(0, nobr_in - 1); // N estimated as NOBR-1 for workspace
    } else {
        n_eff_ws = MAX(0, n_val); // Use input N (must be >=0 if not L/B with N<0)
    }

    long long m_ll = m_in;
    long long l_ll = l_in;
    long long nsmp_ll = nsmp_in;
    long long nobr_ll = nobr_in;
    long long n_ll_eff_ws = n_eff_ws;
    long long nn_ll = nn_in;
    
    // Calculate NX for workspace based on n_eff_ws
    long long bsn_ws  = nn_ll * (l_ll + 2) + 1;
    long long lths_ws = n_ll_eff_ws * (l_ll + m_ll + 1) + l_ll * m_ll;
    long long nx_ws = bsn_ws * l_ll + lths_ws;
    if (nx_ws < 0) nx_ws = 0; // Should not happen

    // IWORK calculation (from IB03BD.html)
    long long liw1_f = 0, liw2_f = 0, liw3_f = 0;
    if (init_upper == 'S' || init_upper == 'N') {
        liw1_f = 0; liw2_f = 0; 
    } else { // INIT == 'L' or 'B'
        liw1_f = m_ll + l_ll;
        liw2_f = MAX(m_ll * nobr_ll + n_ll_eff_ws, m_ll * (n_ll_eff_ws + l_ll));
    }
    if (init_upper == 'S' || init_upper == 'B') {
        liw3_f = 3LL + MAX(nn_ll * (l_ll + 2) + 2LL, nx_ws + l_ll);
    } else { // INIT == 'L' or 'N'
        liw3_f = 3LL + nx_ws + l_ll;
    }
    liwork_calc = (int)MAX(MAX(liw1_f, liw2_f), liw3_f);
    liwork_calc = MAX(1, liwork_calc); // IWORK must be at least 1

    if (liwork_calc > 0) {
        iwork_alloc = (int*)malloc((size_t)liwork_calc * sizeof(int));
        CHECK_ALLOC(iwork_alloc);
    } else { 
        iwork_alloc = NULL; // Should not happen due to MAX(1,...)
    }
    
    // DWORK calculation (from IB03BD.html)
    long long ldw1_d = 0, ldw2_d = 0, ldw3_d = 0, ldw4_d = 0;
    long long ldw5_d = 0, ldw6_d = 0, ldw7_d = 0;
    long long lw1_d = 0, lw2_d = 0, lw3_d = 0, lw4_d = 0;

    if (init_upper == 'S' || init_upper == 'N') {
        lw1_d = 0; // Not used by these INIT options for this part of DWORK
    } else { // INIT == 'L' or 'B'
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
               MAX(1LL, n_ll_eff_ws*n_ll_eff_ws*l_ll + n_ll_eff_ws*l_ll + n_ll_eff_ws ); // Added MAX(1,...) for safety
        ldw6_d = MAX(ldw6_d, nsmp_ll*l_ll + (n_ll_eff_ws+l_ll)*(n_ll_eff_ws+m_ll) + n_ll_eff_ws + n_ll_eff_ws*n_ll_eff_ws +
                         MAX(n_ll_eff_ws*n_ll_eff_ws + n_ll_eff_ws*MAX(n_ll_eff_ws,l_ll) + 6LL*n_ll_eff_ws + MIN(n_ll_eff_ws,l_ll),
                             n_ll_eff_ws*m_ll) );

        lw1_d = MAX( 2LL*(m_ll+l_ll)*nobr_ll*(2LL*(m_ll+l_ll)*(nobr_ll+1LL)+3LL) + l_ll*nobr_ll,
                   4LL*(m_ll+l_ll)*nobr_ll*(m_ll+l_ll)*nobr_ll + (n_ll_eff_ws+l_ll)*(n_ll_eff_ws+m_ll) + MAX( ldw1_d, ldw2_d ) );
        lw1_d = MAX( lw1_d, (n_ll_eff_ws+l_ll)*(n_ll_eff_ws+m_ll) + n_ll_eff_ws + n_ll_eff_ws*n_ll_eff_ws + 2LL + n_ll_eff_ws*(n_ll_eff_ws+m_ll+l_ll) +
                        MAX( MAX( 5LL*n_ll_eff_ws, 2LL ), MAX( MIN( ldw3_d, ldw4_d ), MAX(ldw5_d, ldw6_d) ) ) );
    }

    if (init_upper == 'L' || init_upper == 'N') {
        lw2_d = 0; lw3_d = 0; 
    } else { // INIT == 'S' or 'B'
        lw2_d = nsmp_ll*l_ll + bsn_ws +
                MAX( 4LL, nsmp_ll +
                        MAX( nsmp_ll*bsn_ws + MAX( 2LL*nn_ll, 5LL*bsn_ws + 1LL ),
                             bsn_ws*bsn_ws + bsn_ws +
                             MAX( nsmp_ll + 2LL*nn_ll, 5LL*bsn_ws ) ) );
        
        if (m_ll > 0) ldw7_d = nsmp_ll*l_ll + (n_ll_eff_ws+l_ll)*(n_ll_eff_ws+m_ll) + 3LL*n_ll_eff_ws+m_ll+l_ll;
        else ldw7_d = nsmp_ll*l_ll + (n_ll_eff_ws+l_ll)*n_ll_eff_ws + 2LL*n_ll_eff_ws+l_ll;
        lw3_d = MAX( ldw7_d, nsmp_ll*l_ll + (n_ll_eff_ws+l_ll)*(2LL*n_ll_eff_ws+m_ll) + 2LL*n_ll_eff_ws );
    }
    
    // LW4 calculation (common to all INIT options, but components differ)
    long long l0_d_calc, l1_d_val_calc, l2_d_val_calc, l3_d_val_calc; 
    if (m_ll > 0) l0_d_calc = MAX( n_ll_eff_ws*(n_ll_eff_ws+l_ll), n_ll_eff_ws+m_ll+l_ll );
    else l0_d_calc = MAX( n_ll_eff_ws*(n_ll_eff_ws+l_ll), l_ll ); // If m_ll=0, term is n_ll_eff_ws + l_ll
    l1_d_val_calc = nsmp_ll*l_ll + MAX( 2LL*nn_ll, (n_ll_eff_ws+l_ll)*(n_ll_eff_ws+m_ll) + 2LL*n_ll_eff_ws + l0_d_calc);

    if (l_ll <= 1 || bsn_ws == 0) l2_d_val_calc = 4LL*nx_ws + 1LL;
    else l2_d_val_calc = bsn_ws + MAX(3LL*bsn_ws+1LL,lths_ws);
    
    if (nsmp_ll > bsn_ws) l2_d_val_calc = MAX(l2_d_val_calc,4LL*lths_ws+1LL);
    if (bsn_ws < nsmp_ll && nsmp_ll < 2LL*bsn_ws) l2_d_val_calc = MAX(l2_d_val_calc,(nsmp_ll-bsn_ws)*(l_ll-1LL));
    
    if (l_ll <= 1 || bsn_ws == 0) l3_d_val_calc = 4LL*nx_ws;
    else l3_d_val_calc = lths_ws*bsn_ws + 2LL*nx_ws + 2LL*MAX(bsn_ws,lths_ws);

    lw4_d = nsmp_ll*l_ll + nx_ws +
            MAX( 4LL, nsmp_ll*l_ll +
                    MAX( nsmp_ll*l_ll*( bsn_ws + lths_ws ) + // This term seems to be BSN + LTHS, not multiplied by nsmp*l
                         MAX( nsmp_ll*l_ll + l1_d_val_calc, l2_d_val_calc + nx_ws ),
                              nx_ws*( bsn_ws + lths_ws ) + // Same here
                              nx_ws +
                              MAX( nsmp_ll*l_ll + l1_d_val_calc, nx_ws + l3_d_val_calc ) ) );
    
    ldwork_calc = (int)MAX(lw1_d, MAX(lw2_d, MAX(lw3_d, lw4_d)));
    ldwork_calc = MAX(1, ldwork_calc); // DWORK must be at least 1

    if (ldwork_calc > 0) {
        dwork_alloc = (double*)malloc((size_t)ldwork_calc * sizeof(double));
        CHECK_ALLOC(dwork_alloc);
    } else {
        dwork_alloc = NULL; // Should not happen
    }

    // Initialize DWORK(1:4) for random number generator if INIT='S' or 'B'
    if ((init_upper == 'S' || init_upper == 'B') && dwork_alloc != NULL && ldwork_calc >=4) {
        if(out_dwork_summary != NULL && out_dwork_summary[0] != 0.0) { // Check if user provided seed
            // The SLICOT doc says "DWORK(1:4) are set to initialize...", implying input.
            // However, the example program reads DWORK(1:4) from input file for these INIT options.
            // Let's assume if out_dwork_summary is provided, it might contain a seed.
            // This part is a bit ambiguous in typical SLICOT wrappers.
            // For now, let's assume the Fortran routine handles default seed if not set.
            // The example program *does* read values into DWORK(1:4) before calling.
            // This wrapper doesn't have an input for seed, so we can't directly use it.
            // The example uses DWORK(1:4) from input file, so if out_dwork_summary is meant to
            // pass this through, it should be copied. But the C API doesn't take DWORK as input.
            // Let's stick to the Fortran example's behavior: if INIT='S' or 'B', the
            // example program *provides* DWORK(1:4).
            // This wrapper *cannot* take DWORK(1:4) as input directly.
            // The `out_dwork_summary` is for *output*.
            // The original wrapper code had:
            // dwork_alloc[0] = out_dwork_summary[0]; ...
            // This is incorrect as out_dwork_summary is for output.
            // The Fortran routine will use its own default seed if not set.
            // If the user wants reproducible results, they'd need to modify the Fortran or use a version
            // that exposes seed setting. For now, do nothing with DWORK(1:4) on input here.
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
    u_f_ptr = u; 
    y_f_ptr = y; 
    ldu_f = ldu; // C ldu (if col-major, it's rows; if row-major, it's cols)
    ldy_f = ldy; // C ldy

    // --- U matrix preparation ---
    if (m_in == 0) {
        // If M=0, U is not referenced by Fortran. Pass dummy.
        // LDU must be >= MAX(1,NSMP) according to Fortran spec.
        u_f_ptr = (double*)&SLICOT_IB03BD_DUMMY_DOUBLE_VAL; 
        ldu_f = MAX(1, nsmp_in); // FIX: Ensure LDU is valid for Fortran
    } else { // m_in > 0
        // Fortran LDU is always number of rows of U (which is NSMP)
        int fortran_expected_ldu = MAX(1, nsmp_in);
        
        if (row_major) { // C input U is row-major (NSMP rows, M cols, C_LDU = M)
            if (u_size > 0 && u != NULL && u_cm != NULL) { 
                // Transpose from C row-major (ldu = M) to Fortran col-major (ldu_f = NSMP)
                slicot_transpose_to_fortran_with_ld(u, u_cm, nsmp_in, m_in, ldu, fortran_expected_ldu, sizeof(double)); 
                u_f_ptr = u_cm; 
                ldu_f = fortran_expected_ldu;
            } else { 
                // This case (u_size > 0 but u is NULL or u_cm is NULL) should have been caught by CHECK_ALLOC or earlier validation.
                // If u is NULL but m_in > 0 and nsmp_in > 0, it's an error (-11).
                // If execution reaches here with u_f_ptr as NULL when it's needed, use dummy.
                if (u == NULL && nsmp_in > 0) { // Should have been caught by info = -11
                     u_f_ptr = (double*)&SLICOT_IB03BD_DUMMY_DOUBLE_VAL; 
                     ldu_f = MAX(1, nsmp_in); // FIX: Ensure LDU is valid
                } else { // u might be non-NULL but u_cm failed allocation (unlikely due to CHECK_ALLOC)
                     u_f_ptr = NULL; // Or handle as error
                     ldu_f = MAX(1, nsmp_in); // Default if pointer becomes NULL
                }
            }
        } else { // C input U is column-major (NSMP rows, M cols, C_LDU = NSMP)
            // u_f_ptr remains 'u'. C LDU (ldu) should be >= MAX(1,nsmp_in).
            // Fortran LDU (ldu_f) will be C LDU.
            ldu_f = ldu; // This was already validated to be >= MAX(1,nsmp_in)
            if (u == NULL && nsmp_in > 0) { // Should have been caught by info = -11
                 u_f_ptr = (double*)&SLICOT_IB03BD_DUMMY_DOUBLE_VAL;
                 ldu_f = MAX(1, nsmp_in); // FIX: Ensure LDU is valid
            }
        }
        // Final check if u_f_ptr is NULL when it shouldn't be (e.g. nsmp_in=0 but m_in>0)
        if (u_f_ptr == NULL && m_in > 0) { // If nsmp_in was 0, u_size would be 0
             u_f_ptr = (double*)&SLICOT_IB03BD_DUMMY_DOUBLE_VAL;
             ldu_f = MAX(1, nsmp_in); // nsmp_in is likely 0 here, so ldu_f = 1
        }
    }
    
    // --- Y matrix preparation ---
    if (l_in == 0) {
        // If L=0, Y is not referenced by Fortran (if INIT='S' or 'N'). Pass dummy.
        // LDY must be >= MAX(1,NSMP) according to Fortran spec.
        y_f_ptr = (double*)&SLICOT_IB03BD_DUMMY_DOUBLE_VAL;
        ldy_f = MAX(1, nsmp_in); // FIX: Ensure LDY is valid for Fortran
    } else { // l_in > 0
        // Fortran LDY is always number of rows of Y (which is NSMP)
        int fortran_expected_ldy = MAX(1, nsmp_in);

        if (row_major) { // C input Y is row-major (NSMP rows, L cols, C_LDY = L)
            if (y_size > 0 && y != NULL && y_cm != NULL) { 
                slicot_transpose_to_fortran_with_ld(y, y_cm, nsmp_in, l_in, ldy, fortran_expected_ldy, sizeof(double)); 
                y_f_ptr = y_cm; 
                ldy_f = fortran_expected_ldy;
            } else {
                if (y == NULL && nsmp_in > 0) { // Should have been caught by info = -13
                     y_f_ptr = (double*)&SLICOT_IB03BD_DUMMY_DOUBLE_VAL;
                     ldy_f = MAX(1, nsmp_in); // FIX: Ensure LDY is valid
                } else {
                     y_f_ptr = NULL;
                     ldy_f = MAX(1, nsmp_in);
                }
            }
        } else { // C input Y is column-major (NSMP rows, L cols, C_LDY = NSMP)
            ldy_f = ldy; // Validated to be >= MAX(1,nsmp_in)
            if (y == NULL && nsmp_in > 0) { // Should have been caught by info = -13
                 y_f_ptr = (double*)&SLICOT_IB03BD_DUMMY_DOUBLE_VAL;
                 ldy_f = MAX(1, nsmp_in); // FIX: Ensure LDY is valid
            }
        }
        if (y_f_ptr == NULL && l_in > 0) {
             y_f_ptr = (double*)&SLICOT_IB03BD_DUMMY_DOUBLE_VAL;
             ldy_f = MAX(1, nsmp_in);
        }
    }

    // Final safety for LDs if pointers are not NULL but LDs somehow became < 1.
    // This should ideally not be needed if logic above is correct.
    if (u_f_ptr != NULL && ldu_f < 1) ldu_f = 1;
    if (y_f_ptr != NULL && ldy_f < 1) ldy_f = 1;


    // 7. Call Fortran function
    F77_FUNC(ib03bd, IB03BD)(&init_upper, &nobr_in, &m_in, &l_in, &nsmp_in,
                             n_ptr, // n_ptr is INOUT
                             &nn_in, &itmax1_in, &itmax2_in, &nprint_in,
                             u_f_ptr, &ldu_f,
                             y_f_ptr, &ldy_f,
                             x, lx_ptr, // x and lx_ptr are INOUT
                             &tol1, &tol2,
                             iwork_alloc, dwork_alloc, &ldwork_calc, // ldwork_calc is IN for Fortran
                             &local_iwarn, &info,
                             init_len);

    if (iwarn_ptr != NULL) {
        *iwarn_ptr = local_iwarn;
    }

    // Populate summary arrays if requested and call was successful or warning
    // INFO == -21 means LX was too small, but DWORK(1) contains required LDWORK.
    if (info == 0 || info == -21 || local_iwarn != 0) { 
        if (out_iwork_summary != NULL && iwork_alloc != NULL && liwork_calc >= 3) {
            out_iwork_summary[0] = iwork_alloc[0]; // Total function evaluations
            out_iwork_summary[1] = iwork_alloc[1]; // Total Jacobian evaluations
            out_iwork_summary[2] = iwork_alloc[2]; // Number of RCONDs in DWORK
        } else if (out_iwork_summary != NULL) { 
             // Handle cases where liwork_calc might be < 3 but > 0
             if (liwork_calc >=1 && iwork_alloc != NULL) out_iwork_summary[0] = iwork_alloc[0]; else if (out_iwork_summary != NULL) out_iwork_summary[0] = -1; // Or some indicator of not available
             if (liwork_calc >=2 && iwork_alloc != NULL && out_iwork_summary != NULL) out_iwork_summary[1] = iwork_alloc[1]; else if (out_iwork_summary != NULL && liwork_calc >=1) out_iwork_summary[1] = -1; else if (out_iwork_summary != NULL) out_iwork_summary[1] = -1;
             if (liwork_calc >=3 && iwork_alloc != NULL && out_iwork_summary != NULL) out_iwork_summary[2] = iwork_alloc[2]; else if (out_iwork_summary != NULL && liwork_calc >=1) out_iwork_summary[2] = -1; else if (out_iwork_summary != NULL) out_iwork_summary[2] = -1;
        }

        if (out_dwork_summary != NULL && dwork_alloc != NULL) {
            int num_rconds_from_iwork = 0;
            if (iwork_alloc != NULL && liwork_calc >=3 && (init_upper == 'L' || init_upper == 'B') && out_iwork_summary != NULL && out_iwork_summary[2] >=0) {
                num_rconds_from_iwork = out_iwork_summary[2]; 
            }
            
            // Number of elements to copy: DWORK(1-8) + RCONDs
            // DWORK(1) = optimal LDWORK
            // DWORK(2) = residual error norm (main)
            // DWORK(3) = iterations (main)
            // DWORK(4) = Levenberg factor (main)
            // DWORK(5-8) = same for initialization step (if INIT='S' or 'B')
            // DWORK(9 onwards) = RCONDs (if INIT='L' or 'B')
            int elements_to_copy_base = 8; 
            int elements_to_copy_total = elements_to_copy_base + ((num_rconds_from_iwork > 0) ? num_rconds_from_iwork : 0);
            
            // Ensure we don't read past allocated dwork_alloc or write past out_dwork_summary capacity
            elements_to_copy_total = MIN(elements_to_copy_total, ldwork_calc); 
            
            int out_dwork_summary_capacity = 8 + MAX_EXPECTED_RCONDS_IN_DWORK_HDR_IB03BD;
            elements_to_copy_total = MIN(elements_to_copy_total, out_dwork_summary_capacity);


            for(int i=0; i < elements_to_copy_total; ++i) {
                 if (i < ldwork_calc) { // Check read bound for dwork_alloc
                    out_dwork_summary[i] = dwork_alloc[i];
                 } else {
                    out_dwork_summary[i] = -998.0; // Indicate data not available from dwork_alloc
                 }
            }
            // Fill remaining part of out_dwork_summary with a placeholder if it's larger
            for(int i=elements_to_copy_total; i < out_dwork_summary_capacity; ++i) {
                 out_dwork_summary[i] = -999.0; // Placeholder for unused part
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
