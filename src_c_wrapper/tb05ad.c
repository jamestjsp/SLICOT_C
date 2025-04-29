/**
 * @file tb05ad.c
 * @brief C wrapper for SLICOT routine TB05AD.
 * @details Finds the complex frequency response matrix (transfer matrix)
 * of a given state-space representation (A,B,C).
 * Workspace (IWORK, DWORK, ZWORK) is allocated internally by this wrapper.
 * Input/output matrix format is handled via the row_major parameter.
 */

#include <stdlib.h> // For malloc, free  
#include <ctype.h>  // For toupper  
#include <stddef.h> // For size_t  
#include <stdio.h>  // For error logging

#include "tb05ad.h"     // Public header for this wrapper  
#include "slicot_utils.h" // Provides CHECK_ALLOC, SLICOT_MEMORY_ERROR, MAX/MIN, transpose functions
#include "slicot_f77.h"   // Provides F77_FUNC macro for Fortran name mangling

/* External Fortran routine declaration */
extern void F77_FUNC(tb05ad, TB05AD)(
    const char* baleig, const char* inita,
    const int* n, const int* m, const int* p,
    const void* freq, /* Complex*16 */
    double* a, const int* lda,
    double* b, const int* ldb,
    double* c, const int* ldc,
    double* rcond,
    void* g, /* Complex*16 */
    const int* ldg,
    double* evre, double* evim,
    void* hinvb, /* Complex*16 */
    const int* ldhinv,
    int* iwork, double* dwork, const int* ldwork,
    void* zwork, /* Complex*16 */
    const int* lzwork,
    int* info,
    int baleig_len, int inita_len
);

/* Utility function to create complex*16 value from real and imaginary parts */
static void make_complex(double re, double im, void* complex_ptr) {
    double* ptr = (double*)complex_ptr;
    ptr[0] = re;
    ptr[1] = im;
}

/* C wrapper function definition */
SLICOT_EXPORT
int slicot_tb05ad(char baleig, char inita, int n, int m, int p,
                 double freq_re, double freq_im,
                 double* a, int lda, double* b, int ldb, double* c, int ldc,
                 double* rcond, double* g, int ldg,
                 double* evre, double* evim,
                 double* hinvb, int ldhinv, int row_major)
{
    // Variable declarations
    int info = 0;          // Return status code  
    int *iwork = NULL;     // Internal integer workspace pointer  
    double *dwork = NULL;  // Internal double workspace pointer
    void *zwork = NULL;    // Internal complex workspace pointer
    int ldwork = 0;        // Calculated size for dwork
    int lzwork = 0;        // Calculated size for zwork
    
    // Complex frequency value
    double complex_freq[2]; // Directly define an array for freq
    
    // Pointers for column-major copies (if row_major is used)
    double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL;
    double *g_cm = NULL, *hinvb_cm = NULL;
    
    // Convert character options to uppercase once
    char baleig_upper = toupper(baleig);
    char inita_upper = toupper(inita);

    // Input parameter validation
    // Check scalar parameters first
    if (n < 0) { info = -3; goto cleanup; }
    if (m < 0) { info = -4; goto cleanup; }
    if (p < 0) { info = -5; goto cleanup; }
    
    // Check character options
    if (baleig_upper != 'N' && baleig_upper != 'C' && 
        baleig_upper != 'B' && baleig_upper != 'E' && baleig_upper != 'A') { 
        info = -1; goto cleanup; 
    }
    if (inita_upper != 'G' && inita_upper != 'H') { 
        info = -2; goto cleanup; 
    }

    // Check pointers for required arrays
    if (a == NULL && n > 0) { info = -8; goto cleanup; }
    if (b == NULL && n > 0 && m > 0) { info = -10; goto cleanup; }
    if (c == NULL && p > 0 && n > 0) { info = -12; goto cleanup; }
    if (g == NULL && p > 0 && m > 0) { info = -14; goto cleanup; }
    
    // Check eigenvalue arrays if needed
    if ((inita_upper == 'G') && (baleig_upper == 'B' || baleig_upper == 'E' || baleig_upper == 'A')) {
        if (evre == NULL && n > 0) { info = -16; goto cleanup; }
        if (evim == NULL && n > 0) { info = -17; goto cleanup; }
    }
    
    if (hinvb == NULL && n > 0 && m > 0) { info = -18; goto cleanup; }

    // Check leading dimensions based on row_major flag
    int min_lda_f = MAX(1, n); // Minimum Fortran LDA (rows) for A
    int min_ldb_f = MAX(1, n); // Minimum Fortran LDB (rows) for B
    int min_ldc_f = MAX(1, p); // Minimum Fortran LDC (rows) for C
    int min_ldg_f = MAX(1, p); // Minimum Fortran LDG (rows) for G
    int min_ldhinv_f = MAX(1, n); // Minimum Fortran LDHINV (rows) for HINVB
    
    if (row_major) {
        // C LDs are number of columns
        if (n > 0 && lda < n) { info = -9; goto cleanup; }
        if (m > 0 && ldb < m) { info = -11; goto cleanup; }
        if (n > 0 && ldc < n) { info = -13; goto cleanup; }
        if (m > 0 && ldg < m) { info = -15; goto cleanup; }
        if (m > 0 && ldhinv < m) { info = -19; goto cleanup; }
    } else {
        // C LDs are number of rows (Fortran style)
        if (lda < min_lda_f) { info = -9; goto cleanup; }
        if (ldb < min_ldb_f) { info = -11; goto cleanup; }
        if (ldc < min_ldc_f) { info = -13; goto cleanup; }
        if (ldg < min_ldg_f) { info = -15; goto cleanup; }
        if (ldhinv < min_ldhinv_f) { info = -19; goto cleanup; }
    }

    // Exit if any validation failed before allocating memory
    if (info != 0) { goto cleanup; }

    // Create the complex frequency value
    complex_freq[0] = freq_re;
    complex_freq[1] = freq_im;
    
    // Internal Workspace Allocation
    iwork = (int*)malloc((size_t)n * sizeof(int));
    CHECK_ALLOC(iwork);
    
    // Calculate workspace sizes based on documentation formulas
    if (inita_upper == 'G') {
        if (baleig_upper == 'N' || baleig_upper == 'B' || baleig_upper == 'E') {
            ldwork = MAX(1, n - 1 + MAX(n, MAX(m, p)));
        } else if (baleig_upper == 'C' || baleig_upper == 'A') {
            ldwork = MAX(1, n + MAX(n, MAX(m-1, p-1)));
        }
    } else if (inita_upper == 'H' && (baleig_upper == 'C' || baleig_upper == 'A')) {
        ldwork = MAX(1, 2*n);
    } else {
        ldwork = 1;
    }
    
    dwork = (double*)malloc((size_t)ldwork * sizeof(double));
    CHECK_ALLOC(dwork);
    
    // ZWORK size
    if (baleig_upper == 'C' || baleig_upper == 'A') {
        lzwork = MAX(1, n*n + 2*n);
    } else {
        lzwork = MAX(1, n*n);
    }
    
    zwork = malloc((size_t)lzwork * 2 * sizeof(double)); // Complex*16 = 2 doubles each
    CHECK_ALLOC(zwork);
    
    // Calculate sizes for matrix arrays
    size_t a_size = (size_t)n * n; if (n == 0) a_size = 0;
    size_t b_size = (size_t)n * m; if (n == 0 || m == 0) b_size = 0;
    size_t c_size = (size_t)p * n; if (p == 0 || n == 0) c_size = 0;
    size_t g_size = (size_t)p * m; if (p == 0 || m == 0) g_size = 0;
    size_t hinvb_size = (size_t)n * m; if (n == 0 || m == 0) hinvb_size = 0;
    
    // Allocate memory for column-major copies if row_major
    if (row_major) {
        if (a_size > 0) {
            a_cm = (double*)malloc(a_size * sizeof(double));
            CHECK_ALLOC(a_cm);
        }
        if (b_size > 0) {
            b_cm = (double*)malloc(b_size * sizeof(double));
            CHECK_ALLOC(b_cm);
        }
        if (c_size > 0) {
            c_cm = (double*)malloc(c_size * sizeof(double));
            CHECK_ALLOC(c_cm);
        }
        if (g_size > 0) {
            g_cm = (double*)malloc(2 * g_size * sizeof(double)); // Complex values are twice as large
            CHECK_ALLOC(g_cm);
        }
        if (hinvb_size > 0) {
            hinvb_cm = (double*)malloc(2 * hinvb_size * sizeof(double)); // Complex values are twice as large
            CHECK_ALLOC(hinvb_cm);
        }
    }
    
    // Prepare Fortran parameters and perform conversions
    double* a_ptr = a;
    double* b_ptr = b;
    double* c_ptr = c;
    double* g_ptr = g;
    double* hinvb_ptr = hinvb;
    
    int lda_f = lda;
    int ldb_f = ldb;
    int ldc_f = ldc;
    int ldg_f = ldg;
    int ldhinv_f = ldhinv;
    
    if (row_major) {
        // Set Fortran LDs (number of ROWS)
        lda_f = MAX(1, n);
        ldb_f = MAX(1, n);
        ldc_f = MAX(1, p);
        ldg_f = MAX(1, p);
        ldhinv_f = MAX(1, n);
        
        // Point to column-major copies
        if (a_size > 0) {
            slicot_transpose_to_fortran(a, a_cm, n, n, sizeof(double));
            a_ptr = a_cm;
        } else {
            a_ptr = NULL;
        }
        
        if (b_size > 0) {
            slicot_transpose_to_fortran(b, b_cm, n, m, sizeof(double));
            b_ptr = b_cm;
        } else {
            b_ptr = NULL;
        }
        
        if (c_size > 0) {
            slicot_transpose_to_fortran(c, c_cm, p, n, sizeof(double));
            c_ptr = c_cm;
        } else {
            c_ptr = NULL;
        }
        
        // Output arrays - just use the column-major buffer
        g_ptr = (g_size > 0) ? g_cm : NULL;
        hinvb_ptr = (hinvb_size > 0) ? hinvb_cm : NULL;
        
    } else {
        // Column-major: Ensure NULL is passed if size is 0
        if (a_size == 0) a_ptr = NULL;
        if (b_size == 0) b_ptr = NULL;
        if (c_size == 0) c_ptr = NULL;
        if (g_size == 0) g_ptr = NULL;
        if (hinvb_size == 0) hinvb_ptr = NULL;
    }
    
    // Hidden string lengths
    const int baleig_len = 1;
    const int inita_len = 1;
    
    // Call Fortran function
    F77_FUNC(tb05ad, TB05AD)(
        &baleig_upper, &inita_upper,
        &n, &m, &p,
        complex_freq,
        a_ptr, &lda_f,
        b_ptr, &ldb_f,
        c_ptr, &ldc_f,
        rcond,
        g_ptr, &ldg_f,
        evre, evim,
        hinvb_ptr, &ldhinv_f,
        iwork, dwork, &ldwork,
        zwork, &lzwork,
        &info,
        baleig_len, inita_len
    );
    
    // Convert results back to row-major (if needed)
    if (row_major && info == 0) {
        // Convert modified inputs back if needed
        if (inita_upper == 'G' && a_size > 0) {
            slicot_transpose_to_c(a_cm, a, n, n, sizeof(double));
        }
        if (inita_upper == 'G' && b_size > 0) {
            slicot_transpose_to_c(b_cm, b, n, m, sizeof(double));
        }
        if (inita_upper == 'G' && c_size > 0) {
            slicot_transpose_to_c(c_cm, c, p, n, sizeof(double));
        }
        
        // Convert complex outputs with special handling for complex values
        // G is a complex p-by-m array
        if (g_size > 0) {
            // For complex values we need special handling since transposing interleaved complex values
            // is different from transposing real values. For now we can use standard transpose for interleaved values.
            slicot_transpose_to_c(g_cm, g, p, m, 2 * sizeof(double));
        }
        
        // HINVB is a complex n-by-m array
        if (hinvb_size > 0) {
            slicot_transpose_to_c(hinvb_cm, hinvb, n, m, 2 * sizeof(double));
        }
    }
    
cleanup:
    /* --- Cleanup --- */
    // Free allocated memory
    free(iwork);
    free(dwork);
    free(zwork);
    free(a_cm);
    free(b_cm);
    free(c_cm);
    free(g_cm);
    free(hinvb_cm);

    // Check if info was set by CHECK_ALLOC during workspace/copy allocation
    if (info == SLICOT_MEMORY_ERROR) {
       fprintf(stderr, "Error: Memory allocation failed in slicot_tb05ad.\n");
    }
    
    // Return the info code
    return info;
}