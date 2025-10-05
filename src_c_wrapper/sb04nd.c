#include <stdlib.h>
#include <ctype.h>
#include <limits.h>

#include "sb04nd.h"
#include "slicot_utils.h"
#include "slicot_f77.h"

extern void F77_FUNC(sb04nd, SB04ND)(
    const char* abschu,
    const char* ula,
    const char* ulb,
    const int* n,
    const int* m,
    double* a,
    const int* lda,
    double* b,
    const int* ldb,
    double* c,
    const int* ldc,
    const double* tol,
    int* iwork,
    double* dwork,
    const int* ldwork,
    int* info,
    int abschu_len,
    int ula_len,
    int ulb_len);

SLICOT_EXPORT
int slicot_sb04nd(char abschu, char ula, char ulb,
                  int n, int m,
                  double* a, int lda,
                  double* b, int ldb,
                  double* c, int ldc,
                  double tol,
                  int row_major)
{
    int info = 0;
    int maxmn = (n > m) ? n : m;
    char abschu_u = (char)toupper((unsigned char)abschu);
    char ula_u = (char)toupper((unsigned char)ula);
    char ulb_u = (char)toupper((unsigned char)ulb);

    double dummy_double = 0.0;
    int dummy_int = 0;
    double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL;
    double *a_ptr = NULL, *b_ptr = NULL, *c_ptr = NULL;
    int lda_f = 1, ldb_f = 1, ldc_f = 1;
    int* iwork = NULL;
    double* dwork = NULL;
    int ldwork = 0;

    if (abschu_u != 'A' && abschu_u != 'B' && abschu_u != 'S') { info = -1; goto cleanup; }
    if (ula_u != 'U' && ula_u != 'L') { info = -2; goto cleanup; }
    if (ulb_u != 'U' && ulb_u != 'L') { info = -3; goto cleanup; }
    if (n < 0) { info = -4; goto cleanup; }
    if (m < 0) { info = -5; goto cleanup; }
    if (n > 0 && a == NULL) { info = -6; goto cleanup; }
    if (m > 0 && b == NULL) { info = -8; goto cleanup; }
    if (n > 0 && m > 0 && c == NULL) { info = -10; goto cleanup; }

    if (row_major) {
        if (n > 0 && lda < n) { info = -7; goto cleanup; }
        if (m > 0 && ldb < m) { info = -9; goto cleanup; }
        if (n > 0 && m > 0 && ldc < m) { info = -11; goto cleanup; }
        if (n == 0 && lda < 1) { info = -7; goto cleanup; }
        if (m == 0 && ldb < 1) { info = -9; goto cleanup; }
        if ((n == 0 || m == 0) && ldc < 1) { info = -11; goto cleanup; }
    } else {
        if (lda < ((n > 0) ? n : 1)) { info = -7; goto cleanup; }
        if (ldb < ((m > 0) ? m : 1)) { info = -9; goto cleanup; }
        if (ldc < ((n > 0) ? n : 1)) { info = -11; goto cleanup; }
    }

    if (n == 0 || m == 0) {
        return 0;
    }

    size_t a_size = (n > 0) ? (size_t)n * (size_t)n : 0;
    size_t b_size = (m > 0) ? (size_t)m * (size_t)m : 0;
    size_t c_size = (n > 0 && m > 0) ? (size_t)n * (size_t)m : 0;

    if (row_major) {
        if (a_size > 0) {
            a_cm = (double*)malloc(a_size * sizeof(double));
            CHECK_ALLOC(a_cm);
            slicot_transpose_to_fortran_with_ld(a, a_cm, n, n, lda, n > 0 ? n : 1, sizeof(double));
            a_ptr = a_cm;
            lda_f = (n > 0) ? n : 1;
        } else {
            lda_f = 1;
            a_ptr = &dummy_double;
        }

        if (b_size > 0) {
            b_cm = (double*)malloc(b_size * sizeof(double));
            CHECK_ALLOC(b_cm);
            slicot_transpose_to_fortran_with_ld(b, b_cm, m, m, ldb, m > 0 ? m : 1, sizeof(double));
            b_ptr = b_cm;
            ldb_f = (m > 0) ? m : 1;
        } else {
            ldb_f = 1;
            b_ptr = &dummy_double;
        }

        if (c_size > 0) {
            c_cm = (double*)malloc(c_size * sizeof(double));
            CHECK_ALLOC(c_cm);
            slicot_transpose_to_fortran_with_ld(c, c_cm, n, m, ldc, n > 0 ? n : 1, sizeof(double));
            c_ptr = c_cm;
            ldc_f = (n > 0) ? n : 1;
        } else {
            ldc_f = (n > 0) ? n : 1;
            c_ptr = &dummy_double;
        }
    } else {
        lda_f = (n > 0) ? ((lda > 0) ? lda : n) : 1;
        ldb_f = (m > 0) ? ((ldb > 0) ? ldb : m) : 1;
        ldc_f = (n > 0) ? ((ldc > 0) ? ldc : n) : 1;
        a_ptr = (n > 0) ? a : (a ? a : &dummy_double);
        b_ptr = (m > 0) ? b : (b ? b : &dummy_double);
        c_ptr = (n > 0 && m > 0) ? c : ((c && (n > 0 || m > 0)) ? c : &dummy_double);
    }

    int need_work = !(abschu_u == 'S' && ula_u == 'U' && ulb_u == 'U');

    if (need_work && maxmn > 0) {
        long long maxmn_ll = (long long)maxmn;
        long long iwork_len_ll = 2LL * maxmn_ll;
        if (iwork_len_ll < 1) iwork_len_ll = 1;
        if (iwork_len_ll > INT_MAX) { info = SLICOT_MEMORY_ERROR; goto cleanup; }
        int iwork_len = (int)iwork_len_ll;
        iwork = (int*)malloc((size_t)iwork_len * sizeof(int));
        CHECK_ALLOC(iwork);

        long long dwork_len_ll = 2LL * maxmn_ll * (4LL + 2LL * maxmn_ll);
        if (dwork_len_ll < 1) dwork_len_ll = 1;
        if (dwork_len_ll > INT_MAX) { info = SLICOT_MEMORY_ERROR; goto cleanup; }
        ldwork = (int)dwork_len_ll;
        dwork = (double*)malloc((size_t)ldwork * sizeof(double));
        CHECK_ALLOC(dwork);
    }

    {
        int abschu_len = 1, ula_len = 1, ulb_len = 1;
        int* iwork_ptr = iwork ? iwork : &dummy_int;
        double* dwork_ptr = dwork ? dwork : &dummy_double;
        F77_FUNC(sb04nd, SB04ND)(
            &abschu_u,
            &ula_u,
            &ulb_u,
            &n,
            &m,
            a_ptr,
            &lda_f,
            b_ptr,
            &ldb_f,
            c_ptr,
            &ldc_f,
            &tol,
            iwork_ptr,
            dwork_ptr,
            &ldwork,
            &info,
            abschu_len,
            ula_len,
            ulb_len);
    }

    if (row_major && info >= 0) {
        if (a_cm && a) {
            slicot_transpose_to_c_with_ld(a_cm, a, n, n, lda_f, lda, sizeof(double));
        }
        if (b_cm && b) {
            slicot_transpose_to_c_with_ld(b_cm, b, m, m, ldb_f, ldb, sizeof(double));
        }
        if (c_cm && c) {
            slicot_transpose_to_c_with_ld(c_cm, c, n, m, ldc_f, ldc, sizeof(double));
        }
    }

cleanup:
    free(a_cm);
    free(b_cm);
    free(c_cm);
    free(iwork);
    free(dwork);
    return info;
}
