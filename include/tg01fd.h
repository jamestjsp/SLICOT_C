/**
 * @file tg01fd.h
 * @brief Header for C wrapper of SLICOT routine TG01FD.
 * @details This routine computes for a descriptor system (A-lambda E,B,C)
 * orthogonal transformation matrices Q and Z such that the transformed system
 * (Q'*A*Z-lambda Q'*E*Z, Q'*B, C*Z) is in a SVD-like coordinate form.
 */

#ifndef SLICOT_WRAPPER_TG01FD_H
#define SLICOT_WRAPPER_TG01FD_H

#include "slicot_utils.h" // Provides SLICOT_EXPORT macro

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Computes orthogonal transformations to reduce a descriptor system
 * to a SVD-like coordinate form.
 * @details This is a C wrapper for the SLICOT Fortran routine TG01FD.
 * Workspace (IWORK, DWORK) is allocated internally.
 * Matrices A, E, B, C are modified in-place.
 * Matrices Q and Z are modified in-place if compq/compz are 'U' or 'I'.
 *
 * @param compq (Input) CHARACTER*1. Specifies computation of Q:
 * 'N': Do not compute Q. Q is not referenced by the routine. q_c can be NULL.
 * 'I': Initialize Q to identity and return the computed Q. q_c must point to allocated memory.
 * 'U': Q must contain an orthogonal matrix Q1 on entry; returns Q1*Q. q_c must point to allocated memory.
 * @param compz (Input) CHARACTER*1. Specifies computation of Z:
 * 'N': Do not compute Z. Z is not referenced by the routine. z_c can be NULL.
 * 'I': Initialize Z to identity and return the computed Z. z_c must point to allocated memory.
 * 'U': Z must contain an orthogonal matrix Z1 on entry; returns Z1*Z. z_c must point to allocated memory.
 * @param joba (Input) CHARACTER*1. Specifies reduction of A22:
 * 'N': Do not reduce A22.
 * 'R': Reduce A22 to SVD-like upper triangular form (X=0).
 * 'T': Reduce A22 to upper trapezoidal form.
 * @param l (Input) INTEGER. Number of rows of A, B, E. l >= 0.
 * @param n (Input) INTEGER. Number of columns of A, E, C. n >= 0.
 * @param m (Input) INTEGER. Number of columns of B. m >= 0.
 * @param p (Input) INTEGER. Number of rows of C. p >= 0.
 * @param a (Input/Output) DOUBLE PRECISION array, state matrix A (l x n).
 * On exit, contains transformed Q'*A*Z.
 * Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param lda (Input) INTEGER. Leading dimension of A.
 * If l > 0 and n > 0: If row_major=0, lda >= max(1,l). If row_major=1, lda >= max(1,n). Else lda >= 1.
 * @param e (Input/Output) DOUBLE PRECISION array, descriptor matrix E (l x n).
 * On exit, contains transformed Q'*E*Z.
 * Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param lde (Input) INTEGER. Leading dimension of E.
 * If l > 0 and n > 0: If row_major=0, lde >= max(1,l). If row_major=1, lde >= max(1,n). Else lde >= 1.
 * @param b (Input/Output) DOUBLE PRECISION array, input matrix B (l x m).
 * Not referenced if m=0. On exit, contains transformed Q'*B.
 * Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param ldb (Input) INTEGER. Leading dimension of B.
 * If m>0 and l>0: If row_major=0, ldb >= max(1,l). If row_major=1, ldb >= max(1,m). Else ldb >= 1.
 * @param c (Input/Output) DOUBLE PRECISION array, output matrix C (p x n).
 * Not referenced if p=0. On exit, contains transformed C*Z.
 * Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param ldc (Input) INTEGER. Leading dimension of C.
 * If p>0 and n>0: If row_major=0, ldc >= max(1,p). If row_major=1, ldc >= max(1,n). Else ldc >= 1.
 * @param q (Input/Output) DOUBLE PRECISION array, orthogonal matrix Q (l x l).
 * Referenced according to compq. If compq='N', q can be NULL.
 * Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param ldq (Input) INTEGER. Leading dimension of Q.
 * If compq='N', ldq >= 1. Else (if l>0): if row_major=0, ldq >= max(1,l); if row_major=1, ldq >= max(1,l). Else ldq >=1.
 * @param z (Input/Output) DOUBLE PRECISION array, orthogonal matrix Z (n x n).
 * Referenced according to compz. If compz='N', z can be NULL.
 * Stored column-wise if row_major=0, row-wise if row_major=1.
 * @param ldz (Input) INTEGER. Leading dimension of Z.
 * If compz='N', ldz >= 1. Else (if n>0): if row_major=0, ldz >= max(1,n); if row_major=1, ldz >= max(1,n). Else ldz >= 1.
 * @param ranke (Output) INTEGER pointer. Estimated rank of E.
 * @param rnka22 (Output) INTEGER pointer. Estimated rank of A22 if joba='R' or 'T'. Not referenced if joba='N'.
 * @param tol (Input) DOUBLE PRECISION. Tolerance for rank determination. tol < 1.
 * If tol <= 0, a default value is used.
 * @param row_major (Input) INTEGER. Matrix storage: 0 for column-major, 1 for row-major.
 *
 * @return info Error indicator:
 * = 0: successful exit.
 * < 0: if info = -i, the i-th argument had an illegal value.
 * Fortran argument indices: 1:COMPQ, 2:COMPZ, 3:JOBA, 4:L, 5:N, 6:M, 7:P,
 * 8:A, 9:LDA, 10:E, 11:LDE, 12:B, 13:LDB, 14:C, 15:LDC, 16:Q, 17:LDQ,
 * 18:Z, 19:LDZ, 20:RANKE, 21:RNKA22, 22:TOL, 23:IWORK, 24:DWORK, 25:LDWORK.
 * = SLICOT_MEMORY_ERROR (-1010): internal memory allocation failed.
 */
SLICOT_EXPORT
int slicot_tg01fd(char compq, char compz, char joba,
                  int l, int n, int m, int p,
                  double* a, int lda,
                  double* e, int lde,
                  double* b, int ldb,
                  double* c, int ldc,
                  double* q, int ldq,
                  double* z, int ldz,
                  int* ranke, int* rnka22,
                  double tol, int row_major);

#ifdef __cplusplus
}
#endif

#endif // SLICOT_WRAPPER_TG01FD_H
