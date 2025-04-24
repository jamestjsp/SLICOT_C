/**
 * @file ab05qd.h
 * @brief C wrapper for SLICOT routine AB05QD
 *
 * This file provides a C interface to the SLICOT routine AB05QD,
 * which computes the state-space model (A,B,C,D) corresponding to
 * appending two systems G1 and G2, resulting in G = diag(G1, G2).
 */

 #ifndef AB05QD_H
 #define AB05QD_H
 
 #include <stddef.h> // For size_t
 
 #include "slicot_utils.h" 

 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Computes the appended system G = diag(G1, G2).
  *
  * Calculates the state-space representation (A,B,C,D) of the system
  * formed by appending System 1 (A1,B1,C1,D1) and System 2 (A2,B2,C2,D2)
  * such that the resulting transfer function is block diagonal diag(G1, G2).
  *
  * Note: This wrapper assumes OVER = 'N' (no overlapping arrays).
  *
  * @param[in] over      Specifies array overlapping (Not supported by wrapper, assumed 'N'):
  * = 'N': Do not overlap arrays.
  * = 'O': Overlap pairs of arrays (A1/A, B1/B, C1/C, D1/D).
  * @param[in] n1        Order of system 1 (matrix A1), n1 >= 0.
  * @param[in] m1        Number of inputs for system 1 (columns of B1, D1), m1 >= 0.
  * @param[in] p1        Number of outputs for system 1 (rows of C1, D1), p1 >= 0.
  * @param[in] n2        Order of system 2 (matrix A2), n2 >= 0.
  * @param[in] m2        Number of inputs for system 2 (columns of B2, D2), m2 >= 0.
  * @param[in] p2        Number of outputs for system 2 (rows of C2, D2), p2 >= 0.
  * @param[in] a1        Double array for matrix A1 (System 1). Dimension (lda1, n1) or (n1, lda1).
  * @param[in] lda1      Leading dimension of A1. >= max(1,n1).
  * @param[in] b1        Double array for matrix B1 (System 1). Dimension (ldb1, m1) or (n1, ldb1).
  * @param[in] ldb1      Leading dimension of B1. >= max(1,n1).
  * @param[in] c1        Double array for matrix C1 (System 1). Dimension (ldc1, n1) or (p1, ldc1).
  * @param[in] ldc1      Leading dimension of C1. >= max(1,p1) if n1>0, else >=1.
  * @param[in] d1        Double array for matrix D1 (System 1). Dimension (ldd1, m1) or (p1, ldd1).
  * @param[in] ldd1      Leading dimension of D1. >= max(1,p1).
  * @param[in] a2        Double array for matrix A2 (System 2). Dimension (lda2, n2) or (n2, lda2).
  * @param[in] lda2      Leading dimension of A2. >= max(1,n2).
  * @param[in] b2        Double array for matrix B2 (System 2). Dimension (ldb2, m2) or (n2, ldb2).
  * @param[in] ldb2      Leading dimension of B2. >= max(1,n2).
  * @param[in] c2        Double array for matrix C2 (System 2). Dimension (ldc2, n2) or (p2, ldc2).
  * @param[in] ldc2      Leading dimension of C2. >= max(1,p2) if n2>0, else >=1.
  * @param[in] d2        Double array for matrix D2 (System 2). Dimension (ldd2, m2) or (p2, ldd2).
  * @param[in] ldd2      Leading dimension of D2. >= max(1,p2).
  * @param[out] n        Pointer to integer for the order of the resulting system A (n1 + n2).
  * @param[out] m        Pointer to integer for the number of inputs of the resulting system B (m1 + m2).
  * @param[out] p        Pointer to integer for the number of outputs of the resulting system C (p1 + p2).
  * @param[out] a        Double array for the resulting matrix A. Dimension (lda, n) or (n, lda).
  * @param[in] lda       Leading dimension of A. >= max(1, n1+n2).
  * @param[out] b        Double array for the resulting matrix B. Dimension (ldb, m) or (n, ldb).
  * @param[in] ldb       Leading dimension of B. >= max(1, n1+n2).
  * @param[out] c        Double array for the resulting matrix C. Dimension (ldc, n) or (p, ldc).
  * @param[in] ldc       Leading dimension of C. >= max(1,p1+p2) if n1+n2>0, else >=1.
  * @param[out] d        Double array for the resulting matrix D. Dimension (ldd, m) or (p, ldd).
  * @param[in] ldd       Leading dimension of D. >= max(1,p1+p2).
  * @param[in] row_major Integer flag:
  * = 0: All arrays are column-major (Fortran style).
  * = 1: All arrays are row-major (C style).
  *
  * @return info         Error indicator:
  * = 0: successful exit
  * < 0: if info = -i, the i-th argument had an illegal value.
  * Memory allocation errors may also be returned.
  */
 SLICOT_EXPORT
 int slicot_ab05qd(char over,
                   int n1, int m1, int p1, int n2, int m2, int p2,
                   const double* a1, int lda1, const double* b1, int ldb1,
                   const double* c1, int ldc1, const double* d1, int ldd1,
                   const double* a2, int lda2, const double* b2, int ldb2,
                   const double* c2, int ldc2, const double* d2, int ldd2,
                   int* n, int* m, int* p,
                   double* a, int lda, double* b, int ldb,
                   double* c, int ldc, double* d, int ldd,
                   int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* AB05QD_H */
 