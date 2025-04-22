#include <gtest/gtest.h>
#include <vector>
#include <cmath>

#include "ab05qd.h"

// Column-major test fixture
class Ab05qdTestColMajor : public ::testing::Test {
protected:
    // Test parameters
    int N1 = 3;
    int M1 = 2;
    int P1 = 2;
    int N2 = 3;
    int M2 = 2;
    int P2 = 2;
    int N, M, P;
    char OVER = 'N';
    int ROW_MAJOR = 0;  // Column-major
    double check_tol = 1e-12;
    
    // Input matrices (column-major order)
    // A1 (3×3) - First system state matrix
    std::vector<double> A1 = {
        /* Column 1 */ 1.0, 0.0, 1.0,
        /* Column 2 */ 0.0, -1.0, 1.0,
        /* Column 3 */ -1.0, 1.0, 2.0
    };
    
    // B1 (3×2) - First system input matrix
    std::vector<double> B1 = {
        /* Column 1 */ 1.0, 1.0, 0.0,
        /* Column 2 */ 2.0, 0.0, 1.0
    };
    
    // C1 (2×3) - First system output matrix
    std::vector<double> C1 = {
        /* Column 1 */ 3.0, 0.0,
        /* Column 2 */ -2.0, 1.0,
        /* Column 3 */ 1.0, 0.0
    };
    
    // D1 (2×2) - First system feedthrough matrix
    std::vector<double> D1 = {
        /* Column 1 */ 1.0, 0.0,
        /* Column 2 */ 0.0, 1.0
    };
    
    // A2 (3×3) - Second system state matrix
    std::vector<double> A2 = {
        /* Column 1 */ -3.0, 1.0, 0.0,
        /* Column 2 */ 0.0, 0.0, -1.0,
        /* Column 3 */ 0.0, 1.0, 2.0
    };
    
    // B2 (3×2) - Second system input matrix
    std::vector<double> B2 = {
        /* Column 1 */ 0.0, -1.0, 0.0,
        /* Column 2 */ 1.0, 0.0, 2.0
    };
    
    // C2 (2×3) - Second system output matrix
    std::vector<double> C2 = {
        /* Column 1 */ 1.0, 1.0,
        /* Column 2 */ 1.0, 1.0,
        /* Column 3 */ 0.0, -1.0
    };
    
    // D2 (2×2) - Second system feedthrough matrix
    std::vector<double> D2 = {
        /* Column 1 */ 1.0, 0.0,
        /* Column 2 */ 1.0, 1.0
    };
    
    // Expected outputs for combined system (column-major order)
    // A (6×6) - Combined system state matrix
    std::vector<double> A_expected = {
        /* Column 1 */ 1.0, 0.0, 1.0, 0.0, 0.0, 0.0,
        /* Column 2 */ 0.0, -1.0, 1.0, 0.0, 0.0, 0.0,
        /* Column 3 */ -1.0, 1.0, 2.0, 0.0, 0.0, 0.0,
        /* Column 4 */ 0.0, 0.0, 0.0, -3.0, 1.0, 0.0,
        /* Column 5 */ 0.0, 0.0, 0.0, 0.0, 0.0, -1.0,
        /* Column 6 */ 0.0, 0.0, 0.0, 0.0, 1.0, 2.0
    };
    
    // B (6×4) - Combined system input matrix
    std::vector<double> B_expected = {
        /* Column 1 */ 1.0, 1.0, 0.0, 0.0, 0.0, 0.0,
        /* Column 2 */ 2.0, 0.0, 1.0, 0.0, 0.0, 0.0,
        /* Column 3 */ 0.0, 0.0, 0.0, 0.0, -1.0, 0.0,
        /* Column 4 */ 0.0, 0.0, 0.0, 1.0, 0.0, 2.0
    };
    
    // C (4×6) - Combined system output matrix
    std::vector<double> C_expected = {
        /* Column 1 */ 3.0, 0.0, 0.0, 0.0,
        /* Column 2 */ -2.0, 1.0, 0.0, 0.0,
        /* Column 3 */ 1.0, 0.0, 0.0, 0.0,
        /* Column 4 */ 0.0, 0.0, 1.0, 1.0,
        /* Column 5 */ 0.0, 0.0, 1.0, 1.0,
        /* Column 6 */ 0.0, 0.0, 0.0, -1.0
    };
    
    // D (4×4) - Combined system feedthrough matrix
    std::vector<double> D_expected = {
        /* Column 1 */ 1.0, 0.0, 0.0, 0.0,
        /* Column 2 */ 0.0, 1.0, 0.0, 0.0,
        /* Column 3 */ 0.0, 0.0, 1.0, 0.0,
        /* Column 4 */ 0.0, 0.0, 1.0, 1.0
    };
    
    int expected_info = 0;
    
    void SetUp() override {
        // Calculate expected dimensions for the combined system
        N = N1 + N2;
        M = M1 + M2;
        P = P1 + P2;
    }
};

// Row-major test fixture
class Ab05qdTestRowMajor : public Ab05qdTestColMajor {
public:    Ab05qdTestRowMajor() {
        // Make sure N, M, and P are set since they're used for resizing vectors
        N = N1 + N2;
        M = M1 + M2;
        P = P1 + P2;
        
        ROW_MAJOR = 1;  // Row-major
        
        // Convert input matrices to row-major
        A1_rm.resize(N1 * N1);
        B1_rm.resize(N1 * M1);
        C1_rm.resize(P1 * N1);
        D1_rm.resize(P1 * M1);
        A2_rm.resize(N2 * N2);
        B2_rm.resize(N2 * M2);
        C2_rm.resize(P2 * N2);
        D2_rm.resize(P2 * M2);
        
        // Convert expected output matrices to row-major
        A_expected_rm.resize(N * N);
        B_expected_rm.resize(N * M);
        C_expected_rm.resize(P * N);
        D_expected_rm.resize(P * M);
        
        // Convert A1 from column-major to row-major
        for (int i = 0; i < N1; ++i) {
            for (int j = 0; j < N1; ++j) {
                A1_rm[i*N1 + j] = A1[i + j*N1];
            }
        }
        
        // Convert B1 from column-major to row-major
        for (int i = 0; i < N1; ++i) {
            for (int j = 0; j < M1; ++j) {
                B1_rm[i*M1 + j] = B1[i + j*N1];
            }
        }
        
        // Convert C1 from column-major to row-major
        for (int i = 0; i < P1; ++i) {
            for (int j = 0; j < N1; ++j) {
                C1_rm[i*N1 + j] = C1[i + j*P1];
            }
        }
        
        // Convert D1 from column-major to row-major
        for (int i = 0; i < P1; ++i) {
            for (int j = 0; j < M1; ++j) {
                D1_rm[i*M1 + j] = D1[i + j*P1];
            }
        }
        
        // Convert A2 from column-major to row-major
        for (int i = 0; i < N2; ++i) {
            for (int j = 0; j < N2; ++j) {
                A2_rm[i*N2 + j] = A2[i + j*N2];
            }
        }
        
        // Convert B2 from column-major to row-major
        for (int i = 0; i < N2; ++i) {
            for (int j = 0; j < M2; ++j) {
                B2_rm[i*M2 + j] = B2[i + j*N2];
            }
        }
        
        // Convert C2 from column-major to row-major
        for (int i = 0; i < P2; ++i) {
            for (int j = 0; j < N2; ++j) {
                C2_rm[i*N2 + j] = C2[i + j*P2];
            }
        }
        
        // Convert D2 from column-major to row-major
        for (int i = 0; i < P2; ++i) {
            for (int j = 0; j < M2; ++j) {
                D2_rm[i*M2 + j] = D2[i + j*P2];
            }
        }
        
        // Convert A_expected from column-major to row-major
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                A_expected_rm[i*N + j] = A_expected[i + j*N];
            }
        }
        
        // Convert B_expected from column-major to row-major
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                B_expected_rm[i*M + j] = B_expected[i + j*N];
            }
        }
        
        // Convert C_expected from column-major to row-major
        for (int i = 0; i < P; ++i) {
            for (int j = 0; j < N; ++j) {
                C_expected_rm[i*N + j] = C_expected[i + j*P];
            }
        }
        
        // Convert D_expected from column-major to row-major
        for (int i = 0; i < P; ++i) {
            for (int j = 0; j < M; ++j) {
                D_expected_rm[i*M + j] = D_expected[i + j*P];
            }
        }
    }
    
    std::vector<double> A1_rm, B1_rm, C1_rm, D1_rm;
    std::vector<double> A2_rm, B2_rm, C2_rm, D2_rm;
    std::vector<double> A_expected_rm, B_expected_rm, C_expected_rm, D_expected_rm;
};

// Column-major test
TEST_F(Ab05qdTestColMajor, DocExample) {
    // Set up
    int lda1 = N1, ldb1 = N1, ldc1 = P1, ldd1 = P1;
    int lda2 = N2, ldb2 = N2, ldc2 = P2, ldd2 = P2;
    int lda = N, ldb = N, ldc = P, ldd = P;
    
    std::vector<double> A(N * N);
    std::vector<double> B(N * M);
    std::vector<double> C(P * N);
    std::vector<double> D(P * M);
    
    int n_out, m_out, p_out;
    
    // Call function
    int info = slicot_ab05qd(OVER, N1, M1, P1, N2, M2, P2,
                             A1.data(), lda1, B1.data(), ldb1,
                             C1.data(), ldc1, D1.data(), ldd1,
                             A2.data(), lda2, B2.data(), ldb2,
                             C2.data(), ldc2, D2.data(), ldd2,
                             &n_out, &m_out, &p_out,
                             A.data(), lda, B.data(), ldb,
                             C.data(), ldc, D.data(), ldd,
                             ROW_MAJOR);
    
    // Verify
    ASSERT_EQ(info, expected_info);
    EXPECT_EQ(n_out, N);
    EXPECT_EQ(m_out, M);
    EXPECT_EQ(p_out, P);
    
    // Verify A matrix
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            EXPECT_NEAR(A[i + j*lda], A_expected[i + j*N], check_tol)
                << "A[" << i << "," << j << "] mismatch";
        }
    }
    
    // Verify B matrix
    for (int j = 0; j < M; ++j) {
        for (int i = 0; i < N; ++i) {
            EXPECT_NEAR(B[i + j*ldb], B_expected[i + j*N], check_tol)
                << "B[" << i << "," << j << "] mismatch";
        }
    }
    
    // Verify C matrix
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < P; ++i) {
            EXPECT_NEAR(C[i + j*ldc], C_expected[i + j*P], check_tol)
                << "C[" << i << "," << j << "] mismatch";
        }
    }
    
    // Verify D matrix
    for (int j = 0; j < M; ++j) {
        for (int i = 0; i < P; ++i) {
            EXPECT_NEAR(D[i + j*ldd], D_expected[i + j*P], check_tol)
                << "D[" << i << "," << j << "] mismatch";
        }
    }
}

// Row-major test
TEST_F(Ab05qdTestRowMajor, DocExample) {
    // Set up
    int lda1 = N1, ldb1 = M1, ldc1 = N1, ldd1 = M1;
    int lda2 = N2, ldb2 = M2, ldc2 = N2, ldd2 = M2;
    int lda = N, ldb = M, ldc = N, ldd = M;
    
    std::vector<double> A(N * N);
    std::vector<double> B(N * M);
    std::vector<double> C(P * N);
    std::vector<double> D(P * M);
    
    int n_out, m_out, p_out;
    
    // Call function
    int info = slicot_ab05qd(OVER, N1, M1, P1, N2, M2, P2,
                             A1_rm.data(), lda1, B1_rm.data(), ldb1,
                             C1_rm.data(), ldc1, D1_rm.data(), ldd1,
                             A2_rm.data(), lda2, B2_rm.data(), ldb2,
                             C2_rm.data(), ldc2, D2_rm.data(), ldd2,
                             &n_out, &m_out, &p_out,
                             A.data(), lda, B.data(), ldb,
                             C.data(), ldc, D.data(), ldd,
                             ROW_MAJOR);
    
    // Verify
    ASSERT_EQ(info, expected_info);
    EXPECT_EQ(n_out, N);
    EXPECT_EQ(m_out, M);
    EXPECT_EQ(p_out, P);
    
    // Verify A matrix
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            EXPECT_NEAR(A[i*lda + j], A_expected_rm[i*N + j], check_tol)
                << "A[" << i << "," << j << "] mismatch";
        }
    }
    
    // Verify B matrix
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            EXPECT_NEAR(B[i*ldb + j], B_expected_rm[i*M + j], check_tol)
                << "B[" << i << "," << j << "] mismatch";
        }
    }
    
    // Verify C matrix
    for (int i = 0; i < P; ++i) {
        for (int j = 0; j < N; ++j) {
            EXPECT_NEAR(C[i*ldc + j], C_expected_rm[i*N + j], check_tol)
                << "C[" << i << "," << j << "] mismatch";
        }
    }
    
    // Verify D matrix
    for (int i = 0; i < P; ++i) {
        for (int j = 0; j < M; ++j) {
            EXPECT_NEAR(D[i*ldd + j], D_expected_rm[i*M + j], check_tol)
                << "D[" << i << "," << j << "] mismatch";
        }
    }
}

// Input validation test cases
TEST_F(Ab05qdTestColMajor, InputValidation) {
    int lda1 = N1, ldb1 = N1, ldc1 = P1, ldd1 = P1;
    int lda2 = N2, ldb2 = N2, ldc2 = P2, ldd2 = P2;
    int lda = N, ldb = N, ldc = P, ldd = P;
    
    std::vector<double> A(N * N);
    std::vector<double> B(N * M);
    std::vector<double> C(P * N);
    std::vector<double> D(P * M);
    
    int n_out, m_out, p_out;
    int info;
    
    // Test invalid N1
    info = slicot_ab05qd(OVER, -1, M1, P1, N2, M2, P2,
                         A1.data(), lda1, B1.data(), ldb1,
                         C1.data(), ldc1, D1.data(), ldd1,
                         A2.data(), lda2, B2.data(), ldb2,
                         C2.data(), ldc2, D2.data(), ldd2,
                         &n_out, &m_out, &p_out,
                         A.data(), lda, B.data(), ldb,
                         C.data(), ldc, D.data(), ldd,
                         ROW_MAJOR);
    EXPECT_LT(info, 0);
    
    // Test invalid M1
    info = slicot_ab05qd(OVER, N1, -1, P1, N2, M2, P2,
                         A1.data(), lda1, B1.data(), ldb1,
                         C1.data(), ldc1, D1.data(), ldd1,
                         A2.data(), lda2, B2.data(), ldb2,
                         C2.data(), ldc2, D2.data(), ldd2,
                         &n_out, &m_out, &p_out,
                         A.data(), lda, B.data(), ldb,
                         C.data(), ldc, D.data(), ldd,
                         ROW_MAJOR);
    EXPECT_LT(info, 0);
    
    // Test invalid P1
    info = slicot_ab05qd(OVER, N1, M1, -1, N2, M2, P2,
                         A1.data(), lda1, B1.data(), ldb1,
                         C1.data(), ldc1, D1.data(), ldd1,
                         A2.data(), lda2, B2.data(), ldb2,
                         C2.data(), ldc2, D2.data(), ldd2,
                         &n_out, &m_out, &p_out,
                         A.data(), lda, B.data(), ldb,
                         C.data(), ldc, D.data(), ldd,
                         ROW_MAJOR);
    EXPECT_LT(info, 0);
    
    // Test invalid N2
    info = slicot_ab05qd(OVER, N1, M1, P1, -1, M2, P2,
                         A1.data(), lda1, B1.data(), ldb1,
                         C1.data(), ldc1, D1.data(), ldd1,
                         A2.data(), lda2, B2.data(), ldb2,
                         C2.data(), ldc2, D2.data(), ldd2,
                         &n_out, &m_out, &p_out,
                         A.data(), lda, B.data(), ldb,
                         C.data(), ldc, D.data(), ldd,
                         ROW_MAJOR);
    EXPECT_LT(info, 0);
    
    // Test invalid leading dimensions
    info = slicot_ab05qd(OVER, N1, M1, P1, N2, M2, P2,
                         A1.data(), 0, B1.data(), ldb1,  // lda1 = 0
                         C1.data(), ldc1, D1.data(), ldd1,
                         A2.data(), lda2, B2.data(), ldb2,
                         C2.data(), ldc2, D2.data(), ldd2,
                         &n_out, &m_out, &p_out,
                         A.data(), lda, B.data(), ldb,
                         C.data(), ldc, D.data(), ldd,
                         ROW_MAJOR);
    EXPECT_LT(info, 0);
}

// Zero dimension test case (N1 = 0)
TEST_F(Ab05qdTestColMajor, ZeroDimensionN1) {
    int n1_zero = 0;
    int lda1 = std::max(1, n1_zero), ldb1 = std::max(1, n1_zero), ldc1 = P1, ldd1 = P1;
    int lda2 = N2, ldb2 = N2, ldc2 = P2, ldd2 = P2;
    int lda = std::max(1, n1_zero + N2);
    int ldb = std::max(1, n1_zero + N2);
    int ldc = P1 + P2;
    int ldd = P1 + P2;
    
    int n_expected = n1_zero + N2;
    int m_expected = M1 + M2;
    int p_expected = P1 + P2;
    
    // Create dummy arrays for the zero-dimension system
    std::vector<double> A1_zero(std::max(1, n1_zero * n1_zero), 0.0);
    std::vector<double> B1_zero(std::max(1, n1_zero * M1), 0.0);
    std::vector<double> C1_zero(std::max(1, P1 * n1_zero), 0.0);
    std::vector<double> D1_zero = D1;  // D1 is still P1×M1 when N1=0
    
    // Output arrays
    std::vector<double> A(std::max(1, (n1_zero + N2) * (n1_zero + N2)), 0.0);
    std::vector<double> B(std::max(1, (n1_zero + N2) * (M1 + M2)), 0.0);
    std::vector<double> C(std::max(1, (P1 + P2) * (n1_zero + N2)), 0.0);
    std::vector<double> D(std::max(1, (P1 + P2) * (M1 + M2)), 0.0);
    
    int n_out, m_out, p_out;
    
    // Call function with N1 = 0
    int info = slicot_ab05qd(OVER, n1_zero, M1, P1, N2, M2, P2,
                             A1_zero.data(), lda1, B1_zero.data(), ldb1,
                             C1_zero.data(), ldc1, D1_zero.data(), ldd1,
                             A2.data(), lda2, B2.data(), ldb2,
                             C2.data(), ldc2, D2.data(), ldd2,
                             &n_out, &m_out, &p_out,
                             A.data(), lda, B.data(), ldb,
                             C.data(), ldc, D.data(), ldd,
                             0); // Explicitly use column-major (0)
    
    // Verify
    ASSERT_EQ(info, 0);
    EXPECT_EQ(n_out, n_expected);
    EXPECT_EQ(m_out, m_expected);
    EXPECT_EQ(p_out, p_expected);
    
    // In this case, the result should just be system 2's matrices in the right places
    
    // Verify A matrix (should contain A2)
    for (int j = 0; j < N2; ++j) {
        for (int i = 0; i < N2; ++i) {
            EXPECT_NEAR(A[i + j*lda], A2[i + j*N2], check_tol)
                << "A[" << i << "," << j << "] mismatch";
        }
    }
    
    // Verify B matrix (should have zeros for first M1 columns, then B2)
    for (int j = 0; j < M1; ++j) {
        for (int i = 0; i < N2; ++i) {
            EXPECT_NEAR(B[i + j*ldb], 0.0, check_tol)
                << "B[" << i << "," << j << "] mismatch";
        }
    }
    for (int j = 0; j < M2; ++j) {
        for (int i = 0; i < N2; ++i) {
            EXPECT_NEAR(B[i + (j+M1)*ldb], B2[i + j*N2], check_tol)
                << "B[" << i << "," << (j+M1) << "] mismatch";
        }
    }
    
    // Verify D matrix (should have D1 in top-left and D2 in bottom-right)
    for (int j = 0; j < M1; ++j) {
        for (int i = 0; i < P1; ++i) {
            EXPECT_NEAR(D[i + j*ldd], D1[i + j*P1], check_tol)
                << "D[" << i << "," << j << "] mismatch";
        }
    }
    for (int j = 0; j < M2; ++j) {
        for (int i = 0; i < P2; ++i) {
            EXPECT_NEAR(D[(i+P1) + (j+M1)*ldd], D2[i + j*P2], check_tol)
                << "D[" << (i+P1) << "," << (j+M1) << "] mismatch";
        }
    }
}

// Zero dimension test case (N2 = 0)
TEST_F(Ab05qdTestColMajor, ZeroDimensionN2) {
    int n2_zero = 0;
    int lda1 = N1, ldb1 = N1, ldc1 = P1, ldd1 = P1;
    int lda2 = std::max(1, n2_zero), ldb2 = std::max(1, n2_zero), ldc2 = P2, ldd2 = P2;
    int lda = std::max(1, N1 + n2_zero);
    int ldb = std::max(1, N1 + n2_zero);
    int ldc = P1 + P2;
    int ldd = P1 + P2;
    
    int n_expected = N1 + n2_zero;
    int m_expected = M1 + M2;
    int p_expected = P1 + P2;
    
    // Create dummy arrays for the zero-dimension system
    std::vector<double> A2_zero(std::max(1, n2_zero * n2_zero), 0.0);
    std::vector<double> B2_zero(std::max(1, n2_zero * M2), 0.0);
    std::vector<double> C2_zero(std::max(1, P2 * n2_zero), 0.0);
    std::vector<double> D2_zero = D2;  // D2 is still P2×M2 when N2=0
    
    // Output arrays
    std::vector<double> A(std::max(1, (N1 + n2_zero) * (N1 + n2_zero)), 0.0);
    std::vector<double> B(std::max(1, (N1 + n2_zero) * (M1 + M2)), 0.0);
    std::vector<double> C(std::max(1, (P1 + P2) * (N1 + n2_zero)), 0.0);
    std::vector<double> D(std::max(1, (P1 + P2) * (M1 + M2)), 0.0);
    
    int n_out, m_out, p_out;
    
    // Call function with N2 = 0
    int info = slicot_ab05qd(OVER, N1, M1, P1, n2_zero, M2, P2,
                             A1.data(), lda1, B1.data(), ldb1,
                             C1.data(), ldc1, D1.data(), ldd1,
                             A2_zero.data(), lda2, B2_zero.data(), ldb2,
                             C2_zero.data(), ldc2, D2_zero.data(), ldd2,
                             &n_out, &m_out, &p_out,
                             A.data(), lda, B.data(), ldb,
                             C.data(), ldc, D.data(), ldd,
                             0); // Explicitly use column-major (0)
    
    // Verify
    ASSERT_EQ(info, 0);
    EXPECT_EQ(n_out, n_expected);
    EXPECT_EQ(m_out, m_expected);
    EXPECT_EQ(p_out, p_expected);
    
    // In this case, the result should just be system 1's matrices in the right places
    
    // Verify A matrix (should contain A1)
    for (int j = 0; j < N1; ++j) {
        for (int i = 0; i < N1; ++i) {
            EXPECT_NEAR(A[i + j*lda], A1[i + j*N1], check_tol)
                << "A[" << i << "," << j << "] mismatch";
        }
    }
    
    // Verify B matrix (should have B1 for first M1 columns, then zeros)
    for (int j = 0; j < M1; ++j) {
        for (int i = 0; i < N1; ++i) {
            EXPECT_NEAR(B[i + j*ldb], B1[i + j*N1], check_tol)
                << "B[" << i << "," << j << "] mismatch";
        }
    }
    for (int j = 0; j < M2; ++j) {
        for (int i = 0; i < N1; ++i) {
            EXPECT_NEAR(B[i + (j+M1)*ldb], 0.0, check_tol)
                << "B[" << i << "," << (j+M1) << "] mismatch";
        }
    }
    
    // Verify D matrix (should have D1 in top-left and D2 in bottom-right)
    for (int j = 0; j < M1; ++j) {
        for (int i = 0; i < P1; ++i) {
            EXPECT_NEAR(D[i + j*ldd], D1[i + j*P1], check_tol)
                << "D[" << i << "," << j << "] mismatch";
        }
    }
    for (int j = 0; j < M2; ++j) {
        for (int i = 0; i < P2; ++i) {
            EXPECT_NEAR(D[(i+P1) + (j+M1)*ldd], D2[i + j*P2], check_tol)
                << "D[" << (i+P1) << "," << (j+M1) << "] mismatch";
        }
    }
}
