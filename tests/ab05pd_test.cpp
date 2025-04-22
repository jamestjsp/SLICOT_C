#include <gtest/gtest.h>
#include <vector>
#include <cmath>

#include "ab05pd.h"

// Column-major test fixture for AB05PD
class Ab05pdTestColMajor : public ::testing::Test {
protected:
    // Test parameters from AB05PD.dat example
    int N1 = 3;
    int M = 2;
    int P = 2;
    int N2 = 3;
    double ALPHA = 1.0;
    int ROW_MAJOR = 0;  // Column-major
    int N = N1 + N2;    // Combined system order
    double check_tol = 1e-4;
    
    // Input matrices from documentation (column-major order)
    // A1 (3x3)
    std::vector<double> A1 = {
        /* Column 1 */ 1.0, 0.0, 1.0,
        /* Column 2 */ 0.0, -1.0, 1.0,
        /* Column 3 */ -1.0, 1.0, 2.0
    };
    
    // B1 (3x2)
    std::vector<double> B1 = {
        /* Column 1 */ 1.0, 1.0, 0.0,
        /* Column 2 */ 2.0, 0.0, 1.0
    };
    
    // C1 (2x3)
    std::vector<double> C1 = {
        /* Column 1 */ 3.0, 0.0,
        /* Column 2 */ -2.0, 1.0,
        /* Column 3 */ 1.0, 0.0
    };
    
    // D1 (2x2)
    std::vector<double> D1 = {
        /* Column 1 */ 1.0, 0.0,
        /* Column 2 */ 0.0, 1.0
    };
    
    // A2 (3x3)
    std::vector<double> A2 = {
        /* Column 1 */ -3.0, 1.0, 0.0,
        /* Column 2 */ 0.0, 0.0, -1.0,
        /* Column 3 */ 0.0, 1.0, 2.0
    };
    
    // B2 (3x2)
    std::vector<double> B2 = {
        /* Column 1 */ 0.0, -1.0, 0.0,
        /* Column 2 */ 1.0, 0.0, 2.0
    };
    
    // C2 (2x3)
    std::vector<double> C2 = {
        /* Column 1 */ 1.0, 1.0,
        /* Column 2 */ 1.0, 1.0,
        /* Column 3 */ 0.0, -1.0
    };
    
    // D2 (2x2)
    std::vector<double> D2 = {
        /* Column 1 */ 1.0, 0.0,
        /* Column 2 */ 1.0, 1.0
    };
    
    // Expected outputs (from AB05PD.res)
    // A (6x6) - Combined state transition matrix
    std::vector<double> A_expected = {
        /* Column 1 */ 1.0, 0.0, 1.0, 0.0, 0.0, 0.0,
        /* Column 2 */ 0.0, -1.0, 1.0, 0.0, 0.0, 0.0,
        /* Column 3 */ -1.0, 1.0, 2.0, 0.0, 0.0, 0.0,
        /* Column 4 */ 0.0, 0.0, 0.0, -3.0, 1.0, 0.0,
        /* Column 5 */ 0.0, 0.0, 0.0, 0.0, 0.0, -1.0,
        /* Column 6 */ 0.0, 0.0, 0.0, 0.0, 1.0, 2.0
    };
    
    // B (6x2) - Combined input/state matrix
    std::vector<double> B_expected = {
        /* Column 1 */ 1.0, 1.0, 0.0, 0.0, -1.0, 0.0,
        /* Column 2 */ 2.0, 0.0, 1.0, 1.0, 0.0, 2.0
    };
    
    // C (2x6) - Combined state/output matrix
    std::vector<double> C_expected = {
        /* Column 1 */ 3.0, 0.0,
        /* Column 2 */ -2.0, 1.0,
        /* Column 3 */ 1.0, 0.0,
        /* Column 4 */ 1.0, 1.0,
        /* Column 5 */ 1.0, 1.0,
        /* Column 6 */ 0.0, -1.0
    };
    
    // D (2x2) - Combined input/output matrix
    std::vector<double> D_expected = {
        /* Column 1 */ 2.0, 0.0,
        /* Column 2 */ 1.0, 2.0
    };
    
    int expected_info = 0;
};

// Row-major test fixture
class Ab05pdTestRowMajor : public Ab05pdTestColMajor {
public:
    Ab05pdTestRowMajor() {
        ROW_MAJOR = 1;  // Row-major
        
        // Convert test matrices from column-major to row-major
        A1_rm.resize(N1 * N1);
        B1_rm.resize(N1 * M);
        C1_rm.resize(P * N1);
        D1_rm.resize(P * M);
        A2_rm.resize(N2 * N2);
        B2_rm.resize(N2 * M);
        C2_rm.resize(P * N2);
        D2_rm.resize(P * M);
        A_expected_rm.resize(N * N);
        B_expected_rm.resize(N * M);
        C_expected_rm.resize(P * N);
        D_expected_rm.resize(P * M);
        
        // A1: N1xN1
        for (int i = 0; i < N1; ++i) {
            for (int j = 0; j < N1; ++j) {
                A1_rm[i*N1 + j] = A1[i + j*N1];
            }
        }
        
        // B1: N1xM
        for (int i = 0; i < N1; ++i) {
            for (int j = 0; j < M; ++j) {
                B1_rm[i*M + j] = B1[i + j*N1];
            }
        }
        
        // C1: PxN1
        for (int i = 0; i < P; ++i) {
            for (int j = 0; j < N1; ++j) {
                C1_rm[i*N1 + j] = C1[i + j*P];
            }
        }
        
        // D1: PxM
        for (int i = 0; i < P; ++i) {
            for (int j = 0; j < M; ++j) {
                D1_rm[i*M + j] = D1[i + j*P];
            }
        }
        
        // A2: N2xN2
        for (int i = 0; i < N2; ++i) {
            for (int j = 0; j < N2; ++j) {
                A2_rm[i*N2 + j] = A2[i + j*N2];
            }
        }
        
        // B2: N2xM
        for (int i = 0; i < N2; ++i) {
            for (int j = 0; j < M; ++j) {
                B2_rm[i*M + j] = B2[i + j*N2];
            }
        }
        
        // C2: PxN2
        for (int i = 0; i < P; ++i) {
            for (int j = 0; j < N2; ++j) {
                C2_rm[i*N2 + j] = C2[i + j*P];
            }
        }
        
        // D2: PxM
        for (int i = 0; i < P; ++i) {
            for (int j = 0; j < M; ++j) {
                D2_rm[i*M + j] = D2[i + j*P];
            }
        }
        
        // Expected outputs
        // A: NxN
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                A_expected_rm[i*N + j] = A_expected[i + j*N];
            }
        }
        
        // B: NxM
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                B_expected_rm[i*M + j] = B_expected[i + j*N];
            }
        }
        
        // C: PxN
        for (int i = 0; i < P; ++i) {
            for (int j = 0; j < N; ++j) {
                C_expected_rm[i*N + j] = C_expected[i + j*P];
            }
        }
        
        // D: PxM
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

// Column-major test case from documentation example
TEST_F(Ab05pdTestColMajor, DocExample) {
    // Set up output buffers
    int n_out = 0;
    std::vector<double> A(N * N);
    std::vector<double> B(N * M);
    std::vector<double> C(P * N);
    std::vector<double> D(P * M);
    
    // Leading dimensions for column-major format
    int LDA1 = N1, LDB1 = N1, LDC1 = P, LDD1 = P;
    int LDA2 = N2, LDB2 = N2, LDC2 = P, LDD2 = P;
    int LDA = N, LDB = N, LDC = P, LDD = P;
    
    // Call function under test
    char OVER = 'N'; // No overlap
    int info = slicot_ab05pd(OVER,
                           N1, M, P, N2, ALPHA,
                           A1.data(), LDA1, B1.data(), LDB1,
                           C1.data(), LDC1, D1.data(), LDD1,
                           A2.data(), LDA2, B2.data(), LDB2,
                           C2.data(), LDC2, D2.data(), LDD2,
                           &n_out,
                           A.data(), LDA, B.data(), LDB,
                           C.data(), LDC, D.data(), LDD,
                           ROW_MAJOR);
    
    // Verify function result
    ASSERT_EQ(info, expected_info);
    ASSERT_EQ(n_out, N);
    
    // Verify output matrices (column major order)
    // Verify A
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            EXPECT_NEAR(A[i + j*LDA], A_expected[i + j*N], check_tol)
                << "A[" << i << "," << j << "] mismatch";
        }
    }
    
    // Verify B
    for (int j = 0; j < M; ++j) {
        for (int i = 0; i < N; ++i) {
            EXPECT_NEAR(B[i + j*LDB], B_expected[i + j*N], check_tol)
                << "B[" << i << "," << j << "] mismatch";
        }
    }
    
    // Verify C
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < P; ++i) {
            EXPECT_NEAR(C[i + j*LDC], C_expected[i + j*P], check_tol)
                << "C[" << i << "," << j << "] mismatch";
        }
    }
    
    // Verify D
    for (int j = 0; j < M; ++j) {
        for (int i = 0; i < P; ++i) {
            EXPECT_NEAR(D[i + j*LDD], D_expected[i + j*P], check_tol)
                << "D[" << i << "," << j << "] mismatch";
        }
    }
}

// Row-major test case from documentation example
TEST_F(Ab05pdTestRowMajor, DocExample) {
    // Set up output buffers
    int n_out = 0;
    std::vector<double> A(N * N);
    std::vector<double> B(N * M);
    std::vector<double> C(P * N);
    std::vector<double> D(P * M);
    
    // Leading dimensions for row-major format
    int LDA1 = N1, LDB1 = M, LDC1 = N1, LDD1 = M;
    int LDA2 = N2, LDB2 = M, LDC2 = N2, LDD2 = M;
    int LDA = N, LDB = M, LDC = N, LDD = M;
    
    // Call function under test
    char OVER = 'N'; // No overlap
    int info = slicot_ab05pd(OVER,
                           N1, M, P, N2, ALPHA,
                           A1_rm.data(), LDA1, B1_rm.data(), LDB1,
                           C1_rm.data(), LDC1, D1_rm.data(), LDD1,
                           A2_rm.data(), LDA2, B2_rm.data(), LDB2,
                           C2_rm.data(), LDC2, D2_rm.data(), LDD2,
                           &n_out,
                           A.data(), LDA, B.data(), LDB,
                           C.data(), LDC, D.data(), LDD,
                           ROW_MAJOR);
    
    // Verify function result
    ASSERT_EQ(info, expected_info);
    ASSERT_EQ(n_out, N);
    
    // Verify output matrices (row major order)
    // Verify A
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            EXPECT_NEAR(A[i*LDA + j], A_expected_rm[i*N + j], check_tol)
                << "A[" << i << "," << j << "] mismatch";
        }
    }
    
    // Verify B
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            EXPECT_NEAR(B[i*LDB + j], B_expected_rm[i*M + j], check_tol)
                << "B[" << i << "," << j << "] mismatch";
        }
    }
    
    // Verify C
    for (int i = 0; i < P; ++i) {
        for (int j = 0; j < N; ++j) {
            EXPECT_NEAR(C[i*LDC + j], C_expected_rm[i*N + j], check_tol)
                << "C[" << i << "," << j << "] mismatch";
        }
    }
    
    // Verify D
    for (int i = 0; i < P; ++i) {
        for (int j = 0; j < M; ++j) {
            EXPECT_NEAR(D[i*LDD + j], D_expected_rm[i*M + j], check_tol)
                << "D[" << i << "," << j << "] mismatch";
        }
    }
}

// Test parameter validation
TEST_F(Ab05pdTestColMajor, ParameterValidation) {
    // Set up minimal buffers
    int n_out = 0;
    std::vector<double> A(N * N);
    std::vector<double> B(N * M);
    std::vector<double> C(P * N);
    std::vector<double> D(P * M);
    
    // Leading dimensions
    int LDA1 = N1, LDB1 = N1, LDC1 = P, LDD1 = P;
    int LDA2 = N2, LDB2 = N2, LDC2 = P, LDD2 = P;
    int LDA = N, LDB = N, LDC = P, LDD = P;
    
    // Test invalid parameters
    char OVER = 'N';
    
    // Test negative N1
    int info = slicot_ab05pd(OVER,
                         -1, M, P, N2, ALPHA,
                         A1.data(), LDA1, B1.data(), LDB1,
                         C1.data(), LDC1, D1.data(), LDD1,
                         A2.data(), LDA2, B2.data(), LDB2,
                         C2.data(), LDC2, D2.data(), LDD2,
                         &n_out,
                         A.data(), LDA, B.data(), LDB,
                         C.data(), LDC, D.data(), LDD,
                         ROW_MAJOR);
    ASSERT_EQ(info, -2); // N1 is 2nd parameter
    
    // Test negative M
    info = slicot_ab05pd(OVER,
                      N1, -1, P, N2, ALPHA,
                      A1.data(), LDA1, B1.data(), LDB1,
                      C1.data(), LDC1, D1.data(), LDD1,
                      A2.data(), LDA2, B2.data(), LDB2,
                      C2.data(), LDC2, D2.data(), LDD2,
                      &n_out,
                      A.data(), LDA, B.data(), LDB,
                      C.data(), LDC, D.data(), LDD,
                      ROW_MAJOR);
    ASSERT_EQ(info, -3); // M is 3rd parameter
    
    // Test negative P
    info = slicot_ab05pd(OVER,
                      N1, M, -1, N2, ALPHA,
                      A1.data(), LDA1, B1.data(), LDB1,
                      C1.data(), LDC1, D1.data(), LDD1,
                      A2.data(), LDA2, B2.data(), LDB2,
                      C2.data(), LDC2, D2.data(), LDD2,
                      &n_out,
                      A.data(), LDA, B.data(), LDB,
                      C.data(), LDC, D.data(), LDD,
                      ROW_MAJOR);
    ASSERT_EQ(info, -4); // P is 4th parameter
    
    // Test negative N2
    info = slicot_ab05pd(OVER,
                      N1, M, P, -1, ALPHA,
                      A1.data(), LDA1, B1.data(), LDB1,
                      C1.data(), LDC1, D1.data(), LDD1,
                      A2.data(), LDA2, B2.data(), LDB2,
                      C2.data(), LDC2, D2.data(), LDD2,
                      &n_out,
                      A.data(), LDA, B.data(), LDB,
                      C.data(), LDC, D.data(), LDD,
                      ROW_MAJOR);
    ASSERT_EQ(info, -5); // N2 is 5th parameter
    
    // Test invalid OVER parameter
    info = slicot_ab05pd('X', // Invalid OVER parameter
                      N1, M, P, N2, ALPHA,
                      A1.data(), LDA1, B1.data(), LDB1,
                      C1.data(), LDC1, D1.data(), LDD1,
                      A2.data(), LDA2, B2.data(), LDB2,
                      C2.data(), LDC2, D2.data(), LDD2,
                      &n_out,
                      A.data(), LDA, B.data(), LDB,
                      C.data(), LDC, D.data(), LDD,
                      ROW_MAJOR);
    ASSERT_EQ(info, -1); // OVER is 1st parameter
}

// Test edge case with zero dimensions
TEST_F(Ab05pdTestColMajor, ZeroDimension) {
    // Set up with zero dimensions
    int n1_zero = 0;
    int m_zero = 0;
    int p_zero = 0;
    int n2_zero = 0;
    int n_out = 0;
    
    // Create minimal valid arrays (still need space for 1 element to avoid nullptr)
    std::vector<double> A1_zero(1, 0.0);
    std::vector<double> B1_zero(1, 0.0);
    std::vector<double> C1_zero(1, 0.0);
    std::vector<double> D1_zero(1, 0.0);
    std::vector<double> A2_zero(1, 0.0);
    std::vector<double> B2_zero(1, 0.0);
    std::vector<double> C2_zero(1, 0.0);
    std::vector<double> D2_zero(1, 0.0);
    std::vector<double> A_zero(MAX(1, N1+N2), 0.0); // Ensure enough space for combined system
    std::vector<double> B_zero(MAX(1, (N1+N2)*M), 0.0);
    std::vector<double> C_zero(MAX(1, P*(N1+N2)), 0.0);
    std::vector<double> D_zero(MAX(1, P*M), 0.0);
    
    // For the N1=0 test
    // Test with N1=0, col-major: A2 is N2×N2, B2 is N2×M, C2 is P×N2, D2 is P×M
    int info = slicot_ab05pd('N',
                         n1_zero, M, P, N2, ALPHA,
                         A1_zero.data(), 1, B1_zero.data(), 1,
                         C1_zero.data(), MAX(1,P), D1_zero.data(), MAX(1,P),
                         A2.data(), N2, B2.data(), N2, // Use actual dimensions from A2 and B2
                         C2.data(), P, D2.data(), P,  // Use actual dimensions from C2 and D2
                         &n_out,                         A_zero.data(), MAX(1,N2), B_zero.data(), MAX(1,N2),
                         C_zero.data(), MAX(1,P), D_zero.data(), MAX(1,P),
                         0); // Column-major format (0)
    ASSERT_EQ(info, 0);
    ASSERT_EQ(n_out, N2);    ASSERT_EQ(info, 0);
    ASSERT_EQ(n_out, N2);
    
    // Test with N2=0
    info = slicot_ab05pd('N',
                      N1, M, P, n2_zero, ALPHA,
                      A1.data(), N1, B1.data(), N1,
                      C1.data(), P, D1.data(), P,
                      A2_zero.data(), 1, B2_zero.data(), 1,
                      C2_zero.data(), MAX(1,P), D2_zero.data(), MAX(1,P),
                      &n_out,
                      A_zero.data(), MAX(1,N1), B_zero.data(), MAX(1,N1),
                      C_zero.data(), MAX(1,P), D_zero.data(), MAX(1,P),
                      0); // Column-major format
    ASSERT_EQ(info, 0);
    ASSERT_EQ(n_out, N1);
    
    // Test with M=0
    info = slicot_ab05pd('N',
                      N1, m_zero, P, N2, ALPHA,
                      A1.data(), N1, B1_zero.data(), N1,
                      C1.data(), P, D1_zero.data(), P,
                      A2.data(), N2, B2_zero.data(), N2,
                      C2.data(), P, D2_zero.data(), P,
                      &n_out,
                      A_zero.data(), MAX(1,N1+N2), B_zero.data(), MAX(1,N1+N2),
                      C_zero.data(), MAX(1,P), D_zero.data(), MAX(1,P),
                      0); // Column-major format
    ASSERT_EQ(info, 0);
    
    // Test with P=0
    info = slicot_ab05pd('N',
                      N1, M, p_zero, N2, ALPHA,
                      A1.data(), N1, B1.data(), N1,
                      C1_zero.data(), 1, D1_zero.data(), 1,
                      A2.data(), N2, B2.data(), N2,
                      C2_zero.data(), 1, D2_zero.data(), 1,
                      &n_out,
                      A_zero.data(), MAX(1,N1+N2), B_zero.data(), MAX(1,N1+N2),
                      C_zero.data(), 1, D_zero.data(), 1,
                      0); // Column-major format
    ASSERT_EQ(info, 0);
}
