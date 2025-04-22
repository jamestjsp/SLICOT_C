#include <gtest/gtest.h>
#include <vector>
#include <cmath>

#include "ab07nd.h"

// Column-major test fixture for AB07ND
class AB07NDTestColMajor : public ::testing::Test {
protected:
    // Test parameters from example
    int N = 3;
    int M = 2;
    int ROW_MAJOR = 0;  // Column-major
    double check_tol = 1e-4;
    
    // Input matrices from documentation (column-major order)
    std::vector<double> A_in = {
        /* Column 1 */ 1.0, 4.0, 0.0,
        /* Column 2 */ 2.0, -1.0, 0.0,
        /* Column 3 */ 0.0, 0.0, 1.0
    };
    
    std::vector<double> B_in = {
        /* Column 1 */ 1.0, 0.0, 1.0,
        /* Column 2 */ 0.0, 1.0, 0.0
    };
    
    std::vector<double> C_in = {
        /* Column 1 */ 0.0, 0.0,
        /* Column 2 */ 1.0, 0.0,
        /* Column 3 */ -1.0, 1.0
    };
    
    std::vector<double> D_in = {
        /* Column 1 */ 4.0, 0.0,
        /* Column 2 */ 0.0, 1.0
    };
    
    // Expected outputs (column-major order)
    std::vector<double> A_expected = {
        /* Column 1 */ 1.0000, 4.0000, 0.0000,
        /* Column 2 */ 1.7500, -1.0000, -0.2500,
        /* Column 3 */ 0.2500, -1.0000, 1.2500
    };
    
    std::vector<double> B_expected = {
        /* Column 1 */ -0.2500, 0.0000, -0.2500,
        /* Column 2 */ 0.0000, -1.0000, 0.0000
    };
    
    std::vector<double> C_expected = {
        /* Column 1 */ 0.0000, 0.0000,
        /* Column 2 */ 0.2500, 0.0000,
        /* Column 3 */ -0.2500, 1.0000
    };
    
    std::vector<double> D_expected = {
        /* Column 1 */ 0.2500, 0.0000,
        /* Column 2 */ 0.0000, 1.0000
    };
    
    int expected_info = 0;
};

// Row-major test fixture
class AB07NDTestRowMajor : public AB07NDTestColMajor {
public:
    AB07NDTestRowMajor() {
        ROW_MAJOR = 1;  // Row-major
        
        // Convert test data to row-major
        A_rm.resize(N * N);
        B_rm.resize(N * M);
        C_rm.resize(M * N);
        D_rm.resize(M * M);
        
        A_expected_rm.resize(N * N);
        B_expected_rm.resize(N * M);
        C_expected_rm.resize(M * N);
        D_expected_rm.resize(M * M);
        
        // Convert A from column-major to row-major
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                A_rm[i*N + j] = A_in[i + j*N];
                A_expected_rm[i*N + j] = A_expected[i + j*N];
            }
        }
        
        // Convert B from column-major to row-major
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                B_rm[i*M + j] = B_in[i + j*N];
                B_expected_rm[i*M + j] = B_expected[i + j*N];
            }
        }
        
        // Convert C from column-major to row-major
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < N; ++j) {
                C_rm[i*N + j] = C_in[i + j*M];
                C_expected_rm[i*N + j] = C_expected[i + j*M];
            }
        }
        
        // Convert D from column-major to row-major
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < M; ++j) {
                D_rm[i*M + j] = D_in[i + j*M];
                D_expected_rm[i*M + j] = D_expected[i + j*M];
            }
        }
    }
    
    std::vector<double> A_rm, B_rm, C_rm, D_rm;
    std::vector<double> A_expected_rm, B_expected_rm, C_expected_rm, D_expected_rm;
};

// Column-major test
TEST_F(AB07NDTestColMajor, DocExample) {
    // Set up
    int LDA = N;
    int LDB = N;
    int LDC = M;
    int LDD = M;
    
    std::vector<double> A = A_in;
    std::vector<double> B = B_in;
    std::vector<double> C = C_in;
    std::vector<double> D = D_in;
    double rcond;
    
    // Call function
    int info = slicot_ab07nd(N, M, A.data(), LDA, B.data(), LDB, 
                            C.data(), LDC, D.data(), LDD, &rcond, ROW_MAJOR);
    
    // Verify
    ASSERT_EQ(info, expected_info);
    
    // Check that rcond is reasonable (non-zero for invertible D)
    EXPECT_GT(rcond, 0.0);
    
    // Check output matrices
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            EXPECT_NEAR(A[i + j*LDA], A_expected[i + j*N], check_tol)
                << "A[" << i << "," << j << "] mismatch";
        }
    }
    
    for (int j = 0; j < M; ++j) {
        for (int i = 0; i < N; ++i) {
            EXPECT_NEAR(B[i + j*LDB], B_expected[i + j*N], check_tol)
                << "B[" << i << "," << j << "] mismatch";
        }
    }
    
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < M; ++i) {
            EXPECT_NEAR(C[i + j*LDC], C_expected[i + j*M], check_tol)
                << "C[" << i << "," << j << "] mismatch";
        }
    }
    
    for (int j = 0; j < M; ++j) {
        for (int i = 0; i < M; ++i) {
            EXPECT_NEAR(D[i + j*LDD], D_expected[i + j*M], check_tol)
                << "D[" << i << "," << j << "] mismatch";
        }
    }
}

// Row-major test
TEST_F(AB07NDTestRowMajor, DocExample) {
    // Set up
    int LDA = N;
    int LDB = M;
    int LDC = N;
    int LDD = M;
    
    std::vector<double> A = A_rm;
    std::vector<double> B = B_rm;
    std::vector<double> C = C_rm;
    std::vector<double> D = D_rm;
    double rcond;
    
    // Call function
    int info = slicot_ab07nd(N, M, A.data(), LDA, B.data(), LDB, 
                            C.data(), LDC, D.data(), LDD, &rcond, ROW_MAJOR);
    
    // Verify
    ASSERT_EQ(info, expected_info);
    
    // Check that rcond is reasonable (non-zero for invertible D)
    EXPECT_GT(rcond, 0.0);
    
    // Check output matrices
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            EXPECT_NEAR(A[i*LDA + j], A_expected_rm[i*N + j], check_tol)
                << "A[" << i << "," << j << "] mismatch";
        }
    }
    
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            EXPECT_NEAR(B[i*LDB + j], B_expected_rm[i*M + j], check_tol)
                << "B[" << i << "," << j << "] mismatch";
        }
    }
    
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            EXPECT_NEAR(C[i*LDC + j], C_expected_rm[i*N + j], check_tol)
                << "C[" << i << "," << j << "] mismatch";
        }
    }
    
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M; ++j) {
            EXPECT_NEAR(D[i*LDD + j], D_expected_rm[i*M + j], check_tol)
                << "D[" << i << "," << j << "] mismatch";
        }
    }
}

// Test for edge case with N=0
TEST_F(AB07NDTestColMajor, ZeroN) {
    // Set up
    int n_zero = 0;
    
    std::vector<double> A_zero(1, 0.0); // At least one element
    std::vector<double> B_zero(1, 0.0);
    std::vector<double> C_zero(M, 0.0);
    std::vector<double> D = D_in; // D remains the same size
    
    int LDA = 1;
    int LDB = 1;
    int LDC = M;
    int LDD = M;
    double rcond;
    
    // Call function
    int info = slicot_ab07nd(n_zero, M, A_zero.data(), LDA, B_zero.data(), LDB, 
                            C_zero.data(), LDC, D.data(), LDD, &rcond, ROW_MAJOR);
    
    // Verify
    ASSERT_EQ(info, 0);
    
    // Only D should be inverted when N=0
    for (int j = 0; j < M; ++j) {
        for (int i = 0; i < M; ++i) {
            if (i == j) {
                // D inverse diagonal should be reciprocal of original diagonal values
                EXPECT_NEAR(D[i + j*LDD], 1.0/D_in[i + j*M], check_tol)
                    << "D[" << i << "," << j << "] mismatch";
            } else {
                // D inverse off-diagonal should be -D_orig[i,j]/D_orig[i,i]/D_orig[j,j]
                double expected = -D_in[i + j*M] / (D_in[i + i*M] * D_in[j + j*M]);
                EXPECT_NEAR(D[i + j*LDD], expected, check_tol)
                    << "D[" << i << "," << j << "] mismatch";
            }
        }
    }
}

// Test for edge case with M=0
TEST_F(AB07NDTestColMajor, ZeroM) {
    // Set up
    int m_zero = 0;
    
    std::vector<double> A = A_in; // A remains the same size
    std::vector<double> B_zero(N, 0.0);
    std::vector<double> C_zero(1, 0.0);
    std::vector<double> D_zero(1, 0.0);
    
    int LDA = N;
    int LDB = N;
    int LDC = 1;
    int LDD = 1;
    double rcond;
    
    // Call function
    int info = slicot_ab07nd(N, m_zero, A.data(), LDA, B_zero.data(), LDB, 
                            C_zero.data(), LDC, D_zero.data(), LDD, &rcond, ROW_MAJOR);
    
    // Verify
    ASSERT_EQ(info, 0);
    
    // With M=0, A should remain unchanged
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            EXPECT_NEAR(A[i + j*LDA], A_in[i + j*N], check_tol)
                << "A[" << i << "," << j << "] mismatch";
        }
    }
}

// Test for error case (singular D matrix)
TEST_F(AB07NDTestColMajor, SingularD) {
    // Set up
    std::vector<double> A = A_in;
    std::vector<double> B = B_in;
    std::vector<double> C = C_in;
    // Create a singular D matrix (zero determinant)
    std::vector<double> D_singular = {
        /* Column 1 */ 1.0, 0.0,
        /* Column 2 */ 0.0, 0.0  // Second column all zeros makes D singular
    };
    
    int LDA = N;
    int LDB = N;
    int LDC = M;
    int LDD = M;
    double rcond;
    
    // Call function
    int info = slicot_ab07nd(N, M, A.data(), LDA, B.data(), LDB, 
                            C.data(), LDC, D_singular.data(), LDD, &rcond, ROW_MAJOR);
    
    // Expect error info = 2 (the second diagonal element is zero)
    EXPECT_EQ(info, 2);
    // rcond should be 0 for singular matrix
    EXPECT_DOUBLE_EQ(rcond, 0.0);
}

// Input validation test
TEST_F(AB07NDTestColMajor, InputValidation) {
    // Set up
    std::vector<double> A = A_in;
    std::vector<double> B = B_in;
    std::vector<double> C = C_in;
    std::vector<double> D = D_in;
    double rcond;
    
    // Test negative N
    int info = slicot_ab07nd(-1, M, A.data(), N, B.data(), N, 
                           C.data(), M, D.data(), M, &rcond, ROW_MAJOR);
    EXPECT_EQ(info, -1);
    
    // Test negative M
    info = slicot_ab07nd(N, -1, A.data(), N, B.data(), N, 
                        C.data(), M, D.data(), M, &rcond, ROW_MAJOR);
    EXPECT_EQ(info, -2);
    
    // Test LDA < N (column-major)
    info = slicot_ab07nd(N, M, A.data(), N-1, B.data(), N, 
                        C.data(), M, D.data(), M, &rcond, 0);
    EXPECT_EQ(info, -4);
}
