#include <gtest/gtest.h>
#include <vector>
#include <cmath>

#include "ab08md.h"

// Column-major test fixture for AB08MD
class AB08MDTestColMajor : public ::testing::Test {
protected:
    // Test parameters
    char EQUIL = 'N';  // No balancing in base test
    int N = 3;         // Order of state matrix A
    int M = 2;         // Number of inputs
    int P = 2;         // Number of outputs
    int ROW_MAJOR = 0; // Column-major
    double TOL = 0.0;  // Default tolerance
    double check_tol = 1e-14; // Tolerance for result verification    // Input matrices (column-major order)
    // System with a specific rank of 2 for the transfer function matrix G(s) = C(sI-A)^(-1)B + D
    // Design rationale:
    // 1. A has a specific structure with one zero eigenvalue to test rank calculation with singularities
    // 2. B and C are designed to create two independent modes in the transfer function
    // 3. Careful combination of values in A, B, C, D to ensure the system has rank 2
    std::vector<double> A_in = {
        /* Column 1 */ 1.0, 0.0, 0.0,  // First eigenvalue = 1
        /* Column 2 */ 0.0, 2.0, 0.0,  // Second eigenvalue = 2
        /* Column 3 */ 0.0, 0.0, 0.0   // Third eigenvalue = 0 (singular)
    };
    
    std::vector<double> B_in = {
        /* Column 1 */ 1.0, 0.0, 0.0,  // First input affects only first state
        /* Column 2 */ 0.0, 1.0, 1.0   // Second input affects second and third states
    };
    
    std::vector<double> C_in = {
        /* Column 1 */ 1.0, 0.0,       // First state visible to first output
        /* Column 2 */ 0.0, 1.0,       // Second state visible to second output
        /* Column 3 */ 0.0, 0.0        // Third state not visible (important for rank calculation)
    };
    
    std::vector<double> D_in = {
        /* Column 1 */ 0.0, 1.0,       // Direct feedthrough from first input to second output
        /* Column 2 */ 0.0, 0.0        // No direct feedthrough from second input
    };
    
    // Expected output
    int expected_rank = 2;  // Expected rank of transfer function matrix
    int expected_info = 0;  // Expected return code
};

// Row-major test fixture
class AB08MDTestRowMajor : public AB08MDTestColMajor {
public:
    AB08MDTestRowMajor() {
        ROW_MAJOR = 1;  // Row-major
        
        // Convert test data to row-major
        A_rm.resize(N * N);
        B_rm.resize(N * M);
        C_rm.resize(P * N);
        D_rm.resize(P * M);
        
        // Convert A from column-major to row-major
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                A_rm[i*N + j] = A_in[i + j*N];
            }
        }
        
        // Convert B from column-major to row-major
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                B_rm[i*M + j] = B_in[i + j*N];
            }
        }
        
        // Convert C from column-major to row-major
        for (int i = 0; i < P; ++i) {
            for (int j = 0; j < N; ++j) {
                C_rm[i*N + j] = C_in[i + j*P];
            }
        }
        
        // Convert D from column-major to row-major
        for (int i = 0; i < P; ++i) {
            for (int j = 0; j < M; ++j) {
                D_rm[i*M + j] = D_in[i + j*P];
            }
        }
    }
    
    std::vector<double> A_rm, B_rm, C_rm, D_rm;
};

// Basic test with column-major data
TEST_F(AB08MDTestColMajor, BasicTest) {
    // Set up
    int LDA = N;
    int LDB = N;
    int LDC = P;
    int LDD = P;
    
    std::vector<double> A = A_in;
    std::vector<double> B = B_in;
    std::vector<double> C = C_in;
    std::vector<double> D = D_in;
    int rank = -1;
    
    // Call function
    int info = slicot_ab08md(EQUIL, N, M, P, A.data(), LDA, B.data(), LDB,
                             C.data(), LDC, D.data(), LDD, &rank, TOL, ROW_MAJOR);
    
    // Verify
    ASSERT_EQ(info, expected_info);
    ASSERT_EQ(rank, expected_rank);
}

// Test with row-major data
TEST_F(AB08MDTestRowMajor, BasicTest) {
    // Set up
    int LDA = N;
    int LDB = M;
    int LDC = N;
    int LDD = M;
    
    std::vector<double> A = A_rm;
    std::vector<double> B = B_rm;
    std::vector<double> C = C_rm;
    std::vector<double> D = D_rm;
    int rank = -1;
    
    // Call function
    int info = slicot_ab08md(EQUIL, N, M, P, A.data(), LDA, B.data(), LDB,
                             C.data(), LDC, D.data(), LDD, &rank, TOL, ROW_MAJOR);
    
    // Verify
    ASSERT_EQ(info, expected_info);
    ASSERT_EQ(rank, expected_rank);
}

// Test with equil='S' (perform balancing)
TEST_F(AB08MDTestColMajor, WithBalancing) {
    // Set up
    int LDA = N;
    int LDB = N;
    int LDC = P;
    int LDD = P;
    
    std::vector<double> A = A_in;
    std::vector<double> B = B_in;
    std::vector<double> C = C_in;
    std::vector<double> D = D_in;
    int rank = -1;
    char equil = 'S';  // With balancing
    
    // Call function
    int info = slicot_ab08md(equil, N, M, P, A.data(), LDA, B.data(), LDB,
                             C.data(), LDC, D.data(), LDD, &rank, TOL, ROW_MAJOR);
    
    // Verify
    ASSERT_EQ(info, expected_info);
    ASSERT_EQ(rank, expected_rank);
    // Note: A, B, C, D might be modified due to scaling
}

// Test with equil='S' in row-major
TEST_F(AB08MDTestRowMajor, WithBalancing) {
    // Set up
    int LDA = N;
    int LDB = M;
    int LDC = N;
    int LDD = M;
    
    std::vector<double> A = A_rm;
    std::vector<double> B = B_rm;
    std::vector<double> C = C_rm;
    std::vector<double> D = D_rm;
    int rank = -1;
    char equil = 'S';  // With balancing
    
    // Call function
    int info = slicot_ab08md(equil, N, M, P, A.data(), LDA, B.data(), LDB,
                             C.data(), LDC, D.data(), LDD, &rank, TOL, ROW_MAJOR);
    
    // Verify
    ASSERT_EQ(info, expected_info);
    ASSERT_EQ(rank, expected_rank);
    // Note: A, B, C, D might be modified due to scaling
}

// Test with specific tolerance
TEST_F(AB08MDTestColMajor, SpecificTolerance) {
    // Set up
    int LDA = N;
    int LDB = N;
    int LDC = P;
    int LDD = P;
    
    std::vector<double> A = A_in;
    std::vector<double> B = B_in;
    std::vector<double> C = C_in;
    std::vector<double> D = D_in;
    int rank = -1;
    double tol = 1e-10;  // Specific tolerance
    
    // Call function
    int info = slicot_ab08md(EQUIL, N, M, P, A.data(), LDA, B.data(), LDB,
                             C.data(), LDC, D.data(), LDD, &rank, tol, ROW_MAJOR);
    
    // Verify
    ASSERT_EQ(info, expected_info);
    ASSERT_EQ(rank, expected_rank);
}

// Test with full-rank system
TEST_F(AB08MDTestColMajor, FullRankSystem) {
    // Set up
    int LDA = N;
    int LDB = N;
    int LDC = P;
    int LDD = P;
    
    // Create full-rank A matrix
    std::vector<double> A_full = {
        /* Column 1 */ 1.0, 0.0, 0.0,
        /* Column 2 */ 0.0, 2.0, 0.0,
        /* Column 3 */ 0.0, 0.0, 3.0  // All eigenvalues non-zero
    };
    
    std::vector<double> B = B_in;
    std::vector<double> C = C_in;
    
    // D matrix with full rank
    std::vector<double> D_full = {
        /* Column 1 */ 1.0, 0.0,
        /* Column 2 */ 0.0, 1.0
    };
    
    int rank = -1;
    int expected_full_rank = 2; // Full rank for this system
    
    // Call function
    int info = slicot_ab08md(EQUIL, N, M, P, A_full.data(), LDA, B.data(), LDB,
                             C.data(), LDC, D_full.data(), LDD, &rank, TOL, ROW_MAJOR);
    
    // Verify
    ASSERT_EQ(info, expected_info);
    ASSERT_EQ(rank, expected_full_rank);
}

// Zero dimension tests
TEST_F(AB08MDTestColMajor, ZeroDimensionN) {
    // Test with N=0 (zero state variables)
    // For a system with N=0, the transfer function G(s) = C(sI-A)^(-1)B + D reduces to just D
    // So the rank is determined entirely by the D matrix
    int n_zero = 0;
    int LDA = 1;
    int LDB = 1;
    int LDC = P;
    int LDD = P;
    
    std::vector<double> A_empty(1, 0.0); // At least one element (not used)
    std::vector<double> B_empty(1, 0.0); // At least one element (not used)
    std::vector<double> C_empty(1, 0.0); // At least one element (not used)
    std::vector<double> D = D_in;        // Using the same D matrix which has rank 1
    
    int rank = -1;
    int expected_zero_rank = 1;  // D matrix has rank 1, so system has rank 1
    
    int info = slicot_ab08md(EQUIL, n_zero, M, P, A_empty.data(), LDA, B_empty.data(), LDB,
                             C_empty.data(), LDC, D.data(), LDD, &rank, TOL, ROW_MAJOR);
    
    ASSERT_EQ(info, expected_info);
    ASSERT_EQ(rank, expected_zero_rank);
}

TEST_F(AB08MDTestColMajor, ZeroDimensionM) {
    // Test with M=0
    int m_zero = 0;
    int LDA = N;
    int LDB = N;
    int LDC = P;
    int LDD = P;
    
    std::vector<double> A = A_in;
    std::vector<double> B_empty(1, 0.0);
    std::vector<double> C = C_in;
    std::vector<double> D_empty(1, 0.0);
    
    int rank = -1;
    int expected_zero_rank = 0;
    
    int info = slicot_ab08md(EQUIL, N, m_zero, P, A.data(), LDA, B_empty.data(), LDB,
                             C.data(), LDC, D_empty.data(), LDD, &rank, TOL, ROW_MAJOR);
    
    ASSERT_EQ(info, expected_info);
    ASSERT_EQ(rank, expected_zero_rank);
}

TEST_F(AB08MDTestColMajor, ZeroDimensionP) {
    // Test with P=0
    int p_zero = 0;
    int LDA = N;
    int LDB = N;
    int LDC = 1;
    int LDD = 1;
    
    std::vector<double> A = A_in;
    std::vector<double> B = B_in;
    std::vector<double> C_empty(1, 0.0);
    std::vector<double> D_empty(1, 0.0);
    
    int rank = -1;
    int expected_zero_rank = 0;
    
    int info = slicot_ab08md(EQUIL, N, M, p_zero, A.data(), LDA, B.data(), LDB,
                             C_empty.data(), LDC, D_empty.data(), LDD, &rank, TOL, ROW_MAJOR);
    
    ASSERT_EQ(info, expected_info);
    ASSERT_EQ(rank, expected_zero_rank);
}

// Test error cases
TEST_F(AB08MDTestColMajor, ErrorCases) {
    int LDA = N;
    int LDB = N;
    int LDC = P;
    int LDD = P;
    
    std::vector<double> A = A_in;
    std::vector<double> B = B_in;
    std::vector<double> C = C_in;
    std::vector<double> D = D_in;
    int rank = -1;
    
    // Test with invalid equil
    int info = slicot_ab08md('X', N, M, P, A.data(), LDA, B.data(), LDB,
                             C.data(), LDC, D.data(), LDD, &rank, TOL, ROW_MAJOR);
    EXPECT_EQ(info, -1);
    
    // Test with negative N
    info = slicot_ab08md(EQUIL, -1, M, P, A.data(), LDA, B.data(), LDB,
                         C.data(), LDC, D.data(), LDD, &rank, TOL, ROW_MAJOR);
    EXPECT_EQ(info, -2);
    
    // Test with negative M
    info = slicot_ab08md(EQUIL, N, -1, P, A.data(), LDA, B.data(), LDB,
                         C.data(), LDC, D.data(), LDD, &rank, TOL, ROW_MAJOR);
    EXPECT_EQ(info, -3);
    
    // Test with negative P
    info = slicot_ab08md(EQUIL, N, M, -1, A.data(), LDA, B.data(), LDB,
                         C.data(), LDC, D.data(), LDD, &rank, TOL, ROW_MAJOR);
    EXPECT_EQ(info, -4);
    
    // Test with invalid LDA (column-major)
    info = slicot_ab08md(EQUIL, N, M, P, A.data(), N-1, B.data(), LDB,
                         C.data(), LDC, D.data(), LDD, &rank, TOL, 0);
    EXPECT_EQ(info, -6);
    
    // Test with invalid LDB (column-major)
    info = slicot_ab08md(EQUIL, N, M, P, A.data(), LDA, B.data(), N-1,
                         C.data(), LDC, D.data(), LDD, &rank, TOL, 0);
    EXPECT_EQ(info, -8);
    
    // Test with invalid LDC (column-major)
    info = slicot_ab08md(EQUIL, N, M, P, A.data(), LDA, B.data(), LDB,
                         C.data(), P-1, D.data(), LDD, &rank, TOL, 0);
    EXPECT_EQ(info, -10);
    
    // Test with invalid LDD (column-major)
    info = slicot_ab08md(EQUIL, N, M, P, A.data(), LDA, B.data(), LDB,
                         C.data(), LDC, D.data(), P-1, &rank, TOL, 0);
    EXPECT_EQ(info, -12);
}
