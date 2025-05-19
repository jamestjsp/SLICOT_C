#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm> // For std::min, std::max
#include <iostream>  // For debugging
#include <iomanip>   // For std::setprecision

#include "mb03vd.h"
#include "mb03vy.h"
#include "slicot_utils.h" // For MAX, MIN if not in <algorithm> for C++
#include "test_utils.h"

// Helper function to get 3D array index for column-major storage
// Fortran: A(row, col, page) (1-indexed)
// C: a_flat[ (page_idx * lda2 * lda1) + (col_idx * lda1) + row_idx ] (0-indexed)
size_t get_idx(int r, int c, int pg, int lda1, int lda2) {
    return static_cast<size_t>(pg) * lda1 * lda2 +
           static_cast<size_t>(c) * lda1 +
           static_cast<size_t>(r);
}

// --- Column-Major Test Fixture ---
class MB03VYTestColMajor : public ::testing::Test {
protected:
    // Test parameters from example
    int N = 4;
    int P = 2;
    int ILO = 1;
    int IHI = 4;
    
    // Verification tolerance
    double check_tol = 1e-10;
    
    // Input matrices (column-major)
    std::vector<double> A = {
        // A1 (n×n)
        1.5, 1.0, 1.5, 1.0,    // First column
        -0.7, 0.0, -0.7, 0.0,  // Second column
        3.5, 2.0, 2.5, 2.0,    // Third column
        -0.7, 3.0, -0.3, 1.0,  // Fourth column
        
        // A2 (n×n)
        1.5, 1.0, 1.5, 1.0,    // First column
        -0.7, 0.0, -0.7, 0.0,  // Second column
        3.5, 2.0, 2.5, 2.0,    // Third column
        -0.7, 3.0, -0.3, 1.0   // Fourth column
    };
    
    // Storage for TAU from MB03VD
    std::vector<double> TAU;
    
    // Storage for matrices from MB03VD
    std::vector<double> A_MB03VD;
    
    // Result info
    int info_result = -999;
    
    // Leading dimensions
    int LDA1;
    int LDA2;
    int LDTAU;
    
    void SetUp() override {
        // Size output arrays
        TAU.resize((N-1) * P);
        A_MB03VD = A; // Copy for MB03VD input
        
        // Calculate leading dimensions
        LDA1 = std::max(1, N);
        LDA2 = std::max(1, N);
        LDTAU = std::max(1, N-1);
        
        // Run MB03VD first to get the elementary reflectors
        info_result = slicot_mb03vd(N, P, ILO, IHI,
                                  A_MB03VD.data(), LDA1, LDA2,
                                  TAU.data(), LDTAU,
                                  0 /* row_major = false */);
        ASSERT_EQ(info_result, 0);
    }
};

// --- Row-Major Test Fixture ---
class MB03VYTestRowMajor : public ::testing::Test {
protected:
    // Test parameters from example
    int N = 4;
    int P = 2;
    int ILO = 1;
    int IHI = 4;
    
    // Verification tolerance
    double check_tol = 1e-10;
    
    // Input matrices (row-major format)
    std::vector<double> A_rm = {
        // A1 (n×n) in row-major format
        1.5, -0.7, 3.5, -0.7,  // First row
        1.0, 0.0, 2.0, 3.0,    // Second row
        1.5, -0.7, 2.5, -0.3,  // Third row
        1.0, 0.0, 2.0, 1.0,    // Fourth row
        
        // A2 (n×n) in row-major format
        1.5, -0.7, 3.5, -0.7,  // First row
        1.0, 0.0, 2.0, 3.0,    // Second row
        1.5, -0.7, 2.5, -0.3,  // Third row
        1.0, 0.0, 2.0, 1.0     // Fourth row
    };
    
    // Storage for TAU from MB03VD
    std::vector<double> TAU_rm;
    
    // Storage for matrices from MB03VD
    std::vector<double> A_MB03VD_rm;
    
    // Result info
    int info_result = -999;
    
    // Leading dimensions
    int LDA1;
    int LDA2;
    int LDTAU;
    
    void SetUp() override {
        // Size output arrays
        TAU_rm.resize((N-1) * P);
        A_MB03VD_rm = A_rm; // Copy for MB03VD input
        
        // Calculate leading dimensions for row-major
        LDA1 = std::max(1, N);  // Number of rows in each 2D slice
        LDA2 = std::max(1, N);  // Number of columns in 2D slice
        LDTAU = std::max(1, N-1); // Stride between columns in TAU
        
        // Run MB03VD first to get the elementary reflectors
        info_result = slicot_mb03vd(N, P, ILO, IHI,
                                  A_MB03VD_rm.data(), LDA1, LDA2,
                                  TAU_rm.data(), LDTAU,
                                  1 /* row_major = true */);
        ASSERT_EQ(info_result, 0);
    }
};

// --- Test Cases ---

// Test: Generate Q matrices (Column-Major)
TEST_F(MB03VYTestColMajor, GenerateQMatrices) {
    // Copy the output of MB03VD as input for MB03VY
    std::vector<double> A_MB03VY = A_MB03VD;
    
    // Call MB03VY to generate the orthogonal matrices Q
    info_result = slicot_mb03vy(N, P, ILO, IHI,
                              A_MB03VY.data(), LDA1, LDA2,
                              TAU.data(), LDTAU,
                              0 /* row_major = false */);
    
    // Verify return code
    ASSERT_EQ(info_result, 0);
    
    // Verify orthogonality property for each Q matrix by computing Q'*Q = I
    for (int k = 0; k < P; ++k) {
        double* Q = &A_MB03VY[k * N * N];
        std::vector<double> QTQ(N * N, 0.0); // Initialize Q'*Q to zeros
        
        // Compute Q'*Q using matrix multiplication
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                double sum = 0.0;
                for (int l = 0; l < N; ++l) {
                    sum += Q[l * N + i] * Q[l * N + j]; // Q[l,i] * Q[l,j]
                }
                QTQ[i * N + j] = sum;
            }
        }
        
        // Check that Q'*Q is close to the identity matrix
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (i == j) {
                    EXPECT_NEAR(QTQ[i * N + j], 1.0, check_tol) 
                        << "Non-identity diagonal in Q'*Q at (" << i << "," << j << ") for k=" << k;
                } else {
                    EXPECT_NEAR(QTQ[i * N + j], 0.0, check_tol)
                        << "Non-zero off-diagonal in Q'*Q at (" << i << "," << j << ") for k=" << k;
                }
            }
        }
    }
}

// Test: Generate Q matrices (Row-Major)
TEST_F(MB03VYTestRowMajor, GenerateQMatrices) {
    // Copy the output of MB03VD as input for MB03VY
    std::vector<double> A_MB03VY_rm = A_MB03VD_rm;
    
    // Call MB03VY to generate the orthogonal matrices Q with row_major=true
    info_result = slicot_mb03vy(N, P, ILO, IHI,
                              A_MB03VY_rm.data(), LDA1, LDA2,
                              TAU_rm.data(), LDTAU,
                              1 /* row_major = true */);
    
    // Verify return code
    ASSERT_EQ(info_result, 0);
    
    // Verify orthogonality property for each Q matrix by computing Q*Q' = I
    // Note: In row-major, we're computing Q*Q' since matrices are transposed compared to column-major
    for (int k = 0; k < P; ++k) {
        double* Q = &A_MB03VY_rm[k * N * N];
        std::vector<double> QQT(N * N, 0.0); // Initialize Q*Q' to zeros
        
        // Compute Q*Q' using matrix multiplication
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                double sum = 0.0;
                for (int l = 0; l < N; ++l) {
                    sum += Q[i * N + l] * Q[j * N + l]; // Q[i,l] * Q[j,l]
                }
                QQT[i * N + j] = sum;
            }
        }
        
        // Check that Q*Q' is close to the identity matrix
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (i == j) {
                    EXPECT_NEAR(QQT[i * N + j], 1.0, check_tol) 
                        << "Non-identity diagonal in Q*Q' at (" << i << "," << j << ") for k=" << k;
                } else {
                    EXPECT_NEAR(QQT[i * N + j], 0.0, check_tol)
                        << "Non-zero off-diagonal in Q*Q' at (" << i << "," << j << ") for k=" << k;
                }
            }
        }
    }
}

// Test: Parameter Validation
TEST_F(MB03VYTestColMajor, ParameterValidation) {
    // Use properly sized arrays for all tests to avoid segfaults
    std::vector<double> dummy_A(N*N*P, 0.0); 
    std::vector<double> dummy_tau(std::max(1, N-1) * P, 0.0);
    
    // Test invalid N
    info_result = slicot_mb03vy(-1, P, ILO, IHI,
                               dummy_A.data(), N, N,
                               dummy_tau.data(), N-1, 0);
    EXPECT_EQ(info_result, -1);
    
    // Test invalid P
    info_result = slicot_mb03vy(N, 0, ILO, IHI,
                               dummy_A.data(), N, N,
                               dummy_tau.data(), N-1, 0);
    EXPECT_EQ(info_result, -2);
    
    // Test invalid ILO (< 1)
    info_result = slicot_mb03vy(N, P, 0, IHI,
                               dummy_A.data(), N, N,
                               dummy_tau.data(), N-1, 0);
    EXPECT_EQ(info_result, -3);
    
    // Test invalid IHI (> N)
    info_result = slicot_mb03vy(N, P, ILO, N+1,
                               dummy_A.data(), N, N,
                               dummy_tau.data(), N-1, 0);
    EXPECT_EQ(info_result, -4);
    
    // Test invalid LDA1
    info_result = slicot_mb03vy(N, P, ILO, IHI,
                               dummy_A.data(), 0, N,
                               dummy_tau.data(), N-1, 0);
    EXPECT_EQ(info_result, -6);
    
    // Test invalid LDA2
    info_result = slicot_mb03vy(N, P, ILO, IHI,
                               dummy_A.data(), N, 0,
                               dummy_tau.data(), N-1, 0);
    EXPECT_EQ(info_result, -7);
    
    // Test invalid LDTAU
    info_result = slicot_mb03vy(N, P, ILO, IHI,
                               dummy_A.data(), N, N,
                               dummy_tau.data(), 0, 0);
    EXPECT_EQ(info_result, -9);
    
    // Test NULL pointer for A when N > 0
    info_result = slicot_mb03vy(N, P, ILO, IHI,
                               nullptr, N, N,
                               dummy_tau.data(), N-1, 0);
    EXPECT_EQ(info_result, -5);
    
    // Test NULL pointer for TAU when N > 1 and IHI > ILO
    if (N > 1 && IHI > ILO) {
        info_result = slicot_mb03vy(N, P, ILO, IHI,
                                  dummy_A.data(), N, N,
                                  nullptr, N-1, 0);
        EXPECT_EQ(info_result, -8);
    }
}

// Test: Zero Dimension Case
TEST_F(MB03VYTestColMajor, ZeroDimension) {
    int zero_n = 0;
    std::vector<double> dummy_tau(1);
    
    // Call wrapper with N=0
    // We expect IHI=0 since the constraint is min(ILO,N) <= IHI <= N
    // When N=0, this means IHI should be at most 0
    info_result = slicot_mb03vy(zero_n, P, 1, 0,
                               nullptr, 1, 1,
                               dummy_tau.data(), 1, 0);
    EXPECT_EQ(info_result, 0);
}

