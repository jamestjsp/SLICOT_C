#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm>

#include "mb03vd.h"
#include "test_utils.h"

// --- Column-Major Test Fixture ---
class MB03VDTestColMajor : public ::testing::Test {
protected:
    // Test parameters from example in documentation
    int N = 4;
    int P = 2;
    int ILO = 1;
    int IHI = 4;
    
    // Verification tolerance - increased to account for numerical precision
    double check_tol = 1e-4;
    
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
    
    // Expected results (from documentation example)
    std::vector<double> A_expected = {
        // H1 (n×n, upper Hessenberg)
        -2.3926, 4.1417, 0.0000, 0.0000,  // First column
        2.7042, -1.7046, -1.6247, 0.0000,  // Second column
        -0.9598, 1.3001, -0.2534, -0.0169,  // Third column
        -1.2335, -1.3120, 1.6453, -0.4451,  // Fourth column
        
        // H2 (n×n, upper triangular)
        -2.5495, 0.0000, 0.0000, 0.0000,  // First column
        2.3402, 1.9725, 0.0000, 0.0000,   // Second column
        4.7021, -0.2483, -0.6290, 0.0000,  // Third column
        0.2329, -2.3493, -0.5975, -0.4426  // Fourth column
    };
    
    // Output Tau (ldtau×p)
    std::vector<double> TAU;
    
    // Result info
    int info_result = -999;
    
    // Leading dimensions
    int LDA1;
    int LDA2;
    int LDTAU;
    
    void SetUp() override {
        // Size output arrays
        TAU.resize((N-1) * P);
        
        // Calculate leading dimensions
        LDA1 = std::max(1, N);
        LDA2 = std::max(1, N);
        LDTAU = std::max(1, N-1);
    }
};

// --- Row-Major Test Fixture ---
class MB03VDTestRowMajor : public ::testing::Test {
protected:
    // Test parameters from example in documentation
    int N = 4;
    int P = 2;
    int ILO = 1;
    int IHI = 4;
    
    // Verification tolerance - increased to account for numerical precision
    double check_tol = 1e-4;
    
    // Input matrices (row-major format)
    // For row-major, we need to transpose the data compared to column-major
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
    
    // Expected results (from documentation example, converted to row-major)
    std::vector<double> A_expected_rm = {
        // H1 (n×n, upper Hessenberg) in row-major format
        -2.3926, 2.7042, -0.9598, -1.2335,  // First row
        4.1417, -1.7046, 1.3001, -1.3120,   // Second row
        0.0000, -1.6247, -0.2534, 1.6453,   // Third row
        0.0000, 0.0000, -0.0169, -0.4451,   // Fourth row
        
        // H2 (n×n, upper triangular) in row-major format
        -2.5495, 2.3402, 4.7021, 0.2329,    // First row
        0.0000, 1.9725, -0.2483, -2.3493,   // Second row
        0.0000, 0.0000, -0.6290, -0.5975,   // Third row
        0.0000, 0.0000, 0.0000, -0.4426     // Fourth row
    };
    
    // Output Tau (ldtau×p)
    std::vector<double> TAU_rm;
    
    // Result info
    int info_result = -999;
    
    // Leading dimensions
    int LDA1;
    int LDA2;
    int LDTAU;
    
    void SetUp() override {
        // Size output arrays
        TAU_rm.resize((N-1) * P);
        
        // Calculate leading dimensions for row-major
        LDA1 = std::max(1, N);  // Number of rows in each 2D slice
        LDA2 = std::max(1, N);  // Stride between columns
        LDTAU = std::max(1, N-1); // Stride between columns in TAU
    }
};

// --- Test Cases ---

// Test: Documentation Example (Column-Major)
TEST_F(MB03VDTestColMajor, DocExample) {
    // Make a copy of the input data since it will be modified
    std::vector<double> A_copy = A;
    
    // Call C wrapper function
    info_result = slicot_mb03vd(N, P, ILO, IHI,
                               A_copy.data(), LDA1, LDA2,
                               TAU.data(), LDTAU,
                               0 /* row_major = false */);
    
    // Verify return code
    ASSERT_EQ(info_result, 0);
    
    // Verify output matrices against expected column-major values
    // For A1 (upper Hessenberg)
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            // Only check upper triangle and first subdiagonal for H1
            if (i <= j+1) {
                EXPECT_NEAR(A_copy[j*N + i], A_expected[j*N + i], check_tol)
                    << "H1 mismatch at position (" << i << "," << j << ")";
            }
        }
    }
    
    // For A2 (upper triangular)
    int offset = N*N;  // Offset to second matrix in 3D array
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            // Only check upper triangle for H2
            if (i <= j) {
                EXPECT_NEAR(A_copy[offset + j*N + i], A_expected[offset + j*N + i], check_tol)
                    << "H2 mismatch at position (" << i << "," << j << ")";
            }
        }
    }
}

// Test: Documentation Example (Row-Major)
TEST_F(MB03VDTestRowMajor, DocExample) {
    // Make a copy of the input data since it will be modified
    std::vector<double> A_rm_copy = A_rm;
    
    // Call C wrapper function with row_major = true
    // Now with proper row-major support
    info_result = slicot_mb03vd(N, P, ILO, IHI,
                               A_rm_copy.data(), LDA1, LDA2,
                               TAU_rm.data(), LDTAU,
                               1 /* row_major = true */);
    
    // Verify return code
    ASSERT_EQ(info_result, 0);
    
    // Verify output directly against row-major expected values
    // For A1 (upper Hessenberg)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            // Only check upper triangle and first subdiagonal for H1
            if (i <= j+1) {
                EXPECT_NEAR(A_rm_copy[i*N + j], A_expected_rm[i*N + j], check_tol)
                    << "H1 mismatch at position (" << i << "," << j << ")";
            }
        }
    }
    
    // For A2 (upper triangular)
    int offset = N*N;  // Offset to second matrix in 3D array
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            // Only check upper triangle for H2
            if (i <= j) {
                EXPECT_NEAR(A_rm_copy[offset + i*N + j], A_expected_rm[offset + i*N + j], check_tol)
                    << "H2 mismatch at position (" << i << "," << j << ")";
            }
        }
    }
}

// Test: Parameter Validation
TEST_F(MB03VDTestColMajor, ParameterValidation) {
    // Use properly sized arrays for all tests to avoid segfaults
    std::vector<double> dummy_A(N*N*P, 0.0); 
    std::vector<double> dummy_tau(std::max(1, N-1) * P, 0.0);
    
    // Test invalid N
    info_result = slicot_mb03vd(-1, P, ILO, IHI,
                               dummy_A.data(), N, N,
                               dummy_tau.data(), N-1, 0);
    EXPECT_EQ(info_result, -1);
    
    // Test invalid P
    info_result = slicot_mb03vd(N, 0, ILO, IHI,
                               dummy_A.data(), N, N,
                               dummy_tau.data(), N-1, 0);
    EXPECT_EQ(info_result, -2);
    
    // Test invalid ILO (< 1)
    info_result = slicot_mb03vd(N, P, 0, IHI,
                               dummy_A.data(), N, N,
                               dummy_tau.data(), N-1, 0);
    EXPECT_EQ(info_result, -3);
    
    // Test invalid IHI (> N)
    info_result = slicot_mb03vd(N, P, ILO, N+1,
                               dummy_A.data(), N, N,
                               dummy_tau.data(), N-1, 0);
    EXPECT_EQ(info_result, -4);
    
    // Test invalid LDA1
    info_result = slicot_mb03vd(N, P, ILO, IHI,
                               dummy_A.data(), 0, N,
                               dummy_tau.data(), N-1, 0);
    EXPECT_EQ(info_result, -6);
    
    // Test invalid LDA2
    info_result = slicot_mb03vd(N, P, ILO, IHI,
                               dummy_A.data(), N, 0,
                               dummy_tau.data(), N-1, 0);
    EXPECT_EQ(info_result, -7);
    
    // Test invalid LDTAU
    info_result = slicot_mb03vd(N, P, ILO, IHI,
                               dummy_A.data(), N, N,
                               dummy_tau.data(), 0, 0);
    EXPECT_EQ(info_result, -9);
    
    // Test NULL pointer for A when N > 0
    info_result = slicot_mb03vd(N, P, ILO, IHI,
                               nullptr, N, N,
                               dummy_tau.data(), N-1, 0);
    EXPECT_EQ(info_result, -5);
    
    // Test NULL pointer for TAU when N > 1
    if (N > 1) {
        info_result = slicot_mb03vd(N, P, ILO, IHI,
                                  dummy_A.data(), N, N,
                                  nullptr, N-1, 0);
        EXPECT_EQ(info_result, -8);
    }
}

// Test: Zero Dimension Case
TEST_F(MB03VDTestColMajor, ZeroDimension) {
    int zero_n = 0;
    std::vector<double> dummy_tau(1);
    
    // Call wrapper with N=0
    // We expect IHI=0 since the constraint is min(ILO,N) <= IHI <= N
    // When N=0, this means IHI should be at most 0
    info_result = slicot_mb03vd(zero_n, P, 1, 0,
                               nullptr, 1, 1,
                               dummy_tau.data(), 1, 0);
    EXPECT_EQ(info_result, 0);
}
