#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm>

#include "ab13ed.h"
#include "slicot_utils.h"
#include "test_utils.h"

// --- Column-Major Test Fixture ---
class AB13EDTestColMajor : public ::testing::Test {
protected:
    // Test parameters
    int N = 5;              // Matrix order
    double TOL = 9.0;       // Recommended tolerance value
    double check_tol = 1e-4; // Tolerance for test verification

    // Input data vectors (column-major)
    std::vector<double> A;  // Matrix A

    // Expected results
    double expected_low;    // Expected lower bound
    double expected_high;   // Expected upper bound
    int expected_info = 0;  // Expected return code

    // Result variables
    double low;             // Computed lower bound
    double high;            // Computed upper bound
    int info_result = -999; // Initialize to invalid

    // Leading dimension
    int LDA = 0;

    void SetUp() override {
        // Initialize matrix A with values from the example in the documentation
        A = {
            0.1, 0.0, 0.0, 0.0, 0.0,  // First column
            1.0, 0.1, 0.0, 0.0, 0.0,  // Second column
            0.0, 1.0, 0.1, 0.0, 0.0,  // Third column
            0.0, 0.0, 1.0, 0.1, 0.0,  // Fourth column
            0.0, 0.0, 0.0, 1.0, 0.1   // Fifth column
        };

        // Set expected bounds from the example in the documentation
        // These values are approximate and may vary with different implementations
        expected_low  = 2.0e-6;
        expected_high = 2.1e-5;

        // Calculate leading dimension for column-major
        LDA = std::max(1, N);
    }
};

// --- Row-Major Test Fixture ---
class AB13EDTestRowMajor : public AB13EDTestColMajor {
protected:
    // Input data vector in row-major format
    std::vector<double> A_rm;

    void SetUp() override {
        // Call the base class SetUp to initialize the column-major data
        AB13EDTestColMajor::SetUp();
        
        // Convert column-major matrix to row-major
        A_rm.resize(N * N);
        slicot_transpose_to_c(A.data(), A_rm.data(), N, N, sizeof(double));
        
        // In row-major, LDA is the number of columns
        LDA = N;
    }
};

// --- Test Cases ---

// Test: Basic functionality (Column-Major)
TEST_F(AB13EDTestColMajor, BasicFunctionality) {
    // Call the wrapper function
    info_result = slicot_ab13ed(
        N, A.data(), LDA,
        &low, &high, TOL,
        0 /* column-major */
    );

    // Verify return code
    ASSERT_EQ(info_result, expected_info);
    
    // Verify the bounds are within expected ranges
    // Using wider tolerance since exact values depend on machine precision
    EXPECT_GT(low, 0.0);
    EXPECT_GT(high, low);
    EXPECT_LT(low/high, 1.0);

    // Check if values are in the expected range
    EXPECT_NEAR(low, expected_low, 1e-5);
    EXPECT_NEAR(high, expected_high, 1e-4);
    
    // Output the results for verification
    std::cout << "Low bound: " << low << std::endl;
    std::cout << "High bound: " << high << std::endl;
}

// Test: TOL influence (Column-Major)
TEST_F(AB13EDTestColMajor, TolInfluence) {
    double low1, high1, low2, high2;
    double tol1 = 1.0;  // Tighter tolerance
    double tol2 = 20.0; // Looser tolerance
    
    // Call with tighter tolerance
    info_result = slicot_ab13ed(
        N, A.data(), LDA,
        &low1, &high1, tol1,
        0 /* column-major */
    );
    ASSERT_EQ(info_result, expected_info);
    
    // Call with looser tolerance
    info_result = slicot_ab13ed(
        N, A.data(), LDA,
        &low2, &high2, tol2,
        0 /* column-major */
    );
    ASSERT_EQ(info_result, expected_info);
    
    // Check that bounds with tighter tolerance are better (closer)
    EXPECT_LT(high1 - low1, high2 - low2);
    
    std::cout << "TOL1=" << tol1 << ": Low=" << low1 << ", High=" << high1 << std::endl;
    std::cout << "TOL2=" << tol2 << ": Low=" << low2 << ", High=" << high2 << std::endl;
}

// Test: Row-Major format
TEST_F(AB13EDTestRowMajor, BasicFunctionality) {
    // Call the wrapper function with row_major=1
    info_result = slicot_ab13ed(
        N, A_rm.data(), LDA,
        &low, &high, TOL,
        1 /* row-major */
    );

    // Verify return code
    ASSERT_EQ(info_result, expected_info);
    
    // Verify the bounds are within expected ranges
    EXPECT_GT(low, 0.0);
    EXPECT_GT(high, low);
    EXPECT_LT(low/high, 1.0);

    // Check if values are in the expected range - should match column-major
    EXPECT_NEAR(low, expected_low, 1e-5);
    EXPECT_NEAR(high, expected_high, 1e-4);
    
    // Output the results for verification
    std::cout << "Row-major Low bound: " << low << std::endl;
    std::cout << "Row-major High bound: " << high << std::endl;
}

// Test: Zero dimension
TEST_F(AB13EDTestColMajor, ZeroDimension) {
    int zero_n = 0;
    double local_low, local_high;
    
    // Call with N=0
    info_result = slicot_ab13ed(
        zero_n, nullptr, 1,
        &local_low, &local_high, TOL,
        0 /* column-major */
    );
    
    // Should succeed with zero dimension
    EXPECT_EQ(info_result, 0);
    
    // For a zero-dimensional matrix, the bounds should both be zero
    EXPECT_DOUBLE_EQ(local_low, 0.0);
    EXPECT_DOUBLE_EQ(local_high, 0.0);
}

// Test: Different matrix types
TEST_F(AB13EDTestColMajor, DifferentMatrixTypes) {
    // Create a matrix close to instability (eigenvalues near imaginary axis)
    std::vector<double> A_near_unstable = {
        -0.01,  1.0,  0.0,  0.0,  0.0,
        -1.0, -0.01,  0.0,  0.0,  0.0,
         0.0,   0.0, -0.5,  0.0,  0.0,
         0.0,   0.0,  0.0, -0.5,  0.0,
         0.0,   0.0,  0.0,  0.0, -0.5
    };
    
    double near_low, near_high;
    
    // Call with matrix close to instability
    info_result = slicot_ab13ed(
        N, A_near_unstable.data(), LDA,
        &near_low, &near_high, TOL,
        0 /* column-major */
    );
    
    // Should succeed
    EXPECT_EQ(info_result, 0);
    
    // The bounds should be small (matrix close to instability)
    EXPECT_LT(near_low, 0.1);
    
    std::cout << "Near-unstable matrix: Low=" << near_low << ", High=" << near_high << std::endl;
    
    // Compare with original matrix - should have different bounds
    EXPECT_NE(near_low, low);
    EXPECT_NE(near_high, high);
}

// Test: Parameter validation
TEST_F(AB13EDTestColMajor, ParameterValidation) {
    double local_low, local_high;
    
    // Test negative N
    info_result = slicot_ab13ed(
        -1, A.data(), LDA,
        &local_low, &local_high, TOL,
        0 /* column-major */
    );
    EXPECT_EQ(info_result, -1); // N had an illegal value
    
    // Test invalid LDA (column major)
    info_result = slicot_ab13ed(
        N, A.data(), 0,
        &local_low, &local_high, TOL,
        0 /* column-major */
    );
    EXPECT_EQ(info_result, -3); // LDA had an illegal value
    
    // Negative tolerance is allowed, but let's verify it works
    info_result = slicot_ab13ed(
        N, A.data(), LDA,
        &local_low, &local_high, -1.0,
        0 /* column-major */
    );
    EXPECT_EQ(info_result, 0); // Should still succeed with default tolerance
}
