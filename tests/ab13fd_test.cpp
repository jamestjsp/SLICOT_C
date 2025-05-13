#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm>

#include "ab13fd.h"
#include "slicot_utils.h"
#include "test_utils.h"

// --- Column-Major Test Fixture ---
class AB13FDTestColMajor : public ::testing::Test {
protected:
    // Test parameters
    int N = 4;              // Matrix order
    double TOL = 0.0;       // Default tolerance 
    double check_tol = 1e-4; // Tolerance for test verification

    // Input data vectors (column-major)
    std::vector<double> A;  // Matrix A

    // Expected results from documentation example
    double expected_beta;   // Expected distance beta(A)
    double expected_omega;  // Expected frequency omega
    int expected_info = 0;  // Expected return code

    // Result variables
    double beta;            // Computed distance
    double omega;           // Computed frequency
    int info_result = -999; // Initialize to invalid

    // Leading dimension
    int LDA = 0;

    void SetUp() override {
        // Initialize matrix A with values from the example in the documentation
        A = {
            246.500,  -252.500,  -302.500,  -307.500,  // First column
            242.500,  -248.500,  -297.500,  -302.500,  // Second column
            202.500,  -207.500,  -248.500,  -252.500,  // Third column
           -197.500,   202.500,   242.500,   246.500   // Fourth column
        };

        // Set expected results from the example in the documentation
        expected_beta = 0.0039196472317;
        expected_omega = 0.98966520430;

        // Calculate leading dimension for column-major
        LDA = std::max(1, N);
    }
};

// --- Row-Major Test Fixture ---
class AB13FDTestRowMajor : public AB13FDTestColMajor {
protected:
    // Input data vector in row-major format
    std::vector<double> A_rm;

    void SetUp() override {
        // Call the base class SetUp to initialize the column-major data
        AB13FDTestColMajor::SetUp();
        
        // Convert column-major matrix to row-major
        A_rm.resize(N * N);
        slicot_transpose_to_c(A.data(), A_rm.data(), N, N, sizeof(double));
        
        // In row-major, LDA is the number of columns
        LDA = N;
    }
};

// --- Test Cases ---

// Test: Basic functionality (Column-Major)
TEST_F(AB13FDTestColMajor, BasicFunctionality) {
    // Call the wrapper function
    info_result = slicot_ab13fd(
        N, A.data(), LDA,
        &beta, &omega, TOL,
        0 /* column-major */
    );

    // Verify return code
    ASSERT_EQ(info_result, expected_info);
    
    // Verify the computed values match expected values
    EXPECT_NEAR(beta, expected_beta, 1e-4);
    EXPECT_NEAR(omega, expected_omega, 1e-4);
    
    // Output the results for verification
    std::cout << "Beta (stability radius): " << beta << std::endl;
    std::cout << "Omega (frequency): " << omega << std::endl;
}

// Test: TOL influence (Column-Major)
TEST_F(AB13FDTestColMajor, TolInfluence) {
    double beta1, beta2, omega1, omega2;
    double tol1 = 1e-6;  // Tighter tolerance
    double tol2 = 1e-1;  // Looser tolerance
    
    // Call with tighter tolerance
    info_result = slicot_ab13fd(
        N, A.data(), LDA,
        &beta1, &omega1, tol1,
        0 /* column-major */
    );
    ASSERT_EQ(info_result, expected_info);
    
    // Call with looser tolerance
    info_result = slicot_ab13fd(
        N, A.data(), LDA,
        &beta2, &omega2, tol2,
        0 /* column-major */
    );
    ASSERT_EQ(info_result, expected_info);
    
    // Results should be close but may differ slightly due to convergence criteria
    std::cout << "TOL1=" << tol1 << ": beta=" << beta1 << ", omega=" << omega1 << std::endl;
    std::cout << "TOL2=" << tol2 << ": beta=" << beta2 << ", omega=" << omega2 << std::endl;
    
    // Both results should be close to the expected value
    EXPECT_NEAR(beta1, expected_beta, 1e-3);
    EXPECT_NEAR(beta2, expected_beta, 1e-2);
}

// Test: Row-Major format
TEST_F(AB13FDTestRowMajor, BasicFunctionality) {
    // Call the wrapper function with row_major=1
    info_result = slicot_ab13fd(
        N, A_rm.data(), LDA,
        &beta, &omega, TOL,
        1 /* row-major */
    );

    // Verify return code
    ASSERT_EQ(info_result, expected_info);
    
    // Verify the results match the column-major case
    EXPECT_NEAR(beta, expected_beta, 1e-4);
    EXPECT_NEAR(omega, expected_omega, 1e-4);
    
    // Output the results for verification
    std::cout << "Row-major Beta: " << beta << std::endl;
    std::cout << "Row-major Omega: " << omega << std::endl;
}

// Test: Zero dimension
TEST_F(AB13FDTestColMajor, ZeroDimension) {
    int zero_n = 0;
    double local_beta, local_omega;
    
    // Call with N=0
    info_result = slicot_ab13fd(
        zero_n, nullptr, 1,
        &local_beta, &local_omega, TOL,
        0 /* column-major */
    );
    
    // Should succeed with zero dimension
    EXPECT_EQ(info_result, 0);
    
    // For a zero-dimensional matrix, both values should be zero
    EXPECT_DOUBLE_EQ(local_beta, 0.0);
}

// Test: Different Matrix Types
TEST_F(AB13FDTestColMajor, DifferentMatrixTypes) {
    // Matrix with eigenvalues close to imaginary axis
    std::vector<double> near_unstable = {
        -0.01,  1.0,  0.0,  0.0,
        -1.0, -0.01,  0.0,  0.0,
         0.0,   0.0, -0.5,  0.0,
         0.0,   0.0,  0.0, -0.5
    };
    
    double unstable_beta, unstable_omega;
    
    // Call with matrix that's close to unstable
    info_result = slicot_ab13fd(
        N, near_unstable.data(), LDA,
        &unstable_beta, &unstable_omega, TOL,
        0 /* column-major */
    );
    
    // Should succeed 
    EXPECT_EQ(info_result, 0);
    
    // Should have a small beta
    EXPECT_LT(unstable_beta, 0.5);
    
    std::cout << "Near-unstable matrix: beta=" << unstable_beta
              << ", omega=" << unstable_omega << std::endl;
}

// Test: Parameter validation
TEST_F(AB13FDTestColMajor, ParameterValidation) {
    double local_beta, local_omega;
    
    // Test negative N
    info_result = slicot_ab13fd(
        -1, A.data(), LDA,
        &local_beta, &local_omega, TOL,
        0 /* column-major */
    );
    EXPECT_EQ(info_result, -1); // N had an illegal value
    
    // Test invalid LDA (column major)
    info_result = slicot_ab13fd(
        N, A.data(), 0,
        &local_beta, &local_omega, TOL,
        0 /* column-major */
    );
    EXPECT_EQ(info_result, -3); // LDA had an illegal value
    
    // Test negative tolerance (should still work but use default)
    info_result = slicot_ab13fd(
        N, A.data(), LDA,
        &local_beta, &local_omega, -1.0,
        0 /* column-major */
    );
    EXPECT_EQ(info_result, -6); // TOL had an illegal value or
    // EXPECT_EQ(info_result, 0);  // Should still succeed with default tolerance
}
