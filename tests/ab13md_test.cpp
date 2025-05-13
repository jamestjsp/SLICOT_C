#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>

#include "ab13md.h"
#include "slicot_utils.h"
#include "test_utils.h"

// --- Column-Major Test Fixture ---
class AB13MDTestColMajor : public ::testing::Test {
protected:
    // Test parameters
    int N = 6;              // Matrix order
    int M = 5;              // Number of diagonal blocks
    double check_tol = 1e-8; // Tolerance for verification

    // Block structure
    std::vector<int> nblock;  // Block sizes
    std::vector<int> itype;   // Block types (1=real, 2=complex)

    // Input/output data vectors
    std::vector<slicot_complex_double> Z;  // Complex matrix Z
    std::vector<double> X;                 // State vector for reuse
    double bound;                          // Output: computed upper bound
    std::vector<double> D;                 // Output: scaling matrix D
    std::vector<double> G;                 // Output: scaling matrix G

    // Expected results (from example in documentation)
    double expected_bound;
    int expected_info = 0;

    // Result variable
    int info_result = -999;

    // Leading dimension
    int LDZ = 0;

    void SetUp() override {
        // Initialize the block structure
        nblock = {1, 1, 2, 1, 1};  // Block sizes
        itype = {1, 1, 2, 2, 2};   // Block types (1=real, 2=complex)

        // Allocate result arrays
        D.resize(N, 0.0);
        G.resize(N, 0.0);
        X.resize(M + std::count(itype.begin(), itype.end(), 1) - 1, 0.0);

        // Initialize matrix Z with values from the example in the documentation
        // Using std::complex<double> for clarity, will convert to slicot_complex_double
        using std::complex;
        std::vector<complex<double>> Z_complex = {
            // First column
            complex<double>(-1.0, 6.0), complex<double>(3.0, 8.0), 
            complex<double>(4.0, 2.0), complex<double>(-4.0, 11.0),
            complex<double>(5.0, -4.0), complex<double>(-6.0, 14.0),

            // Second column
            complex<double>(2.0, -3.0), complex<double>(-5.0, -9.0),
            complex<double>(-2.0, 5.0), complex<double>(8.0, -7.0),
            complex<double>(-4.0, -8.0), complex<double>(2.0, -5.0),

            // Third column
            complex<double>(3.0, 8.0), complex<double>(-6.0, 2.0),
            complex<double>(-6.0, -7.0), complex<double>(12.0, -1.0),
            complex<double>(1.0, -3.0), complex<double>(4.0, 16.0),

            // Fourth column
            complex<double>(-1.0, 6.0), complex<double>(3.0, 8.0),
            complex<double>(4.0, 2.0), complex<double>(-4.0, 11.0),
            complex<double>(5.0, -4.0), complex<double>(-6.0, 14.0),

            // Fifth column
            complex<double>(2.0, -3.0), complex<double>(-5.0, -9.0),
            complex<double>(-2.0, 5.0), complex<double>(8.0, -7.0),
            complex<double>(-4.0, -8.0), complex<double>(2.0, -5.0),

            // Sixth column
            complex<double>(3.0, 8.0), complex<double>(-6.0, 2.0),
            complex<double>(-6.0, -7.0), complex<double>(12.0, -1.0),
            complex<double>(1.0, -3.0), complex<double>(4.0, 16.0)
        };

        // Convert to slicot_complex_double
        Z.resize(N * N);
        for (size_t i = 0; i < Z_complex.size(); i++) {
            Z[i] = reinterpret_cast<slicot_complex_double&>(Z_complex[i]);
        }

        // Set expected bound from the documentation example
        expected_bound = 41.7475340800; // Approximate value from example

        // Calculate leading dimension for column-major
        LDZ = std::max(1, N);
    }
};

// --- Row-Major Test Fixture ---
class AB13MDTestRowMajor : public AB13MDTestColMajor {
protected:
    // Input data vector in row-major format
    std::vector<slicot_complex_double> Z_rm;

    void SetUp() override {
        // Call the base class SetUp to initialize the column-major data
        AB13MDTestColMajor::SetUp();
        
        // Convert column-major matrix to row-major
        Z_rm.resize(N * N);
        slicot_transpose_to_c(Z.data(), Z_rm.data(), N, N, sizeof(slicot_complex_double));
        
        // In row-major, LDZ is the number of columns
        LDZ = N;
    }
};

// --- Test Cases ---

// Test: Basic functionality (Column-Major)
TEST_F(AB13MDTestColMajor, BasicFunctionality) {
    // First call with fact='N' (no reuse)
    info_result = slicot_ab13md(
        'N', N, Z.data(), LDZ, 
        M, nblock.data(), itype.data(), 
        X.data(), &bound, D.data(), G.data(),
        0 /* column-major */
    );

    // Verify return code
    ASSERT_EQ(info_result, expected_info);
    
    // Verify the computed bound is close to expected
    EXPECT_NEAR(bound, expected_bound, check_tol);
    
    // Output the results for verification
    std::cout << "Structured singular value bound: " << bound << std::endl;

    // Verify that D and G matrices are valid (should be non-zero)
    double d_sum = 0.0, g_sum = 0.0;
    for (int i = 0; i < N; i++) {
        d_sum += std::abs(D[i]);
        g_sum += std::abs(G[i]);
    }
    EXPECT_GT(d_sum, 0.0);
    
    // Test calling again with fact='F' (reuse previous result)
    double bound2;
    std::vector<double> D2(N, 0.0);
    std::vector<double> G2(N, 0.0);

    info_result = slicot_ab13md(
        'F', N, Z.data(), LDZ, 
        M, nblock.data(), itype.data(), 
        X.data(), &bound2, D2.data(), G2.data(),
        0 /* column-major */
    );

    // Verify return code
    ASSERT_EQ(info_result, expected_info);
    
    // Results should be the same or at least close
    EXPECT_NEAR(bound, bound2, check_tol);
}

// Test: Row-Major format
TEST_F(AB13MDTestRowMajor, BasicFunctionality) {
    // Call the wrapper function with row_major=1
    info_result = slicot_ab13md(
        'N', N, Z_rm.data(), LDZ, 
        M, nblock.data(), itype.data(), 
        X.data(), &bound, D.data(), G.data(),
        1 /* row-major */
    );

    // Verify return code
    ASSERT_EQ(info_result, expected_info);
    
    // Verify the results match the expected value
    EXPECT_NEAR(bound, expected_bound, check_tol);
    
    // Output the results for verification
    std::cout << "Row-major structured singular value bound: " << bound << std::endl;
}

// Test: Different Block Structures
TEST_F(AB13MDTestColMajor, DifferentBlockStructures) {
    // Test with a simpler block structure
    std::vector<int> simple_nblock = {2, 4};  // Block sizes
    std::vector<int> simple_itype = {2, 2};   // Block types (both complex)
    
    // Resize X for the new structure
    X.resize(simple_nblock.size() + std::count(simple_itype.begin(), simple_itype.end(), 1) - 1, 0.0);
    
    info_result = slicot_ab13md(
        'N', N, Z.data(), LDZ, 
        static_cast<int>(simple_nblock.size()), simple_nblock.data(), simple_itype.data(), 
        X.data(), &bound, D.data(), G.data(),
        0 /* column-major */
    );

    // Should succeed with different block structure
    EXPECT_EQ(info_result, 0);
    
    // The bound should be different but positive
    EXPECT_GT(bound, 0.0);
    
    std::cout << "Bound with simpler block structure: " << bound << std::endl;
}

// Test: Parameter validation
TEST_F(AB13MDTestColMajor, ParameterValidation) {
    // Test invalid FACT
    info_result = slicot_ab13md(
        'X', N, Z.data(), LDZ, 
        M, nblock.data(), itype.data(), 
        X.data(), &bound, D.data(), G.data(),
        0 /* column-major */
    );
    EXPECT_EQ(info_result, -1);  // First argument had illegal value
    
    // Test negative N
    info_result = slicot_ab13md(
        'N', -1, Z.data(), LDZ, 
        M, nblock.data(), itype.data(), 
        X.data(), &bound, D.data(), G.data(),
        0 /* column-major */
    );
    EXPECT_EQ(info_result, -2);  // N had an illegal value
    
    // Test invalid M
    info_result = slicot_ab13md(
        'N', N, Z.data(), LDZ, 
        0, nblock.data(), itype.data(), 
        X.data(), &bound, D.data(), G.data(),
        0 /* column-major */
    );
    EXPECT_EQ(info_result, -5);  // M must be >= 1
    
    // Test invalid LDZ (column major)
    info_result = slicot_ab13md(
        'N', N, Z.data(), 0, 
        M, nblock.data(), itype.data(), 
        X.data(), &bound, D.data(), G.data(),
        0 /* column-major */
    );
    EXPECT_EQ(info_result, -4);  // LDZ had an illegal value
    
    // Test invalid block type
    std::vector<int> invalid_itype = itype;
    invalid_itype[0] = 3;  // Only 1 or 2 are valid
    
    info_result = slicot_ab13md(
        'N', N, Z.data(), LDZ, 
        M, nblock.data(), invalid_itype.data(), 
        X.data(), &bound, D.data(), G.data(),
        0 /* column-major */
    );
    EXPECT_EQ(info_result, 4);  // Block type must be 1 or 2
    
    // Test invalid real block size (must be 1)
    std::vector<int> invalid_nblock = nblock;
    std::vector<int> real_itype = {1, 1, 1, 1, 1};  // All real blocks
    invalid_nblock[2] = 2;  // Size 2 for a real block
    
    info_result = slicot_ab13md(
        'N', N, Z.data(), LDZ, 
        M, invalid_nblock.data(), real_itype.data(), 
        X.data(), &bound, D.data(), G.data(),
        0 /* column-major */
    );
    EXPECT_EQ(info_result, 3);  // Real block size must be 1
}

// Test: Zero dimension handling
TEST_F(AB13MDTestColMajor, ZeroDimensions) {
    int zero_n = 0;
    
    // Set up parameters for zero matrix
    std::vector<int> zero_nblock = {0};  // Zero block
    std::vector<int> zero_itype = {2};   // Complex type
    std::vector<double> zero_x(1, 0.0);
    std::vector<double> zero_d(1, 0.0);
    std::vector<double> zero_g(1, 0.0);
    double zero_bound = 0.0;
    
    info_result = slicot_ab13md(
        'N', zero_n, nullptr, 1, 
        1, zero_nblock.data(), zero_itype.data(), 
        zero_x.data(), &zero_bound, zero_d.data(), zero_g.data(),
        0 /* column-major */
    );
    
    // Either should succeed with zero dimension or return an appropriate error
    // The FORTRAN routine might have different handling for this case
    if (info_result == 0) {
        EXPECT_EQ(zero_bound, 0.0);  // For zero dimension, bound should be zero
    } else {
        // If it fails, it should be because sum of blocks != N (info=2)
        EXPECT_EQ(info_result, 2);
    }
}
