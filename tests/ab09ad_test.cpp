#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>

#include "ab09ad.h"
#include "slicot_utils.h"
#include "test_config.h"

// Using the example from AB09AD documentation
class AB09ADTestColMajor : public ::testing::Test {
protected:
    // Test parameters
    char DICO = 'C';      // Continuous-time system
    char JOB = 'N';       // Use balancing-free square-root method
    char EQUIL = 'N';     // Do not perform equilibration
    char ORDSEL = 'A';    // Order is automatically determined by TOL
    int N = 7;            // Order of the original state-space
    int M = 2;            // Number of inputs
    int P = 3;            // Number of outputs
    int NR = 0;           // Desired reduced order (0 for auto with ORDSEL='A')
    double TOL = 1.0e-1;  // Tolerance for determining the order
    int IWARN = 0;        // Warning indicator

    // Verification tolerance
    double check_tol = 1e-4;

    // Input data vectors (column-major format)
    std::vector<double> A;  // Original state matrix
    std::vector<double> B;  // Original input matrix
    std::vector<double> C;  // Original output matrix
    std::vector<double> HSV; // Hankel singular values

    // Expected results (from documentation)
    std::vector<double> HSV_expected; // Expected Hankel singular values
    std::vector<double> A_expected;   // Expected reduced A matrix
    std::vector<double> B_expected;   // Expected reduced B matrix
    std::vector<double> C_expected;   // Expected reduced C matrix
    int NR_expected = 5;  // Expected reduced order
    int info_expected = 0;

    // Result variable
    int info_result = -999;

    // Leading dimensions
    int LDA = 0;
    int LDB = 0;
    int LDC = 0;

    void SetUp() override {
        // Set leading dimensions for column-major format
        LDA = std::max(1, N);
        LDB = std::max(1, N);
        LDC = std::max(1, P);

        // Initialize matrices based on example from documentation
        A = {
            // Column 1
            -0.04165, -5.2100, 0.0, 0.545, 0.0, 0.0, 0.0,
            // Column 2
            0.0, -12.500, 3.3300, 0.0, 0.0, 0.0, 0.0,
            // Column 3
            4.9200, 0.0, -3.3300, 0.0, 0.0, 0.0, 0.0,
            // Column 4
            -4.9200, 0.0, 0.0, 0.0, 4.9200, 0.0, 0.0,
            // Column 5
            0.0, 0.0, 0.0, -0.5450, -0.04165, -5.2100, 0.0,
            // Column 6
            0.0, 0.0, 0.0, 0.0, 0.0, -12.500, 3.3300,
            // Column 7
            0.0, 0.0, 0.0, 0.0, 4.9200, 0.0, -3.3300
        };

        B = {
            // Column 1
            0.0, 12.5, 0.0, 0.0, 0.0, 0.0, 0.0, 
            // Column 2
            0.0, 0.0, 0.0, 0.0, 0.0, 12.5, 0.0
        };

        C = {
            // Column 1
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            // Column 2
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
            // Column 3
            0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0
        };

        HSV.resize(N);

        // Expected results for reduced system (order 5)
        // Hankel singular values from documentation example
        HSV_expected = {2.5139, 2.0846, 1.9178, 0.7666, 0.5473, 0.0253, 0.0246};
        
        // Expected reduced A matrix (5x5)
        A_expected = {
            1.3451, -4.0214, 0.0, 0.0, 1.2402,
            5.0399, -3.6604, 0.0, 0.0, 1.6416,
            0.0, 0.0, 0.5124, -4.2167, 0.0,
            0.0, 0.0, 1.7910, -2.9900, 0.0,
            4.5315, -0.9056, 0.0, 0.0, -0.0586
        };

        // Expected reduced B matrix (5x2)
        B_expected = {
            -0.3857, 0.3857,
            -3.1753, 3.1753,
            -0.7447, -0.7447,
            -3.6872, -3.6872,
            1.8197, -1.8197
        };

        // Expected reduced C matrix (3x5)
        C_expected = {
            -0.6704, 0.1828, -0.6582, 0.2222, -0.0104,
            0.1089, 0.4867, 0.0, 0.0, 0.8651,
            0.6704, -0.1828, -0.6582, 0.2222, 0.0104
        };
    }
};

// Row-major test fixture derived from the column-major fixture
class AB09ADTestRowMajor : public AB09ADTestColMajor {
protected:
    // Input data vectors (row-major format)
    std::vector<double> A_rm;
    std::vector<double> B_rm;
    std::vector<double> C_rm;
    std::vector<double> HSV_rm; // Will be the same as HSV (1D array)

    // Expected results in row-major format
    std::vector<double> A_expected_rm;
    std::vector<double> B_expected_rm;
    std::vector<double> C_expected_rm;

    void SetUp() override {
        // Call base class SetUp to initialize column-major data
        AB09ADTestColMajor::SetUp();

        // Convert from column-major to row-major
        A_rm.resize(N * N);
        B_rm.resize(N * M);
        C_rm.resize(P * N);
        HSV_rm.resize(N);  // 1D array, no conversion needed

        // Convert A matrix
        slicot_transpose_to_c(A.data(), A_rm.data(), N, N, sizeof(double));
        // Convert B matrix
        slicot_transpose_to_c(B.data(), B_rm.data(), N, M, sizeof(double));
        // Convert C matrix
        slicot_transpose_to_c(C.data(), C_rm.data(), P, N, sizeof(double));

        // Convert expected results to row-major
        A_expected_rm.resize(NR_expected * NR_expected);
        B_expected_rm.resize(NR_expected * M);
        C_expected_rm.resize(P * NR_expected);
        
        slicot_transpose_to_c(A_expected.data(), A_expected_rm.data(), NR_expected, NR_expected, sizeof(double));
        slicot_transpose_to_c(B_expected.data(), B_expected_rm.data(), NR_expected, M, sizeof(double));
        slicot_transpose_to_c(C_expected.data(), C_expected_rm.data(), P, NR_expected, sizeof(double));

        // For row-major, LDA/LDB/LDC refer to columns
        LDA = std::max(1, N);
        LDB = std::max(1, M);
        LDC = std::max(1, N);
    }
};

// Test: Documentation Example (Column-Major)
TEST_F(AB09ADTestColMajor, DocExample) {
    // Make copies of arrays that will be modified
    std::vector<double> A_copy = A;
    std::vector<double> B_copy = B;
    std::vector<double> C_copy = C;
    
    // Set NR to 0 for automatic determination by TOL
    int nr_copy = NR;
    
    // Call C wrapper
    info_result = slicot_ab09ad(DICO, JOB, EQUIL, ORDSEL,
                                N, M, P, &nr_copy,
                                A_copy.data(), LDA,
                                B_copy.data(), LDB,
                                C_copy.data(), LDC,
                                HSV.data(), TOL, &IWARN,
                                0 /* column-major */);

    // Verify return code
    ASSERT_EQ(info_result, info_expected);
    
    // Verify reduced order
    EXPECT_EQ(nr_copy, NR_expected);

    // Verify Hankel singular values
    for (int i = 0; i < N; i++) {
        EXPECT_NEAR(HSV[i], HSV_expected[i], check_tol) << "HSV mismatch at index " << i;
    }

    // Verify reduced A matrix
    for (int i = 0; i < NR_expected * NR_expected; i++) {
        EXPECT_NEAR(A_copy[i], A_expected[i], check_tol) << "A mismatch at index " << i;
    }

    // Verify reduced B matrix
    for (int i = 0; i < NR_expected * M; i++) {
        EXPECT_NEAR(B_copy[i], B_expected[i], check_tol) << "B mismatch at index " << i;
    }

    // Verify reduced C matrix
    for (int i = 0; i < P * NR_expected; i++) {
        EXPECT_NEAR(C_copy[i], C_expected[i], check_tol) << "C mismatch at index " << i;
    }
}

// Test: Documentation Example (Row-Major)
TEST_F(AB09ADTestRowMajor, DocExample) {
    // Make copies of arrays that will be modified
    std::vector<double> A_copy = A_rm;
    std::vector<double> B_copy = B_rm;
    std::vector<double> C_copy = C_rm;
    std::vector<double> HSV_copy = HSV_rm;
    
    // Set NR to 0 for automatic determination by TOL
    int nr_copy = NR;
    
    // Call C wrapper (row-major)
    info_result = slicot_ab09ad(DICO, JOB, EQUIL, ORDSEL,
                                N, M, P, &nr_copy,
                                A_copy.data(), LDA,
                                B_copy.data(), LDB,
                                C_copy.data(), LDC,
                                HSV_copy.data(), TOL, &IWARN,
                                1 /* row-major */);

    // Verify return code
    ASSERT_EQ(info_result, info_expected);
    
    // Verify reduced order
    EXPECT_EQ(nr_copy, NR_expected);

    // Verify Hankel singular values 
    for (int i = 0; i < N; i++) {
        EXPECT_NEAR(HSV_copy[i], HSV_expected[i], check_tol) << "HSV mismatch at index " << i;
    }

    // Verify reduced A matrix (row-major)
    for (int i = 0; i < NR_expected * NR_expected; i++) {
        EXPECT_NEAR(A_copy[i], A_expected_rm[i], check_tol) << "A mismatch at index " << i;
    }

    // Verify reduced B matrix (row-major)
    for (int i = 0; i < NR_expected * M; i++) {
        EXPECT_NEAR(B_copy[i], B_expected_rm[i], check_tol) << "B mismatch at index " << i;
    }

    // Verify reduced C matrix (row-major)
    for (int i = 0; i < P * NR_expected; i++) {
        EXPECT_NEAR(C_copy[i], C_expected_rm[i], check_tol) << "C mismatch at index " << i;
    }
}

// Test: Fixed Order Selection (Column-Major)
TEST_F(AB09ADTestColMajor, FixedOrderSelection) {
    // Make copies of arrays that will be modified
    std::vector<double> A_copy = A;
    std::vector<double> B_copy = B;
    std::vector<double> C_copy = C;
    std::vector<double> HSV_copy = HSV;
    
    // Set NR to fixed value
    int nr_copy = 3; // Ask for order 3
    char ordsel_fixed = 'F'; // Fixed order selection
    
    // Call C wrapper
    info_result = slicot_ab09ad(DICO, JOB, EQUIL, ordsel_fixed,
                                N, M, P, &nr_copy,
                                A_copy.data(), LDA,
                                B_copy.data(), LDB,
                                C_copy.data(), LDC,
                                HSV_copy.data(), TOL, &IWARN,
                                0 /* column-major */);

    // Verify return code
    ASSERT_EQ(info_result, info_expected);
    
    // Verify reduced order is 3 (as requested)
    EXPECT_EQ(nr_copy, 3);

    // Verify Hankel singular values
    for (int i = 0; i < N; i++) {
        EXPECT_NEAR(HSV_copy[i], HSV_expected[i], check_tol) << "HSV mismatch at index " << i;
    }
}

// Test: Zero Dimension Case
TEST_F(AB09ADTestColMajor, ZeroDimension) {
    // Set parameters for zero system
    int n_zero = 0;
    int m_zero = 0;
    int p_zero = 0;
    int nr_zero = 0;
    
    // Call C wrapper with zero-dimensional system
    info_result = slicot_ab09ad(DICO, JOB, EQUIL, ORDSEL,
                               n_zero, m_zero, p_zero, &nr_zero,
                               nullptr, 1, // A
                               nullptr, 1, // B
                               nullptr, 1, // C
                               nullptr, TOL, &IWARN,
                               0 /* column-major */);

    // Should succeed with zero-dimension system
    ASSERT_EQ(info_result, 0);
    EXPECT_EQ(nr_zero, 0);
}

// Test: Parameter Validation
TEST_F(AB09ADTestColMajor, ParameterValidation) {
    // Temporary parameters
    std::vector<double> dummy_A(1);
    std::vector<double> dummy_B(1);
    std::vector<double> dummy_C(1);
    std::vector<double> dummy_HSV(1);
    int dummy_nr = 1;
    int dummy_iwarn = 0;
    
    // Test invalid DICO
    info_result = slicot_ab09ad('X', JOB, EQUIL, ORDSEL, 
                                N, M, P, &dummy_nr,
                                dummy_A.data(), N, 
                                dummy_B.data(), N, 
                                dummy_C.data(), P, 
                                dummy_HSV.data(), TOL, &dummy_iwarn, 
                                0);
    EXPECT_EQ(info_result, -1);
    
    // Test invalid JOB
    info_result = slicot_ab09ad(DICO, 'X', EQUIL, ORDSEL, 
                                N, M, P, &dummy_nr,
                                dummy_A.data(), N, 
                                dummy_B.data(), N, 
                                dummy_C.data(), P, 
                                dummy_HSV.data(), TOL, &dummy_iwarn, 
                                0);
    EXPECT_EQ(info_result, -2);
    
    // Test invalid EQUIL
    info_result = slicot_ab09ad(DICO, JOB, 'X', ORDSEL, 
                                N, M, P, &dummy_nr,
                                dummy_A.data(), N, 
                                dummy_B.data(), N, 
                                dummy_C.data(), P, 
                                dummy_HSV.data(), TOL, &dummy_iwarn, 
                                0);
    EXPECT_EQ(info_result, -3);
    
    // Test invalid ORDSEL
    info_result = slicot_ab09ad(DICO, JOB, EQUIL, 'X', 
                                N, M, P, &dummy_nr,
                                dummy_A.data(), N, 
                                dummy_B.data(), N, 
                                dummy_C.data(), P, 
                                dummy_HSV.data(), TOL, &dummy_iwarn, 
                                0);
    EXPECT_EQ(info_result, -4);
    
    // Test invalid N
    info_result = slicot_ab09ad(DICO, JOB, EQUIL, ORDSEL, 
                                -1, M, P, &dummy_nr,
                                dummy_A.data(), N, 
                                dummy_B.data(), N, 
                                dummy_C.data(), P, 
                                dummy_HSV.data(), TOL, &dummy_iwarn, 
                                0);
    EXPECT_EQ(info_result, -5);
    
    // Test invalid M
    info_result = slicot_ab09ad(DICO, JOB, EQUIL, ORDSEL, 
                                N, -1, P, &dummy_nr,
                                dummy_A.data(), N, 
                                dummy_B.data(), N, 
                                dummy_C.data(), P, 
                                dummy_HSV.data(), TOL, &dummy_iwarn, 
                                0);
    EXPECT_EQ(info_result, -6);
    
    // Test invalid P
    info_result = slicot_ab09ad(DICO, JOB, EQUIL, ORDSEL, 
                                N, M, -1, &dummy_nr,
                                dummy_A.data(), N, 
                                dummy_B.data(), N, 
                                dummy_C.data(), P, 
                                dummy_HSV.data(), TOL, &dummy_iwarn, 
                                0);
    EXPECT_EQ(info_result, -7);
    
    // Test invalid fixed NR
    dummy_nr = -1;
    info_result = slicot_ab09ad(DICO, JOB, EQUIL, 'F', 
                                N, M, P, &dummy_nr,
                                dummy_A.data(), N, 
                                dummy_B.data(), N, 
                                dummy_C.data(), P, 
                                dummy_HSV.data(), TOL, &dummy_iwarn, 
                                0);
    EXPECT_EQ(info_result, -8);
    
    // Test missing A matrix
    info_result = slicot_ab09ad(DICO, JOB, EQUIL, ORDSEL, 
                                N, M, P, &dummy_nr,
                                nullptr, N, 
                                dummy_B.data(), N, 
                                dummy_C.data(), P, 
                                dummy_HSV.data(), TOL, &dummy_iwarn, 
                                0);
    EXPECT_EQ(info_result, -9 /* For A */);
    
    // Test invalid LDA
    info_result = slicot_ab09ad(DICO, JOB, EQUIL, ORDSEL, 
                                N, M, P, &dummy_nr,
                                dummy_A.data(), 0 /* invalid LDA */, 
                                dummy_B.data(), N, 
                                dummy_C.data(), P, 
                                dummy_HSV.data(), TOL, &dummy_iwarn, 
                                0);
    EXPECT_EQ(info_result, -10);
}
