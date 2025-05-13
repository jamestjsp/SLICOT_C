#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>

#include "ab09bd.h"
#include "slicot_utils.h"
#include "test_utils.h"

// --- Column-Major Test Fixture ---
class AB09BDTestColMajor : public ::testing::Test {
protected:
    // Test parameters
    char DICO = 'C';   // Continuous-time system
    char JOB = 'B';    // Balance & Truncate method
    char EQUIL = 'S';  // Scale original state-space
    char ORDSEL = 'F'; // Fixed order
    int N = 4;         // Original system order
    int M = 2;         // Number of inputs
    int P = 3;         // Number of outputs
    int NR = 2;        // Desired reduced order
    double TOL1 = 0.0; // First tolerance (default)
    double TOL2 = 0.0; // Second tolerance (default)
    int IWARN = 0;     // Warning indicator

    // Verification tolerance
    double check_tol = 1e-12;

    // Input/output data vectors (column-major)
    std::vector<double> A; // State matrix
    std::vector<double> B; // Input matrix
    std::vector<double> C; // Output matrix
    std::vector<double> D; // Direct transmission matrix
    std::vector<double> HSV; // Hankel singular values

    // Expected results
    std::vector<double> A_expected; // Reduced A
    std::vector<double> B_expected; // Reduced B
    std::vector<double> C_expected; // Reduced C
    std::vector<double> HSV_expected; // Expected HSVs
    int expected_info = 0;
    int expected_iwarn = 0;

    // Result variables
    int info_result = -999;

    // Leading dimensions
    int LDA = 0;
    int LDB = 0;
    int LDC = 0;
    int LDD = 0;

    void SetUp() override {
        // Initialize system matrices in column-major format
        // Create a stable 4x4 matrix
        A = {
            -1.0, 0.5, 0.0, 0.0,   // First column
             0.5, -2.0, 0.3, 0.0,   // Second column
             0.0, 0.3, -3.0, 0.7,   // Third column
             0.0, 0.0, 0.7, -4.0    // Fourth column
        };
        
        // Input matrix B
        B = {
            1.0, 0.0,   // First column
            0.0, 1.0,   // Second column
            1.0, 0.5,   // Third column
            0.0, 1.0    // Fourth column
        };
        
        // Output matrix C
        C = {
            1.0, 0.0, 0.5, 0.0,   // First row
            0.0, 1.0, 0.0, 0.2,   // Second row
            0.5, 0.3, 0.0, 0.0    // Third row
        };
        
        // Direct transmission matrix D
        D = {
            0.0, 0.0,   // First row
            0.0, 0.0,   // Second row
            0.0, 0.0    // Third row
        };
        
        // Expected results for reduced order model (NR=2)
        // These values would be computed using a reference implementation
        A_expected = {
            -1.5, 0.3,
             0.3, -2.0
        };
        
        B_expected = {
            0.9, 0.1,
            0.1, 0.8
        };
        
        C_expected = {
            0.8, 0.1,
            0.1, 0.9,
            0.4, 0.2
        };
        
        // Allocate HSV vector
        HSV.resize(N, 0.0);

        // Expected Hankel singular values (example values)
        HSV_expected = {2.5, 1.8, 0.4, 0.1};

        // Calculate leading dimensions for column-major
        LDA = std::max(1, N);
        LDB = std::max(1, N);
        LDC = std::max(1, P);
        LDD = std::max(1, P);
    }
};

// --- Row-Major Test Fixture ---
class AB09BDTestRowMajor : public AB09BDTestColMajor {
protected:
    // Input data vectors in row-major format
    std::vector<double> A_rm;
    std::vector<double> B_rm;
    std::vector<double> C_rm;
    std::vector<double> D_rm;

    void SetUp() override {
        // Call the base class SetUp to initialize the column-major data
        AB09BDTestColMajor::SetUp();
        
        // Convert column-major matrices to row-major
        A_rm.resize(N * N);
        B_rm.resize(N * M);
        C_rm.resize(P * N);
        D_rm.resize(P * M);
        
        // Convert A (N x N)
        slicot_transpose_to_c(A.data(), A_rm.data(), N, N, sizeof(double));
        
        // Convert B (N x M)
        slicot_transpose_to_c(B.data(), B_rm.data(), N, M, sizeof(double));
        
        // Convert C (P x N)
        slicot_transpose_to_c(C.data(), C_rm.data(), P, N, sizeof(double));
        
        // Convert D (P x M)
        slicot_transpose_to_c(D.data(), D_rm.data(), P, M, sizeof(double));
        
        // In row-major, the leading dimensions are different (number of columns)
        LDA = N;  // Number of columns of A
        LDB = M;  // Number of columns of B
        LDC = N;  // Number of columns of C
        LDD = M;  // Number of columns of D
    }
};

// --- Test Cases ---

// Test: Fixed order reduction (Column-Major)
TEST_F(AB09BDTestColMajor, FixedOrderReduction) {
    // Make a copy of matrices for testing since they'll be overwritten
    std::vector<double> A_copy = A;
    std::vector<double> B_copy = B;
    std::vector<double> C_copy = C;
    std::vector<double> D_copy = D;
    
    int nr_copy = NR;
    
    // Call the wrapper function - now with both TOL1 and TOL2
    info_result = slicot_ab09bd(
        DICO, JOB, EQUIL, ORDSEL, 
        N, M, P, &nr_copy,
        A_copy.data(), LDA,
        B_copy.data(), LDB,
        C_copy.data(), LDC,
        D_copy.data(), LDD,
        HSV.data(),
        TOL1, TOL2, &IWARN,
        0 /* column-major */
    );

    // Verify return code and NR
    ASSERT_EQ(info_result, expected_info);
    ASSERT_EQ(nr_copy, NR);
    ASSERT_EQ(IWARN, expected_iwarn);

    // Verify Hankel singular values are in decreasing order
    for (int i = 0; i < N-1; i++) {
        EXPECT_GE(HSV[i], HSV[i+1]) << "HSV not in decreasing order at index " << i;
    }
    
    // Verify reduced matrices - Fix indexing for column-major format
    // For A (NR x NR)
    for (int i = 0; i < NR; i++) {
        for (int j = 0; j < NR; j++) {
            EXPECT_NEAR(A_copy[j*LDA + i], A_expected[j*NR + i], check_tol) 
                << "A mismatch at position (" << i << ", " << j << ")";
        }
    }
    
    // For B (NR x M)
    for (int i = 0; i < NR; i++) {
        for (int j = 0; j < M; j++) {
            EXPECT_NEAR(B_copy[j*LDB + i], B_expected[j*NR + i], check_tol)
                << "B mismatch at position (" << i << ", " << j << ")";
        }
    }
    
    // For C (P x NR)
    for (int i = 0; i < P; i++) {
        for (int j = 0; j < NR; j++) {
            EXPECT_NEAR(C_copy[j*LDC + i], C_expected[j*NR + i], check_tol)
                << "C mismatch at position (" << i << ", " << j << ")";
        }
    }
    
    // D is not modified by the function
}

// Test: Auto order selection (Column-Major)
TEST_F(AB09BDTestColMajor, AutoOrderSelection) {
    // Make a copy of matrices for testing
    std::vector<double> A_copy = A;
    std::vector<double> B_copy = B;
    std::vector<double> C_copy = C;
    std::vector<double> D_copy = D;
    
    int nr_copy = N;  // Will be automatically determined
    char ordsel_copy = 'A'; // Automatic order selection
    double tol1_copy = 0.5; // Set tolerance to get specific reduction
    double tol2_copy = 0.0; // Use default for second tolerance
    
    // Call the wrapper function - now with both tolerances
    info_result = slicot_ab09bd(
        DICO, JOB, EQUIL, ordsel_copy,
        N, M, P, &nr_copy,
        A_copy.data(), LDA,
        B_copy.data(), LDB,
        C_copy.data(), LDC,
        D_copy.data(), LDD,
        HSV.data(),
        tol1_copy, tol2_copy, &IWARN,
        0 /* column-major */
    );

    // Verify return code
    EXPECT_EQ(info_result, expected_info);
    
    // Expect reduced order based on HSV and tolerance
    EXPECT_GT(nr_copy, 0);
    EXPECT_LE(nr_copy, N);
}

// Test: Fixed order reduction (Row-Major)
TEST_F(AB09BDTestRowMajor, FixedOrderReduction) {
    // Make a copy of row-major matrices for testing
    std::vector<double> A_rm_copy = A_rm;
    std::vector<double> B_rm_copy = B_rm;
    std::vector<double> C_rm_copy = C_rm;
    std::vector<double> D_rm_copy = D_rm;
    
    int nr_copy = NR;
    
    // Call the wrapper function with row_major=1 - now with both tolerances
    info_result = slicot_ab09bd(
        DICO, JOB, EQUIL, ORDSEL,
        N, M, P, &nr_copy,
        A_rm_copy.data(), LDA,
        B_rm_copy.data(), LDB,
        C_rm_copy.data(), LDC,
        D_rm_copy.data(), LDD,
        HSV.data(),
        TOL1, TOL2, &IWARN,
        1 /* row-major */
    );

    // Verify return code and NR
    ASSERT_EQ(info_result, expected_info);
    ASSERT_EQ(nr_copy, NR);
    ASSERT_EQ(IWARN, expected_iwarn);

    // Verify Hankel singular values
    for (int i = 0; i < N-1; i++) {
        EXPECT_GE(HSV[i], HSV[i+1]) << "HSV not in decreasing order at index " << i;
    }
    
    // Convert expected results to row-major for comparison
    std::vector<double> A_expected_rm(NR * NR);
    std::vector<double> B_expected_rm(NR * M);
    std::vector<double> C_expected_rm(P * NR);
    
    slicot_transpose_to_c(A_expected.data(), A_expected_rm.data(), NR, NR, sizeof(double));
    slicot_transpose_to_c(B_expected.data(), B_expected_rm.data(), NR, M, sizeof(double));
    slicot_transpose_to_c(C_expected.data(), C_expected_rm.data(), P, NR, sizeof(double));
    
    // Verify reduced matrices in row-major format
    // For A (NR x NR)
    for (int i = 0; i < NR; i++) {
        for (int j = 0; j < NR; j++) {
            EXPECT_NEAR(A_rm_copy[i*LDA + j], A_expected_rm[i*NR + j], check_tol) 
                << "A_rm mismatch at position (" << i << ", " << j << ")";
        }
    }
    
    // For B (NR x M)
    for (int i = 0; i < NR; i++) {
        for (int j = 0; j < M; j++) {
            EXPECT_NEAR(B_rm_copy[i*LDB + j], B_expected_rm[i*M + j], check_tol)
                << "B_rm mismatch at position (" << i << ", " << j << ")";
        }
    }
    
    // For C (P x NR)
    for (int i = 0; i < P; i++) {
        for (int j = 0; j < NR; j++) {
            EXPECT_NEAR(C_rm_copy[i*LDC + j], C_expected_rm[i*NR + j], check_tol)
                << "C_rm mismatch at position (" << i << ", " << j << ")";
        }
    }
}

// Test: Zero dimensions
TEST_F(AB09BDTestColMajor, ZeroDimensions) {
    int zero_n = 0;
    int zero_m = 0;
    int zero_p = 0;
    int zero_nr = 0;
    int local_iwarn = 0;  // Need local variable for IWARN
    
    // Call with N=0, M=0, P=0 - now with both tolerances
    info_result = slicot_ab09bd(
        DICO, JOB, EQUIL, ORDSEL,
        zero_n, zero_m, zero_p, &zero_nr,
        nullptr, 1,
        nullptr, 1,
        nullptr, 1,
        nullptr, 1,
        nullptr,
        TOL1, TOL2, &local_iwarn,
        0 /* column-major */
    );
    
    // Should succeed with zero dimensions
    EXPECT_EQ(info_result, 0);
    EXPECT_EQ(zero_nr, 0);
}

// Test: Parameter Validation
TEST_F(AB09BDTestColMajor, ParameterValidation) {
    int nr_copy = NR;
    int local_iwarn = 0;  // Need local variable for IWARN
    
    // Test invalid N - now with both tolerances
    info_result = slicot_ab09bd(
        DICO, JOB, EQUIL, ORDSEL,
        -1, M, P, &nr_copy,
        A.data(), LDA, B.data(), LDB, C.data(), LDC, D.data(), LDD,
        HSV.data(), TOL1, TOL2, &local_iwarn,
        0 /* column-major */
    );
    EXPECT_EQ(info_result, -5); // N had an illegal value
    
    // Update all remaining test calls to include both tolerances
    info_result = slicot_ab09bd(
        'X', JOB, EQUIL, ORDSEL,
        N, M, P, &nr_copy,
        A.data(), LDA, B.data(), LDB, C.data(), LDC, D.data(), LDD,
        HSV.data(), TOL1, TOL2, &local_iwarn,
        0 /* column-major */
    );
    EXPECT_EQ(info_result, -1);
    
    info_result = slicot_ab09bd(
        DICO, 'X', EQUIL, ORDSEL,
        N, M, P, &nr_copy,
        A.data(), LDA, B.data(), LDB, C.data(), LDC, D.data(), LDD,
        HSV.data(), TOL1, TOL2, &local_iwarn,
        0 /* column-major */
    );
    EXPECT_EQ(info_result, -2);
    
    info_result = slicot_ab09bd(
        DICO, JOB, 'X', ORDSEL,
        N, M, P, &nr_copy,
        A.data(), LDA, B.data(), LDB, C.data(), LDC, D.data(), LDD,
        HSV.data(), TOL1, TOL2, &local_iwarn,
        0 /* column-major */
    );
    EXPECT_EQ(info_result, -3);
    
    info_result = slicot_ab09bd(
        DICO, JOB, EQUIL, 'X',
        N, M, P, &nr_copy,
        A.data(), LDA, B.data(), LDB, C.data(), LDC, D.data(), LDD,
        HSV.data(), TOL1, TOL2, &local_iwarn,
        0 /* column-major */
    );
    EXPECT_EQ(info_result, -4);
    
    nr_copy = N + 1;
    info_result = slicot_ab09bd(
        DICO, JOB, EQUIL, 'F',
        N, M, P, &nr_copy,
        A.data(), LDA, B.data(), LDB, C.data(), LDC, D.data(), LDD,
        HSV.data(), TOL1, TOL2, &local_iwarn,
        0 /* column-major */
    );
    EXPECT_EQ(info_result, -8);
    
    info_result = slicot_ab09bd(
        DICO, JOB, EQUIL, ORDSEL,
        N, M, P, &NR,
        A.data(), 0, B.data(), LDB, C.data(), LDC, D.data(), LDD,
        HSV.data(), TOL1, TOL2, &local_iwarn,
        0 /* column-major */
    );
    EXPECT_EQ(info_result, -10);
}
