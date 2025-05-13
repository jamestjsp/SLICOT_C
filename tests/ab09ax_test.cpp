#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>

#include "ab09ax.h"
#include "slicot_utils.h"
#include "test_utils.h"

// --- Column-Major Test Fixture ---
class AB09AXTestColMajor : public ::testing::Test {
protected:
    // Test parameters
    char DICO = 'C';   // Continuous-time system
    char JOB = 'B';    // Balance & Truncate method
    char ORDSEL = 'F'; // Fixed order
    int N = 4;         // Original system order
    int M = 2;         // Number of inputs
    int P = 3;         // Number of outputs
    int NR = 2;        // Desired reduced order
    double TOL = 0.0;  // Tolerance (default)
    int IWARN = 0;     // Warning indicator

    // Verification tolerance
    double check_tol = 1e-12;

    // Input/output data vectors (column-major)
    std::vector<double> A; // State matrix (must be in Schur form)
    std::vector<double> B; // Input matrix
    std::vector<double> C; // Output matrix
    std::vector<double> HSV; // Hankel singular values
    std::vector<double> T;   // Right truncation matrix
    std::vector<double> TI;  // Left truncation matrix

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
    int LDT = 0;
    int LDTI = 0;

    void SetUp() override {
        // Initialize system matrices in column-major format
        // Create a stable 4x4 matrix in Schur form (upper quasi-triangular)
        A = {
            -1.0, 1.0, 0.0, 0.0,   // First column
             0.0, -2.0, 0.0, 0.0,   // Second column
             0.0, 0.0, -3.0, 1.0,   // Third column
             0.0, 0.0, 0.0, -4.0    // Fourth column
        };
        
        // Input matrix B
        B = {
            1.0, 0.0,   // First column
            0.0, 1.0,   // Second column
            1.0, 1.0,   // Third column
            0.0, 1.0    // Fourth column
        };
        
        // Output matrix C
        C = {
            1.0, 0.0, 1.0, 0.0,   // First row
            0.0, 1.0, 0.0, 1.0,   // Second row
            1.0, 1.0, 0.0, 0.0    // Third row
        };
        
        // Expected results for reduced order model (NR=2)
        // These should be calculated based on the actual algorithm output
        // For now using placeholders - real test would use verified reduced model
        A_expected = {
            -1.0, 1.0,
             0.0, -2.0
        };
        
        B_expected = {
            1.0, 0.0,
            0.0, 1.0
        };
        
        C_expected = {
            1.0, 0.0,
            0.0, 1.0,
            1.0, 1.0
        };
        
        // Allocate HSV, T, TI vectors
        HSV.resize(N, 0.0);
        T.resize(N * N, 0.0);      // N x N
        TI.resize(N * N, 0.0);     // N x N

        // Expected Hankel singular values (example values)
        HSV_expected = {3.0, 2.0, 0.5, 0.1};  // Should be computed from the actual algorithm

        // Calculate leading dimensions for column-major
        LDA = std::max(1, N);
        LDB = std::max(1, N);
        LDC = std::max(1, P);
        LDT = std::max(1, N);
        LDTI = std::max(1, N);
    }
};

// --- Row-Major Test Fixture ---
class AB09AXTestRowMajor : public AB09AXTestColMajor {
protected:
    // Input data vectors in row-major format
    std::vector<double> A_rm;
    std::vector<double> B_rm;
    std::vector<double> C_rm;
    std::vector<double> T_rm;
    std::vector<double> TI_rm;

    void SetUp() override {
        // Call the base class SetUp to initialize the column-major data
        AB09AXTestColMajor::SetUp();
        
        // Convert column-major matrices to row-major
        A_rm.resize(N * N);
        B_rm.resize(N * M);
        C_rm.resize(P * N);
        T_rm.resize(N * N);
        TI_rm.resize(N * N);
        
        // Convert A (N x N)
        slicot_transpose_to_c(A.data(), A_rm.data(), N, N, sizeof(double));
        
        // Convert B (N x M)
        slicot_transpose_to_c(B.data(), B_rm.data(), N, M, sizeof(double));
        
        // Convert C (P x N)
        slicot_transpose_to_c(C.data(), C_rm.data(), P, N, sizeof(double));
        
        // In row-major, the leading dimensions are different (number of columns)
        LDA = N;  // Number of columns of A
        LDB = M;  // Number of columns of B
        LDC = N;  // Number of columns of C
        LDT = N;  // Number of columns of T
        LDTI = N; // Number of columns of TI
    }
};

// --- Test Cases ---

// Test: Fixed order reduction (Column-Major)
TEST_F(AB09AXTestColMajor, FixedOrderReduction) {
    // Make a copy of A, B, C for testing since they'll be overwritten
    std::vector<double> A_copy = A;
    std::vector<double> B_copy = B;
    std::vector<double> C_copy = C;
    
    int nr_copy = NR;
    
    // Call the wrapper function
    info_result = slicot_ab09ax(
        DICO, JOB, ORDSEL, N, M, P, &nr_copy,
        A_copy.data(), LDA,
        B_copy.data(), LDB,
        C_copy.data(), LDC,
        HSV.data(),
        T.data(), LDT,
        TI.data(), LDTI,
        TOL, &IWARN,
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
    
    // Verify reduced matrices (first NR columns and rows where applicable)
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
}

// Test: Auto order selection (Column-Major)
TEST_F(AB09AXTestColMajor, AutoOrderSelection) {
    // Make a copy of A, B, C for testing
    std::vector<double> A_copy = A;
    std::vector<double> B_copy = B;
    std::vector<double> C_copy = C;
    
    int nr_copy = N;  // Will be automatically determined
    char ordsel_copy = 'A'; // Automatic order selection
    double tol_copy = 1.0; // Set tolerance to get specific reduction
    
    // Call the wrapper function
    info_result = slicot_ab09ax(
        DICO, JOB, ordsel_copy, N, M, P, &nr_copy,
        A_copy.data(), LDA,
        B_copy.data(), LDB,
        C_copy.data(), LDC,
        HSV.data(),
        T.data(), LDT,
        TI.data(), LDTI,
        tol_copy, &IWARN,
        0 /* column-major */
    );

    // Verify return code
    EXPECT_EQ(info_result, expected_info);
    // Expect reduced order based on HSV and tolerance
    EXPECT_GT(nr_copy, 0);
    EXPECT_LE(nr_copy, N);
}

// Test: Fixed order reduction (Row-Major)
TEST_F(AB09AXTestRowMajor, FixedOrderReduction) {
    // Make a copy of row-major A, B, C for testing
    std::vector<double> A_rm_copy = A_rm;
    std::vector<double> B_rm_copy = B_rm;
    std::vector<double> C_rm_copy = C_rm;
    
    int nr_copy = NR;
    
    // Call the wrapper function with row_major=1
    info_result = slicot_ab09ax(
        DICO, JOB, ORDSEL, N, M, P, &nr_copy,
        A_rm_copy.data(), LDA,
        B_rm_copy.data(), LDB,
        C_rm_copy.data(), LDC,
        HSV.data(), // HSV is 1D and not affected by row/column major
        T_rm.data(), LDT,
        TI_rm.data(), LDTI,
        TOL, &IWARN,
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
TEST_F(AB09AXTestColMajor, ZeroDimensions) {
    int zero_n = 0;
    int zero_m = 0;
    int zero_p = 0;
    int zero_nr = 0;
    
    // Call with N=0, M=0, P=0
    info_result = slicot_ab09ax(
        DICO, JOB, ORDSEL, zero_n, zero_m, zero_p, &zero_nr,
        nullptr, 1, // A with LD=1
        nullptr, 1, // B with LD=1
        nullptr, 1, // C with LD=1
        nullptr,    // No HSV needed
        nullptr, 1, // T with LD=1
        nullptr, 1, // TI with LD=1
        TOL, &IWARN,
        0 /* column-major */
    );
    
    // Should succeed with zero dimensions
    EXPECT_EQ(info_result, 0);
    EXPECT_EQ(zero_nr, 0);
}

// Test: Parameter Validation
TEST_F(AB09AXTestColMajor, ParameterValidation) {
    int nr_copy = NR;
    
    // Test invalid N
    info_result = slicot_ab09ax(
        DICO, JOB, ORDSEL, -1, M, P, &nr_copy,
        A.data(), LDA, B.data(), LDB, C.data(), LDC,
        HSV.data(), T.data(), LDT, TI.data(), LDTI, TOL, &IWARN,
        0 /* column-major */
    );
    EXPECT_EQ(info_result, -4); // N had an illegal value
    
    // Test invalid DICO
    info_result = slicot_ab09ax(
        'X', JOB, ORDSEL, N, M, P, &nr_copy,
        A.data(), LDA, B.data(), LDB, C.data(), LDC,
        HSV.data(), T.data(), LDT, TI.data(), LDTI, TOL, &IWARN,
        0 /* column-major */
    );
    EXPECT_EQ(info_result, -1); // DICO had an illegal value
    
    // Test invalid JOB
    info_result = slicot_ab09ax(
        DICO, 'X', ORDSEL, N, M, P, &nr_copy,
        A.data(), LDA, B.data(), LDB, C.data(), LDC,
        HSV.data(), T.data(), LDT, TI.data(), LDTI, TOL, &IWARN,
        0 /* column-major */
    );
    EXPECT_EQ(info_result, -2); // JOB had an illegal value
    
    // Test invalid ORDSEL
    info_result = slicot_ab09ax(
        DICO, JOB, 'X', N, M, P, &nr_copy,
        A.data(), LDA, B.data(), LDB, C.data(), LDC,
        HSV.data(), T.data(), LDT, TI.data(), LDTI, TOL, &IWARN,
        0 /* column-major */
    );
    EXPECT_EQ(info_result, -3); // ORDSEL had an illegal value
    
    // Test invalid NR (when ORDSEL='F')
    nr_copy = N + 1; // NR > N is invalid
    info_result = slicot_ab09ax(
        DICO, JOB, 'F', N, M, P, &nr_copy,
        A.data(), LDA, B.data(), LDB, C.data(), LDC,
        HSV.data(), T.data(), LDT, TI.data(), LDTI, TOL, &IWARN,
        0 /* column-major */
    );
    EXPECT_EQ(info_result, -7); // NR had an illegal value
    
    // Test invalid LDA (column-major)
    info_result = slicot_ab09ax(
        DICO, JOB, ORDSEL, N, M, P, &NR,
        A.data(), 0, B.data(), LDB, C.data(), LDC,
        HSV.data(), T.data(), LDT, TI.data(), LDTI, TOL, &IWARN,
        0 /* column-major */
    );
    EXPECT_EQ(info_result, -9); // LDA had an illegal value
    
    // Test invalid LDB (column-major)
    info_result = slicot_ab09ax(
        DICO, JOB, ORDSEL, N, M, P, &NR,
        A.data(), LDA, B.data(), 0, C.data(), LDC,
        HSV.data(), T.data(), LDT, TI.data(), LDTI, TOL, &IWARN,
        0 /* column-major */
    );
    EXPECT_EQ(info_result, -11); // LDB had an illegal value
}
