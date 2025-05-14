#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm> // For std::max, std::min
#include <stdexcept> // For std::runtime_error
#include <iostream>  // For debugging output if needed

#include "ab09ax.h"
#include "slicot_utils.h" // For transpose functions if needed for setup/verification
// #include "test_config.h" // Assuming this defines TEST_DATA_DIR or similar if loading from CSV

// --- Column-Major Test Fixture ---
class AB09AXTestColMajor : public ::testing::Test {
protected:
    // Test parameters
    char DICO = 'C';   // Continuous-time system
    char JOB = 'B';    // Balance & Truncate method
    char ORDSEL = 'F'; // Fixed order for FixedOrderReduction, 'A' for AutoOrderSelection
    int N = 4;         // Original system order
    int M = 2;         // Number of inputs
    int P = 3;         // Number of outputs
    int NR = 2;        // Desired reduced order for fixed order tests
    double TOL = 0.0;  // Tolerance (default for fixed, specific for auto)
    int IWARN_out = 0;     // Warning indicator output

    // Verification tolerance - increased slightly as example results might be rounded
    double check_tol = 1e-8; 

    // Input/output data vectors (column-major)
    std::vector<double> A_io; // State matrix (must be in Schur form)
    std::vector<double> B_io; // Input matrix
    std::vector<double> C_io; // Output matrix
    std::vector<double> HSV_out; // Hankel singular values output
    std::vector<double> T_out;   // Right truncation matrix output
    std::vector<double> TI_out;  // Left truncation matrix output

    // Expected results (placeholders, actual values depend on the specific A, B, C and method)
    // For a real test with verification, these would come from a trusted source (e.g., MATLAB, example output)
    std::vector<double> A_expected_reduced; 
    std::vector<double> B_expected_reduced; 
    std::vector<double> C_expected_reduced; 
    std::vector<double> HSV_expected;     
    int expected_info = 0;
    int expected_iwarn = 0;
    int expected_nr_fixed = 2; // For fixed order test
    // For auto order, expected_nr_auto would depend on TOL and HSV

    // Result variables
    int info_result = -999;

    // Leading dimensions for C arrays (will be set in SetUp)
    int LDA_io = 0;
    int LDB_io = 0;
    int LDC_io = 0;
    int LDT_out = 0;
    int LDTI_out = 0;

    void SetUp() override {
        // Set leading dimensions for column-major format (Fortran style)
        LDA_io = std::max(1, N);
        LDB_io = std::max(1, N); 
        LDC_io = std::max(1, P); 
        LDT_out = std::max(1, N);  // T is N x NR
        LDTI_out = std::max(1, N); // TI is NR x N (Fortran LDTI >= NR, C wrapper uses N for allocation safety)


        // Initialize A in REAL SCHUR FORM (upper quasi-triangular)
        // Using the transpose of the previous lower-triangular matrix.
        // A = [ -1.0   1.0   0.0   0.0 ]
        //     [  0.0  -2.0   0.0   0.0 ]
        //     [  0.0   0.0  -3.0   1.0 ]
        //     [  0.0   0.0   0.0  -4.0 ]
        // Eigenvalues: -1, -2, -3, -4 (stable for DICO='C')
        A_io = {
            // Column 1
            -1.0,  0.0,  0.0,  0.0,
            // Column 2
             1.0, -2.0,  0.0,  0.0,
            // Column 3
             0.0,  0.0, -3.0,  0.0,
            // Column 4
             0.0,  0.0,  1.0, -4.0
        };
        A_io.resize((size_t)LDA_io * N); // Ensure correct size if LDA_io > N

        // Input matrix B (N x M = 4 x 2)
        B_io = {
            // Column 1 of B
            1.0, 0.0, 1.0, 0.0,
            // Column 2 of B
            0.0, 1.0, 1.0, 1.0
        };
        B_io.resize((size_t)LDB_io * M);

        // Output matrix C (P x N = 3 x 4), stored column-major for Fortran
        // C = [ 1.0  0.0  1.0  0.0 ]
        //     [ 0.0  1.0  0.0  1.0 ]
        //     [ 1.0  1.0  0.0  0.0 ]
        C_io = {
            // Column 1 of C (C(:,0))
            1.0, 0.0, 1.0,
            // Column 2 of C (C(:,1))
            0.0, 1.0, 1.0,
            // Column 3 of C (C(:,2))
            1.0, 0.0, 0.0,
            // Column 4 of C (C(:,3))
            0.0, 1.0, 0.0
        };
        C_io.resize((size_t)LDC_io * N); // Ensure correct size if LDC_io > P

        // Allocate output HSV, T, TI vectors
        HSV_out.resize(N);
        T_out.resize((size_t)LDT_out * N);    // Max possible NR is N for allocation
        TI_out.resize((size_t)LDTI_out * N);  // Max possible NR is N for allocation

        // Placeholder expected HSVs (actual values depend on A, B, C)
        // For a real test, these should be pre-calculated or from a reliable source.
        HSV_expected = {3.0, 2.0, 0.5, 0.1}; // Example, replace with actual expected for new A,B,C

        // Placeholder expected reduced matrices (NR=2)
        // These would need to be re-calculated for the new A, B, C.
        // For now, just size them.
        A_expected_reduced.assign((size_t)expected_nr_fixed * expected_nr_fixed, 0.0);
        B_expected_reduced.assign((size_t)expected_nr_fixed * M, 0.0);
        C_expected_reduced.assign((size_t)P * expected_nr_fixed, 0.0);
    }
};

// --- Row-Major Test Fixture ---
class AB09AXTestRowMajor : public AB09AXTestColMajor {
protected:
    std::vector<double> A_rm_io;
    std::vector<double> B_rm_io;
    std::vector<double> C_rm_io;
    std::vector<double> T_rm_out;  // Row-major T
    std::vector<double> TI_rm_out; // Row-major TI
    
    // Expected reduced matrices in row-major
    std::vector<double> A_expected_reduced_rm;
    std::vector<double> B_expected_reduced_rm;
    std::vector<double> C_expected_reduced_rm;


    void SetUp() override {
        AB09AXTestColMajor::SetUp(); // Initialize column-major base data

        // For row-major C arrays, leading dimensions are number of columns
        LDA_io = std::max(1, N); 
        LDB_io = std::max(1, M); 
        LDC_io = std::max(1, N); 
        LDT_out = std::max(1, N); // T is N x NR, C LD is N (cols)
        LDTI_out = std::max(1, N); // TI is NR x N, C LD is N (cols)


        A_rm_io.resize((size_t)N * LDA_io);
        B_rm_io.resize((size_t)N * LDB_io); 
        C_rm_io.resize((size_t)P * LDC_io); 
        T_rm_out.resize((size_t)N * LDT_out); 
        TI_rm_out.resize((size_t)N * LDTI_out); // Max possible NR is N for TI rows in RM storage

        slicot_transpose_to_c_with_ld(A_io.data(), A_rm_io.data(), N, N, std::max(1,N), LDA_io, sizeof(double));
        slicot_transpose_to_c_with_ld(B_io.data(), B_rm_io.data(), N, M, std::max(1,N), LDB_io, sizeof(double));
        slicot_transpose_to_c_with_ld(C_io.data(), C_rm_io.data(), P, N, std::max(1,P), LDC_io, sizeof(double));
        
        // Expected reduced matrices to row-major for comparison
        A_expected_reduced_rm.assign((size_t)expected_nr_fixed * expected_nr_fixed, 0.0);
        B_expected_reduced_rm.assign((size_t)expected_nr_fixed * M, 0.0);
        C_expected_reduced_rm.assign((size_t)P * expected_nr_fixed, 0.0);
        // Transpose if A_expected_reduced has valid data
        // slicot_transpose_to_c_with_ld(A_expected_reduced.data(), A_expected_reduced_rm.data(), expected_nr_fixed, expected_nr_fixed, std::max(1,expected_nr_fixed), std::max(1,expected_nr_fixed), sizeof(double));
        // slicot_transpose_to_c_with_ld(B_expected_reduced.data(), B_expected_reduced_rm.data(), expected_nr_fixed, M, std::max(1,expected_nr_fixed), std::max(1,M), sizeof(double));
        // slicot_transpose_to_c_with_ld(C_expected_reduced.data(), C_expected_reduced_rm.data(), P, expected_nr_fixed, std::max(1,P), std::max(1,expected_nr_fixed), sizeof(double));
    }
};

// --- Test Cases ---

// Test: Fixed order reduction (Column-Major)
TEST_F(AB09AXTestColMajor, FixedOrderReduction) {
    std::vector<double> A_copy = A_io;
    std::vector<double> B_copy = B_io;
    std::vector<double> C_copy = C_io;
    int nr_run = NR; // Use the NR from fixture for fixed order
    
    info_result = slicot_ab09ax(
        DICO, JOB, 'F', N, M, P, &nr_run, // ORDSEL = 'F'
        A_copy.data(), LDA_io,
        B_copy.data(), LDB_io,
        C_copy.data(), LDC_io,
        HSV_out.data(),
        T_out.data(), LDT_out,
        TI_out.data(), LDTI_out,
        TOL, &IWARN_out,
        0 /* column-major */
    );

    ASSERT_EQ(info_result, expected_info);
    EXPECT_EQ(nr_run, expected_nr_fixed); // Check if reduced order matches expected for fixed
    // EXPECT_EQ(IWARN_out, expected_iwarn); // Add if specific iwarn is expected

    // Basic check: HSV should be non-increasing.
    if (N > 1) {
        for (int i = 0; i < N - 1; ++i) {
            EXPECT_GE(HSV_out[i], HSV_out[i+1] - 1e-9); // Allow for small numerical fuzz
        }
    }
    // NOTE: Detailed numerical verification of A_copy, B_copy, C_copy, T_out, TI_out
    // requires known correct results for the given A, B, C and NR.
    // The current A_expected_reduced etc. are placeholders.
}

// Test: Auto order selection (Column-Major)
TEST_F(AB09AXTestColMajor, AutoOrderSelection) {
    std::vector<double> A_copy = A_io;
    std::vector<double> B_copy = B_io;
    std::vector<double> C_copy = C_io;
    int nr_auto_run = 0; // NR will be determined by the routine
    double tol_auto = 0.1; // Example tolerance for auto selection
    
    info_result = slicot_ab09ax(
        DICO, JOB, 'A', N, M, P, &nr_auto_run, // ORDSEL = 'A'
        A_copy.data(), LDA_io,
        B_copy.data(), LDB_io,
        C_copy.data(), LDC_io,
        HSV_out.data(),
        T_out.data(), LDT_out,
        TI_out.data(), LDTI_out,
        tol_auto, &IWARN_out,
        0 /* column-major */
    );

    EXPECT_EQ(info_result, expected_info);
    EXPECT_GT(nr_auto_run, 0); // Expect some reduction
    EXPECT_LE(nr_auto_run, N);
    // Further checks on nr_auto_run would depend on HSV_out and tol_auto.
}

// Test: Fixed order reduction (Row-Major)
TEST_F(AB09AXTestRowMajor, FixedOrderReduction) {
    std::vector<double> A_rm_copy = A_rm_io;
    std::vector<double> B_rm_copy = B_rm_io;
    std::vector<double> C_rm_copy = C_rm_io;
    int nr_run = NR;
    
    info_result = slicot_ab09ax(
        DICO, JOB, 'F', N, M, P, &nr_run, // ORDSEL = 'F'
        A_rm_copy.data(), LDA_io, // LDA_io is cols for row-major
        B_rm_copy.data(), LDB_io, // LDB_io is cols for row-major
        C_rm_copy.data(), LDC_io, // LDC_io is cols for row-major
        HSV_out.data(),
        T_rm_out.data(), LDT_out,   // LDT_out is cols for row-major
        TI_rm_out.data(), LDTI_out, // LDTI_out is cols for row-major
        TOL, &IWARN_out,
        1 /* row-major */
    );

    ASSERT_EQ(info_result, expected_info);
    EXPECT_EQ(nr_run, expected_nr_fixed);
    // EXPECT_EQ(IWARN_out, expected_iwarn);

    if (N > 1) {
        for (int i = 0; i < N - 1; ++i) {
            EXPECT_GE(HSV_out[i], HSV_out[i+1] - 1e-9);
        }
    }
    // Numerical verification for row-major would require transposing expected results
    // or comparing element-wise carefully.
}

// Test: Zero dimensions
TEST_F(AB09AXTestColMajor, ZeroDimensions) {
    int n_zero = 0, m_zero = 0, p_zero = 0;
    int nr_zero_io = 0;
    int lda_z = 1, ldb_z = 1, ldc_z = 1, ldt_z = 1, ldti_z = 1;
    int iwarn_z = 0;
    
    info_result = slicot_ab09ax(
        DICO, JOB, 'F', n_zero, m_zero, p_zero, &nr_zero_io, // ORDSEL='F', NR=0
        nullptr, lda_z, 
        nullptr, ldb_z, 
        nullptr, ldc_z, 
        nullptr,    
        nullptr, ldt_z, 
        nullptr, ldti_z, 
        TOL, &iwarn_z,
        0 /* column-major */
    );
    
    EXPECT_EQ(info_result, 0); // Should succeed for N=0
    EXPECT_EQ(nr_zero_io, 0);  // Reduced order should be 0
}

// Test: Parameter Validation (selected checks)
TEST_F(AB09AXTestColMajor, ParameterValidation) {
    int nr_copy = NR; // Valid NR for N=4
    std::vector<double> a_val = A_io; // Use valid data for non-tested params
    std::vector<double> b_val = B_io;
    std::vector<double> c_val = C_io;
    std::vector<double> hsv_val = HSV_out;
    std::vector<double> t_val = T_out;
    std::vector<double> ti_val = TI_out;

    // Test invalid N (arg 4)
    info_result = slicot_ab09ax( DICO, JOB, ORDSEL, -1, M, P, &nr_copy, a_val.data(), LDA_io, b_val.data(), LDB_io, c_val.data(), LDC_io, hsv_val.data(), t_val.data(), LDT_out, ti_val.data(), LDTI_out, TOL, &IWARN_out, 0);
    EXPECT_EQ(info_result, -4);
    
    // Test invalid DICO (arg 1)
    info_result = slicot_ab09ax('X', JOB, ORDSEL, N, M, P, &nr_copy, a_val.data(), LDA_io, b_val.data(), LDB_io, c_val.data(), LDC_io,hsv_val.data(), t_val.data(), LDT_out, ti_val.data(), LDTI_out, TOL, &IWARN_out, 0);
    EXPECT_EQ(info_result, -1);
        
    // Test invalid NR for ORDSEL='F' (arg 7)
    int nr_invalid = N + 1; 
    info_result = slicot_ab09ax(DICO, JOB, 'F', N, M, P, &nr_invalid, a_val.data(), LDA_io, b_val.data(), LDB_io, c_val.data(), LDC_io, hsv_val.data(), t_val.data(), LDT_out, ti_val.data(), LDTI_out, TOL, &IWARN_out, 0);
    EXPECT_EQ(info_result, -7); 
    nr_invalid = -1;
    info_result = slicot_ab09ax(DICO, JOB, 'F', N, M, P, &nr_invalid, a_val.data(), LDA_io, b_val.data(), LDB_io, c_val.data(), LDC_io, hsv_val.data(), t_val.data(), LDT_out, ti_val.data(), LDTI_out, TOL, &IWARN_out, 0);
    EXPECT_EQ(info_result, -7);

    // Test invalid LDA (arg 9)
    if (N > 0) {
        info_result = slicot_ab09ax(DICO, JOB, ORDSEL, N, M, P, &nr_copy, a_val.data(), 0, b_val.data(), LDB_io, c_val.data(), LDC_io, hsv_val.data(), t_val.data(), LDT_out, ti_val.data(), LDTI_out, TOL, &IWARN_out, 0);
        EXPECT_EQ(info_result, -9);
    }
}
