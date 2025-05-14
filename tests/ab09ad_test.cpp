#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm> // For std::max, std::min
#include <iostream>  // For debugging output if needed

#include "ab09ad.h"
#include "slicot_utils.h" // For transpose functions if needed for setup/verification
// #include "test_config.h" // Assuming this defines TEST_DATA_DIR or similar if loading from CSV

// Using the example from AB09AD documentation
class AB09ADTestColMajor : public ::testing::Test {
protected:
    // Test parameters
    char DICO_in = 'C';      // Continuous-time system
    char JOB_in = 'N';       // Use balancing-free square-root method
    char EQUIL_in = 'N';     // Do not perform equilibration
    char ORDSEL_in = 'A';    // Order is automatically determined by TOL
    int N_in = 7;            // Order of the original state-space
    int M_in = 2;            // Number of inputs
    int P_in = 3;            // Number of outputs
    int NR_io = 0;           // Desired reduced order (0 for auto with ORDSEL='A')
    double TOL_in = 1.0e-1;  // Tolerance for determining the order
    
    int IWARN_out = 0;       // Warning indicator output
    int info_result = -999;  // Result from slicot_ab09ad call

    // Verification tolerance
    double check_tol = 1e-4;

    // Input data vectors (column-major format)
    std::vector<double> A_io;  // Original state matrix, becomes reduced Ar
    std::vector<double> B_io;  // Original input matrix, becomes reduced Br
    std::vector<double> C_io;  // Original output matrix, becomes reduced Cr
    std::vector<double> HSV_out; // Hankel singular values output

    // Expected results (from documentation)
    std::vector<double> HSV_expected; 
    std::vector<double> A_expected_reduced;   // Expected reduced Ar matrix (NR_expected x NR_expected)
    std::vector<double> B_expected_reduced;   // Expected reduced Br matrix (NR_expected x M_in)
    std::vector<double> C_expected_reduced;   // Expected reduced Cr matrix (P_in x NR_expected)
    int NR_expected = 5;  // Expected reduced order from example
    int info_expected = 0; // Expected info from slicot_ab09ad call

    // Leading dimensions for C arrays
    int LDA_in = 0;
    int LDB_in = 0;
    int LDC_in = 0;

    void SetUp() override {
        // Set leading dimensions for column-major format (Fortran style)
        LDA_in = std::max(1, N_in);
        LDB_in = std::max(1, N_in); // B is N x M, so LDB >= N
        LDC_in = std::max(1, P_in); // C is P x N, so LDC >= P

        // Initialize matrices based on example from documentation (column-major)
        A_io = {
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
        // Resize to LDA_in * N_in if LDA_in > N_in (though for this example N_in=7, LDA_in=7)
        A_io.resize((size_t)LDA_in * N_in);


        B_io = {
            // Column 1
            0.0, 12.5, 0.0, 0.0, 0.0, 0.0, 0.0, 
            // Column 2
            0.0, 0.0, 0.0, 0.0, 0.0, 12.5, 0.0
        };
        B_io.resize((size_t)LDB_in * M_in);


        C_io = {
            // Column 1
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            // Column 2
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
            // Column 3
            0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0
        };
        // Transpose C for storage: C is PxN. Stored column-wise means C(i,j) is C[j*LDC + i]
        // The example data is given row-wise for C. We need to store it column-wise for the test.
        std::vector<double> C_temp_row_major = {
             1.0000,   0.0000,  0.0000,   0.0000,  0.0000,  0.0000,  0.0000, // Row 1
             0.0000,   0.0000,  0.0000,   1.0000,  0.0000,  0.0000,  0.0000, // Row 2
             0.0000,   0.0000,  0.0000,   0.0000,  1.0000,  0.0000,  0.0000  // Row 3
        };
        C_io.resize((size_t)LDC_in * N_in); // LDC_in is P_in for column major C
        slicot_transpose_to_fortran_with_ld(C_temp_row_major.data(), C_io.data(), P_in, N_in, N_in, LDC_in, sizeof(double));


        HSV_out.resize(N_in);

        // Expected results for reduced system (order 5)
        HSV_expected = {2.5139, 2.0846, 1.9178, 0.7666, 0.5473, 0.0253, 0.0246};
        
        // Expected reduced Ar matrix (5x5, column-major)
        A_expected_reduced = {
            1.3451, -4.0214,  0.0000,  0.0000,  1.2402,
            5.0399, -3.6604,  0.0000,  0.0000,  1.6416,
            0.0000,  0.0000,  0.5124, -4.2167,  0.0000,
            0.0000,  0.0000,  1.7910, -2.9900,  0.0000,
            4.5315, -0.9056,  0.0000,  0.0000, -0.0586
        };

        // Expected reduced Br matrix (5x2, column-major)
        B_expected_reduced = {
            -0.3857, -3.1753, -0.7447, -3.6872,  1.8197, // Col 1 of Br
             0.3857,  3.1753, -0.7447, -3.6872, -1.8197  // Col 2 of Br
        };

        // Expected reduced Cr matrix (3x5, column-major)
        C_expected_reduced = {
            -0.6704,  0.1089,  0.6704, // Col 1 of Cr
             0.1828,  0.4867, -0.1828, // Col 2 of Cr
            -0.6582,  0.0000, -0.6582, // Col 3 of Cr
             0.2222,  0.0000,  0.2222, // Col 4 of Cr
            -0.0104,  0.8651,  0.0104  // Col 5 of Cr
        };
    }
};

// Row-major test fixture derived from the column-major fixture
class AB09ADTestRowMajor : public AB09ADTestColMajor {
protected:
    std::vector<double> A_rm_io;
    std::vector<double> B_rm_io;
    std::vector<double> C_rm_io;
    // HSV_out is 1D, same for row/col major

    // Expected results in row-major format for the reduced parts
    std::vector<double> A_expected_reduced_rm;
    std::vector<double> B_expected_reduced_rm;
    std::vector<double> C_expected_reduced_rm;

    void SetUp() override {
        AB09ADTestColMajor::SetUp(); // Initialize column-major base data

        // For row-major C arrays, leading dimensions are number of columns
        LDA_in = std::max(1, N_in); // A is N x N, LDA (cols) >= N
        LDB_in = std::max(1, M_in); // B is N x M, LDB (cols) >= M
        LDC_in = std::max(1, N_in); // C is P x N, LDC (cols) >= N

        A_rm_io.resize((size_t)N_in * LDA_in);
        B_rm_io.resize((size_t)N_in * LDB_in); // N_in rows, LDB_in columns (M_in actual cols)
        C_rm_io.resize((size_t)P_in * LDC_in); // P_in rows, LDC_in columns (N_in actual cols)

        // Transpose initial A_io (N_in x N_in, col-major) to A_rm_io (N_in x N_in, row-major)
        slicot_transpose_to_c_with_ld(A_io.data(), A_rm_io.data(), N_in, N_in, 
                                      std::max(1,N_in) /*LD src (rows)*/, LDA_in /*LD dest (cols)*/, sizeof(double));
        
        // Transpose initial B_io (N_in x M_in, col-major) to B_rm_io (N_in x M_in, row-major)
        slicot_transpose_to_c_with_ld(B_io.data(), B_rm_io.data(), N_in, M_in, 
                                      std::max(1,N_in) /*LD src (rows)*/, LDB_in /*LD dest (cols)*/, sizeof(double));

        // Transpose initial C_io (P_in x N_in, col-major) to C_rm_io (P_in x N_in, row-major)
        slicot_transpose_to_c_with_ld(C_io.data(), C_rm_io.data(), P_in, N_in, 
                                      std::max(1,P_in) /*LD src (rows)*/, LDC_in /*LD dest (cols)*/, sizeof(double));


        // Convert expected *reduced* matrices to row-major for comparison
        A_expected_reduced_rm.resize((size_t)NR_expected * NR_expected);
        slicot_transpose_to_c_with_ld(A_expected_reduced.data(), A_expected_reduced_rm.data(), 
                                      NR_expected, NR_expected, 
                                      std::max(1,NR_expected), std::max(1,NR_expected), sizeof(double));

        B_expected_reduced_rm.resize((size_t)NR_expected * M_in);
        slicot_transpose_to_c_with_ld(B_expected_reduced.data(), B_expected_reduced_rm.data(), 
                                      NR_expected, M_in, 
                                      std::max(1,NR_expected), std::max(1,M_in), sizeof(double));
        
        C_expected_reduced_rm.resize((size_t)P_in * NR_expected);
        slicot_transpose_to_c_with_ld(C_expected_reduced.data(), C_expected_reduced_rm.data(), 
                                      P_in, NR_expected, 
                                      std::max(1,P_in), std::max(1,NR_expected), sizeof(double));
    }
};

// Test: Documentation Example (Column-Major)
TEST_F(AB09ADTestColMajor, DocExample) {
    NR_io = 0; // Auto determine order
    
    info_result = slicot_ab09ad(DICO_in, JOB_in, EQUIL_in, ORDSEL_in,
                                N_in, M_in, P_in, &NR_io,
                                A_io.data(), LDA_in,
                                B_io.data(), LDB_in,
                                C_io.data(), LDC_in,
                                HSV_out.data(), TOL_in, &IWARN_out,
                                0 /* column-major */);

    ASSERT_EQ(info_result, info_expected);
    EXPECT_EQ(NR_io, NR_expected);

    for (int i = 0; i < N_in; i++) {
        EXPECT_NEAR(HSV_out[i], HSV_expected[i], check_tol) << "HSV mismatch at index " << i;
    }

    // Verify the reduced Ar (NR_expected x NR_expected) in the top-left of A_io
    for (int j = 0; j < NR_expected; ++j) { // Iterate over columns
        for (int i = 0; i < NR_expected; ++i) { // Iterate over rows
            EXPECT_NEAR(A_io[j * LDA_in + i], A_expected_reduced[j * NR_expected + i], check_tol)
                << "A_reduced (col-major) mismatch at (" << i << "," << j << ")";
        }
    }

    // Verify the reduced Br (NR_expected x M_in) in the top-left of B_io
    for (int j = 0; j < M_in; ++j) { // Iterate over columns
        for (int i = 0; i < NR_expected; ++i) { // Iterate over rows
            EXPECT_NEAR(B_io[j * LDB_in + i], B_expected_reduced[j * NR_expected + i], check_tol)
                << "B_reduced (col-major) mismatch at (" << i << "," << j << ")";
        }
    }
    
    // Verify the reduced Cr (P_in x NR_expected) in the top-left of C_io
    for (int j = 0; j < NR_expected; ++j) { // Iterate over columns
        for (int i = 0; i < P_in; ++i) { // Iterate over rows
            EXPECT_NEAR(C_io[j * LDC_in + i], C_expected_reduced[j * P_in + i], check_tol)
                << "C_reduced (col-major) mismatch at (" << i << "," << j << ")";
        }
    }
}

// Test: Documentation Example (Row-Major)
TEST_F(AB09ADTestRowMajor, DocExample) {
    NR_io = 0; // Auto determine order
    
    info_result = slicot_ab09ad(DICO_in, JOB_in, EQUIL_in, ORDSEL_in,
                                N_in, M_in, P_in, &NR_io,
                                A_rm_io.data(), LDA_in, // LDA_in is cols for row-major
                                B_rm_io.data(), LDB_in, // LDB_in is cols for row-major
                                C_rm_io.data(), LDC_in, // LDC_in is cols for row-major
                                HSV_out.data(), TOL_in, &IWARN_out,
                                1 /* row-major */);

    ASSERT_EQ(info_result, info_expected);
    EXPECT_EQ(NR_io, NR_expected);

    for (int i = 0; i < N_in; i++) {
        EXPECT_NEAR(HSV_out[i], HSV_expected[i], check_tol) << "HSV mismatch at index " << i;
    }
    
    // Verify the reduced Ar (NR_expected x NR_expected) in the top-left of A_rm_io (row-major)
    for (int i = 0; i < NR_expected; ++i) { // Iterate over rows
        for (int j = 0; j < NR_expected; ++j) { // Iterate over columns
            EXPECT_NEAR(A_rm_io[i * LDA_in + j], A_expected_reduced_rm[i * NR_expected + j], check_tol)
                << "A_reduced (row-major) mismatch at (" << i << "," << j << ")";
        }
    }

    // Verify the reduced Br (NR_expected x M_in) in the top-left of B_rm_io (row-major)
    for (int i = 0; i < NR_expected; ++i) { // Iterate over rows
        for (int j = 0; j < M_in; ++j) { // Iterate over columns
             EXPECT_NEAR(B_rm_io[i * LDB_in + j], B_expected_reduced_rm[i * M_in + j], check_tol)
                << "B_reduced (row-major) mismatch at (" << i << "," << j << ")";
        }
    }
    
    // Verify the reduced Cr (P_in x NR_expected) in the top-left of C_rm_io (row-major)
    for (int i = 0; i < P_in; ++i) { // Iterate over rows
        for (int j = 0; j < NR_expected; ++j) { // Iterate over columns
            EXPECT_NEAR(C_rm_io[i * LDC_in + j], C_expected_reduced_rm[i * NR_expected + j], check_tol)
                << "C_reduced (row-major) mismatch at (" << i << "," << j << ")";
        }
    }
}

// Test: Fixed Order Selection (Column-Major)
TEST_F(AB09ADTestColMajor, FixedOrderSelection) {
    char ordsel_fixed = 'F'; 
    NR_io = 3; // Ask for order 3
    int nr_fixed_expected = 3;
    
    // Make copies as data will be overwritten
    std::vector<double> A_copy = A_io;
    std::vector<double> B_copy = B_io;
    std::vector<double> C_copy = C_io;
    std::vector<double> HSV_fixed_out = HSV_out;

    info_result = slicot_ab09ad(DICO_in, JOB_in, EQUIL_in, ordsel_fixed,
                                N_in, M_in, P_in, &NR_io,
                                A_copy.data(), LDA_in,
                                B_copy.data(), LDB_in,
                                C_copy.data(), LDC_in,
                                HSV_fixed_out.data(), TOL_in, &IWARN_out,
                                0 /* column-major */);

    ASSERT_EQ(info_result, info_expected); // Expect success
    EXPECT_EQ(NR_io, nr_fixed_expected); // Expect order 3

    // HSV should still be the full original system's HSV
    for (int i = 0; i < N_in; i++) {
        EXPECT_NEAR(HSV_fixed_out[i], HSV_expected[i], check_tol) << "HSV mismatch at index " << i;
    }
    // Further checks on A_copy, B_copy, C_copy for order 3 could be added if expected values are known.
}

// Test: Zero Dimension Case (N=0)
TEST_F(AB09ADTestColMajor, ZeroDimension) {
    int n_zero = 0, m_zero = 0, p_zero = 0;
    int nr_zero_io = 0;
    int lda_zero = 1, ldb_zero = 1, ldc_zero = 1; // Min leading dims
    int iwarn_zero = 0;
    
    info_result = slicot_ab09ad(DICO_in, JOB_in, EQUIL_in, ORDSEL_in, // ORDSEL='A'
                               n_zero, m_zero, p_zero, &nr_zero_io,
                               nullptr, lda_zero, 
                               nullptr, ldb_zero, 
                               nullptr, ldc_zero, 
                               nullptr, TOL_in, &iwarn_zero,
                               0 /* column-major */);

    ASSERT_EQ(info_result, 0); // Should succeed for N=0
    EXPECT_EQ(nr_zero_io, 0);  // Reduced order should be 0
}

// Test: Parameter Validation
TEST_F(AB09ADTestColMajor, ParameterValidation) {
    std::vector<double> dummy_A(1, 0.0);
    std::vector<double> dummy_B(1, 0.0);
    std::vector<double> dummy_C(1, 0.0);
    std::vector<double> dummy_HSV(1, 0.0);
    int dummy_nr = 0; // Valid for N=0 or N=1, ORDSEL='A'
    int dummy_iwarn = 0;
    int N_val = 1, M_val = 1, P_val = 1; // Small valid dimensions
    int LDA_val = 1, LDB_val = 1, LDC_val = 1;

    // Test invalid DICO (arg 1)
    info_result = slicot_ab09ad('X', JOB_in, EQUIL_in, ORDSEL_in, N_val, M_val, P_val, &dummy_nr, dummy_A.data(), LDA_val, dummy_B.data(), LDB_val, dummy_C.data(), LDC_val, dummy_HSV.data(), TOL_in, &dummy_iwarn, 0);
    EXPECT_EQ(info_result, -1);
    
    // Test invalid JOB (arg 2)
    info_result = slicot_ab09ad(DICO_in, 'X', EQUIL_in, ORDSEL_in, N_val, M_val, P_val, &dummy_nr, dummy_A.data(), LDA_val, dummy_B.data(), LDB_val, dummy_C.data(), LDC_val, dummy_HSV.data(), TOL_in, &dummy_iwarn, 0);
    EXPECT_EQ(info_result, -2);
    
    // Test invalid EQUIL (arg 3)
    info_result = slicot_ab09ad(DICO_in, JOB_in, 'X', ORDSEL_in, N_val, M_val, P_val, &dummy_nr, dummy_A.data(), LDA_val, dummy_B.data(), LDB_val, dummy_C.data(), LDC_val, dummy_HSV.data(), TOL_in, &dummy_iwarn, 0);
    EXPECT_EQ(info_result, -3);

    // Test invalid ORDSEL (arg 4)
    info_result = slicot_ab09ad(DICO_in, JOB_in, EQUIL_in, 'X', N_val, M_val, P_val, &dummy_nr, dummy_A.data(), LDA_val, dummy_B.data(), LDB_val, dummy_C.data(), LDC_val, dummy_HSV.data(), TOL_in, &dummy_iwarn, 0);
    EXPECT_EQ(info_result, -4);

    // Test invalid N (arg 5)
    info_result = slicot_ab09ad(DICO_in, JOB_in, EQUIL_in, ORDSEL_in, -1, M_val, P_val, &dummy_nr, dummy_A.data(), LDA_val, dummy_B.data(), LDB_val, dummy_C.data(), LDC_val, dummy_HSV.data(), TOL_in, &dummy_iwarn, 0);
    EXPECT_EQ(info_result, -5);

    // Test invalid M (arg 6)
    info_result = slicot_ab09ad(DICO_in, JOB_in, EQUIL_in, ORDSEL_in, N_val, -1, P_val, &dummy_nr, dummy_A.data(), LDA_val, dummy_B.data(), LDB_val, dummy_C.data(), LDC_val, dummy_HSV.data(), TOL_in, &dummy_iwarn, 0);
    EXPECT_EQ(info_result, -6);

    // Test invalid P (arg 7)
    info_result = slicot_ab09ad(DICO_in, JOB_in, EQUIL_in, ORDSEL_in, N_val, M_val, -1, &dummy_nr, dummy_A.data(), LDA_val, dummy_B.data(), LDB_val, dummy_C.data(), LDC_val, dummy_HSV.data(), TOL_in, &dummy_iwarn, 0);
    EXPECT_EQ(info_result, -7);
    
    // Test invalid NR for ORDSEL='F' (arg 8)
    dummy_nr = -1;
    info_result = slicot_ab09ad(DICO_in, JOB_in, EQUIL_in, 'F', N_val, M_val, P_val, &dummy_nr, dummy_A.data(), LDA_val, dummy_B.data(), LDB_val, dummy_C.data(), LDC_val, dummy_HSV.data(), TOL_in, &dummy_iwarn, 0);
    EXPECT_EQ(info_result, -8);
    dummy_nr = N_val + 1; // NR > N
    info_result = slicot_ab09ad(DICO_in, JOB_in, EQUIL_in, 'F', N_val, M_val, P_val, &dummy_nr, dummy_A.data(), LDA_val, dummy_B.data(), LDB_val, dummy_C.data(), LDC_val, dummy_HSV.data(), TOL_in, &dummy_iwarn, 0);
    EXPECT_EQ(info_result, -8);
    dummy_nr = 0; // Reset for next tests

    // Test NULL A with N > 0 (arg 9)
    if (N_val > 0) {
        info_result = slicot_ab09ad(DICO_in, JOB_in, EQUIL_in, ORDSEL_in, N_val, M_val, P_val, &dummy_nr, nullptr, LDA_val, dummy_B.data(), LDB_val, dummy_C.data(), LDC_val, dummy_HSV.data(), TOL_in, &dummy_iwarn, 0);
        EXPECT_EQ(info_result, -9);
    }

    // Test invalid LDA (arg 10)
    if (N_val > 0) {
        info_result = slicot_ab09ad(DICO_in, JOB_in, EQUIL_in, ORDSEL_in, N_val, M_val, P_val, &dummy_nr, dummy_A.data(), 0, dummy_B.data(), LDB_val, dummy_C.data(), LDC_val, dummy_HSV.data(), TOL_in, &dummy_iwarn, 0);
        EXPECT_EQ(info_result, -10);
    }
    
    // Test NULL B with N > 0, M > 0 (arg 11)
    if (N_val > 0 && M_val > 0) {
        info_result = slicot_ab09ad(DICO_in, JOB_in, EQUIL_in, ORDSEL_in, N_val, M_val, P_val, &dummy_nr, dummy_A.data(), LDA_val, nullptr, LDB_val, dummy_C.data(), LDC_val, dummy_HSV.data(), TOL_in, &dummy_iwarn, 0);
        EXPECT_EQ(info_result, -11);
    }

    // Test invalid LDB (arg 12)
    if (N_val > 0 && M_val > 0) { // LDB must be >= N for col-major
        info_result = slicot_ab09ad(DICO_in, JOB_in, EQUIL_in, ORDSEL_in, N_val, M_val, P_val, &dummy_nr, dummy_A.data(), LDA_val, dummy_B.data(), 0, dummy_C.data(), LDC_val, dummy_HSV.data(), TOL_in, &dummy_iwarn, 0);
        EXPECT_EQ(info_result, -12);
    }
     // Test invalid LDB (row-major: LDB must be >= M)
    if (N_val > 0 && M_val > 0) {
        info_result = slicot_ab09ad(DICO_in, JOB_in, EQUIL_in, ORDSEL_in, N_val, M_val, P_val, &dummy_nr, dummy_A.data(), N_val, dummy_B.data(), 0, dummy_C.data(), N_val, dummy_HSV.data(), TOL_in, &dummy_iwarn, 1); // LDB_in = M_val for row-major
        EXPECT_EQ(info_result, -12);
    }


    // Test NULL C with P > 0, N > 0 (arg 13)
    if (P_val > 0 && N_val > 0) {
         info_result = slicot_ab09ad(DICO_in, JOB_in, EQUIL_in, ORDSEL_in, N_val, M_val, P_val, &dummy_nr, dummy_A.data(), LDA_val, dummy_B.data(), LDB_val, nullptr, LDC_val, dummy_HSV.data(), TOL_in, &dummy_iwarn, 0);
        EXPECT_EQ(info_result, -13);
    }

    // Test invalid LDC (arg 14)
    if (P_val > 0 && N_val > 0) { // LDC must be >= P for col-major
         info_result = slicot_ab09ad(DICO_in, JOB_in, EQUIL_in, ORDSEL_in, N_val, M_val, P_val, &dummy_nr, dummy_A.data(), LDA_val, dummy_B.data(), LDB_val, dummy_C.data(), 0, dummy_HSV.data(), TOL_in, &dummy_iwarn, 0);
        EXPECT_EQ(info_result, -14);
    }
    // Test invalid LDC (row-major: LDC must be >= N)
    if (P_val > 0 && N_val > 0) {
         info_result = slicot_ab09ad(DICO_in, JOB_in, EQUIL_in, ORDSEL_in, N_val, M_val, P_val, &dummy_nr, dummy_A.data(), N_val, dummy_B.data(), M_val, dummy_C.data(), 0, dummy_HSV.data(), TOL_in, &dummy_iwarn, 1); // LDC_in = N_val for row-major
        EXPECT_EQ(info_result, -14);
    }


    // Test NULL HSV with N > 0 (arg 15)
    if (N_val > 0) {
        info_result = slicot_ab09ad(DICO_in, JOB_in, EQUIL_in, ORDSEL_in, N_val, M_val, P_val, &dummy_nr, dummy_A.data(), LDA_val, dummy_B.data(), LDB_val, dummy_C.data(), LDC_val, nullptr, TOL_in, &dummy_iwarn, 0);
        EXPECT_EQ(info_result, -15);
    }
    
    // Test NULL IWARN (arg 20)
    info_result = slicot_ab09ad(DICO_in, JOB_in, EQUIL_in, ORDSEL_in, N_val, M_val, P_val, &dummy_nr, dummy_A.data(), LDA_val, dummy_B.data(), LDB_val, dummy_C.data(), LDC_val, dummy_HSV.data(), TOL_in, nullptr, 0);
    EXPECT_EQ(info_result, -20);

}
