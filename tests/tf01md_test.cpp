#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <string>
#include <stdexcept> // For std::runtime_error
#include <algorithm> // For std::max
#include <numeric>   // For std::iota
#include <iostream>  // For debugging output

#include "tf01md.h"       // Include the wrapper header
#include "slicot_utils.h" // For transpose functions if needed for setup/verification
// #include "test_utils.h"   // No longer needed for load_test_data_from_csv
#include "test_config.h"  // Still needed for TEST_DATA_DIR potentially used elsewhere

// --- Column-Major Test Fixture ---
class Tf01mdTestColMajor : public ::testing::Test {
protected:
    // Test parameters (set based on .dat/.res file and HTML doc)
    int N = 3;
    int M = 2;
    int P = 2;
    int NY = 10; // Number of samples is fixed now

    // Verification tolerance
    double check_tol = 1e-4; // Using the same tolerance as before

    // Input data vectors/matrices (column-major unless specified)
    std::vector<double> A;
    std::vector<double> B;
    std::vector<double> C;
    std::vector<double> D;
    std::vector<double> X_initial;
    // U: Input data directly embedded in Fortran column-major order
    std::vector<double> U;
    // Y_expected: Expected output based on Python simulation results, in Fortran column-major order
    std::vector<double> Y_expected;

    // Expected results (Final state X is not given in example, only Y)
    int expected_info = 0;

    // Result variables (initialized in test body)
    int info_result = -999; // Initialize to indicate not run
    std::vector<double> X_final_result; // Final state (N)
    std::vector<double> Y_result;      // Output sequence (P x NY)

    // Leading dimensions (**calculated in SetUp**)
    int LDA = 1;
    int LDB = 1;
    int LDC = 1;
    int LDD = 1;
    int LDU = 1;
    int LDY = 1;

    // SetUp method: Initialize inputs, size outputs
    void SetUp() override {
        // --- Initialize Fixed Input Matrices/Vectors (from HTML Example Data) ---
        // A (N x N) = (3 x 3) Column Major
        A = {0.0, 1.0, 0.0,
             -0.07, 0.8, 0.0,
             0.015, -0.15, 0.5};
        // B (N x M) = (3 x 2) Column Major
        B = {0.0, -1.0, 1.0,
             2.0, -0.1, 1.0};
        // C (P x N) = (2 x 3) Column Major
        C = {1.0, 0.0,
             0.0, 1.0,
             1.0, 0.5};
        // D (P x M) = (2 x 2) Column Major
        D = {0.0, 1.0,
             0.5, 0.0};
        // X_initial (N)
        X_initial = {1.0, 1.0, 1.0};

        // --- Embed Time Series Data (U and Y_expected) in Fortran Column-Major ---
        // U (M x NY) = (2 x 10) Column Major
        // Fortran layout: [u1(1), u2(1), u1(2), u2(2), ..., u1(10), u2(10)]
        U = {-0.6922, 0.2614,  // k=1
             -1.4934, -0.9160, // k=2
              0.3081, -0.6030, // k=3
             -2.7726, 1.2556,  // k=4
              2.0039, 0.2951,  // k=5
             -1.5734, -0.9942, // k=6
              1.5639, 1.8957,  // k=7
             -0.9942, 0.8988,  // k=8
              1.8957, -0.0701, // k=9
              0.8988, 0.4118}; // k=10

        // Y_expected (P x NY) = (2 x 10) Column Major
        // CORRECTED based on python-control / harold simulation results
        // Fortran layout: [y1(1), y2(1), y1(2), y2(2), ..., y1(10), y2(10)]
        Y_expected = {
            2.1307,    0.8078,    // k=1 (yout[:,0])
            0.0790,    0.85726,   // k=2 (yout[:,1])
           -4.6693862, 3.015968,  // k=3 (yout[:,2])
           -2.36879076, -2.2822018,// k=4 (yout[:,3])
            0.29215713, 3.2151578, // k=5 (yout[:,4])
            1.06532945, 1.59189837,// k=6 (yout[:,5])
           -3.18626214, 4.52668565,// k=7 (yout[:,6])
            6.40255438, -0.20432487,// k=8 (yout[:,7])
            2.96941343, 6.12860886,// k=9 (yout[:,8])
            2.22130372, 4.8333377  // k=10 (yout[:,9])
        };

        // --- Size output arrays based on parameters and NY ---
        X_final_result.resize(N);
        Y_result.resize((size_t)P * NY); // P rows, NY columns

        // **Calculate Leading Dimensions** (Column Major C / Fortran)
        // Fortran/Col-major LD is number of ROWS
        LDA = std::max(1, N);
        LDB = std::max(1, N); // LDB is rows = N
        LDC = std::max(1, P); // LDC is rows = P
        LDD = std::max(1, P); // LDD is rows = P
        LDU = std::max(1, M); // U is M x NY, LD is rows = M
        LDY = std::max(1, P); // Y is P x NY, LD is rows = P
    }
};

// --- Row-Major Test Fixture ---
class Tf01mdTestRowMajor : public Tf01mdTestColMajor {
protected:
    // Input data vectors/matrices in row-major format
    std::vector<double> A_rm;
    std::vector<double> B_rm;
    std::vector<double> C_rm;
    std::vector<double> D_rm;
    // U_rm: Fortran-ordered U data transposed to C row-major
    std::vector<double> U_rm;
    // X_initial is 1D, no _rm needed
    std::vector<double> Y_result_rm; // Output Y stored row-major

    // Expected output vectors in row-major format (for comparison)
    std::vector<double> Y_expected_rm;
    // X_final_result is 1D, no _rm needed

    void SetUp() override {
        // --- Load/Initialize Column-Major Data First using base class SetUp ---
        // This initializes A, B, C, D, X_initial, U, Y_expected (col-major with corrected Y_expected)
        Tf01mdTestColMajor::SetUp();

        // --- Convert Fixed Inputs to Row-Major ---
        A_rm.resize(A.size());
        B_rm.resize(B.size());
        C_rm.resize(C.size());
        D_rm.resize(D.size());
        if (!A.empty()) slicot_transpose_to_c(A.data(), A_rm.data(), N, N, sizeof(double));
        if (!B.empty()) slicot_transpose_to_c(B.data(), B_rm.data(), N, M, sizeof(double));
        if (!C.empty()) slicot_transpose_to_c(C.data(), C_rm.data(), P, N, sizeof(double));
        if (!D.empty()) slicot_transpose_to_c(D.data(), D_rm.data(), P, M, sizeof(double));

        // --- Convert Fortran-ordered U and expected Y to Row-Major ---
        // Base::U is already in Fortran column-major order (M rows x NY cols)
        // Base::Y_expected is already in Fortran column-major order (P rows x NY cols)
        U_rm.resize(U.size());
        Y_expected_rm.resize(Y_expected.size());
        // Transpose from Fortran column-major (MxNY) U to C row-major U_rm
        if (M > 0 && NY > 0 && !U.empty()) slicot_transpose_to_c(U.data(), U_rm.data(), M, NY, sizeof(double));
        // Transpose from Fortran column-major (PxNY) Y_expected to C row-major Y_expected_rm
        if (P > 0 && NY > 0 && !Y_expected.empty()) slicot_transpose_to_c(Y_expected.data(), Y_expected_rm.data(), P, NY, sizeof(double));

        // --- Size outputs (wrapper handles output format, store result row-major) ---
        // X_final_result is already sized in base SetUp (1D)
        Y_result_rm.resize((size_t)P * NY);

        // **Calculate Leading Dimensions AFTER NY is set** (Row Major C requires cols)
        LDA = std::max(1, N);  // Cols of A
        LDB = std::max(1, M);  // Cols of B
        LDC = std::max(1, N);  // Cols of C
        LDD = std::max(1, M);  // Cols of D
        LDU = std::max(1, NY); // Cols of U
        LDY = std::max(1, NY); // Cols of Y
    }
};

// --- Test Cases ---

// Test: Documentation Example (Column-Major)
TEST_F(Tf01mdTestColMajor, DocExample) {
    // Copy initial state for the test run
    X_final_result = X_initial;

    // Call C wrapper function (workspace handled internally)
    // Pass the directly embedded, Fortran-ordered U data
    info_result = slicot_tf01md(N, M, P, NY, // Use NY from SetUp
                                (N > 0 ? A.data() : nullptr), LDA,
                                (N > 0 && M > 0 ? B.data() : nullptr), LDB,
                                (P > 0 && N > 0 ? C.data() : nullptr), LDC,
                                (P > 0 && M > 0 ? D.data() : nullptr), LDD,
                                (M > 0 && NY > 0 ? U.data() : nullptr), LDU, // Pass Fortran-ordered U
                                (N > 0 ? X_final_result.data() : nullptr), // Pass initial/final X
                                (P > 0 && NY > 0 ? Y_result.data() : nullptr), LDY, // Pass output Y buffer
                                0 /* row_major = false */);

    // Verify return code
    ASSERT_EQ(info_result, expected_info);

    // Verify output sequence Y against expected column-major values (corrected)
    ASSERT_EQ(Y_result.size(), Y_expected.size());
    for (size_t i = 0; i < Y_expected.size(); ++i) {
        EXPECT_NEAR(Y_result[i], Y_expected[i], check_tol)
            << "Mismatch in Y (Column-Major) at flat index " << i;
    }

    // Verification of X_final_result would go here if X_final_expected was available
}

// Test: Documentation Example (Row-Major)
TEST_F(Tf01mdTestRowMajor, DocExample) {
    // Copy initial state for the test run
    X_final_result = X_initial; // X is 1D, use the same buffer

    // Call C wrapper function (workspace handled internally)
    // Pass the row-major transposed U_rm data
    info_result = slicot_tf01md(N, M, P, NY, // Use NY from SetUp
                                (N > 0 ? A_rm.data() : nullptr), LDA,
                                (N > 0 && M > 0 ? B_rm.data() : nullptr), LDB,
                                (P > 0 && N > 0 ? C_rm.data() : nullptr), LDC,
                                (P > 0 && M > 0 ? D_rm.data() : nullptr), LDD,
                                (M > 0 && NY > 0 ? U_rm.data() : nullptr), LDU, // Pass row-major U_rm
                                (N > 0 ? X_final_result.data() : nullptr), // Pass initial/final X (1D)
                                (P > 0 && NY > 0 ? Y_result_rm.data() : nullptr), LDY, // Pass row-major Y buffer
                                1 /* row_major = true */);

    // Verify return code
    ASSERT_EQ(info_result, expected_info);

    // Verify output sequence Y (now in Y_result_rm) against expected row-major values (Y_expected_rm - corrected)
    ASSERT_EQ(Y_result_rm.size(), Y_expected_rm.size());
     for (size_t i = 0; i < Y_expected_rm.size(); ++i) {
        EXPECT_NEAR(Y_result_rm[i], Y_expected_rm[i], check_tol)
            << "Mismatch in Y (Row-Major) at flat index " << i;
    }

    // Verification of X_final_result would go here if X_final_expected was available
}

// Test: Parameter Validation (specific checks for the wrapper)
TEST_F(Tf01mdTestColMajor, ParameterValidation) {
    // This test checks if the C wrapper correctly validates inputs.
    // It does NOT need valid computational data, just dummy placeholders.
    std::vector<double> dummy_A(1), dummy_B(1), dummy_C(1), dummy_D(1), dummy_U(1), dummy_X(1), dummy_Y(1);
    int valid_N = 1, valid_M = 1, valid_P = 1, valid_NY = 1;
    int valid_LDA = 1, valid_LDB = 1, valid_LDC = 1, valid_LDD = 1, valid_LDU = 1, valid_LDY = 1;

    // Test invalid N
    info_result = slicot_tf01md(-1, valid_M, valid_P, valid_NY, dummy_A.data(), valid_LDA, dummy_B.data(), valid_LDB, dummy_C.data(), valid_LDC, dummy_D.data(), valid_LDD, dummy_U.data(), valid_LDU, dummy_X.data(), dummy_Y.data(), valid_LDY, 0);
    EXPECT_EQ(info_result, -1);

    // Test invalid M
    info_result = slicot_tf01md(valid_N, -1, valid_P, valid_NY, dummy_A.data(), valid_LDA, dummy_B.data(), valid_LDB, dummy_C.data(), valid_LDC, dummy_D.data(), valid_LDD, dummy_U.data(), valid_LDU, dummy_X.data(), dummy_Y.data(), valid_LDY, 0);
    EXPECT_EQ(info_result, -2);

    // Test invalid P
    info_result = slicot_tf01md(valid_N, valid_M, -1, valid_NY, dummy_A.data(), valid_LDA, dummy_B.data(), valid_LDB, dummy_C.data(), valid_LDC, dummy_D.data(), valid_LDD, dummy_U.data(), valid_LDU, dummy_X.data(), dummy_Y.data(), valid_LDY, 0);
    EXPECT_EQ(info_result, -3);

    // Test invalid NY
    info_result = slicot_tf01md(valid_N, valid_M, valid_P, -1, dummy_A.data(), valid_LDA, dummy_B.data(), valid_LDB, dummy_C.data(), valid_LDC, dummy_D.data(), valid_LDD, dummy_U.data(), valid_LDU, dummy_X.data(), dummy_Y.data(), valid_LDY, 0);
    EXPECT_EQ(info_result, -4);

    // Test NULL pointer for required array A (assuming N > 0)
    info_result = slicot_tf01md(valid_N, valid_M, valid_P, valid_NY, nullptr, valid_LDA, dummy_B.data(), valid_LDB, dummy_C.data(), valid_LDC, dummy_D.data(), valid_LDD, dummy_U.data(), valid_LDU, dummy_X.data(), dummy_Y.data(), valid_LDY, 0);
    EXPECT_EQ(info_result, -5);

    // Test invalid LDA (Column-Major)
    info_result = slicot_tf01md(valid_N, valid_M, valid_P, valid_NY, dummy_A.data(), 0, dummy_B.data(), valid_LDB, dummy_C.data(), valid_LDC, dummy_D.data(), valid_LDD, dummy_U.data(), valid_LDU, dummy_X.data(), dummy_Y.data(), valid_LDY, 0);
    EXPECT_EQ(info_result, -6);

     // Test NULL pointer for required array B (assuming N > 0, M > 0)
    info_result = slicot_tf01md(valid_N, valid_M, valid_P, valid_NY, dummy_A.data(), valid_LDA, nullptr, valid_LDB, dummy_C.data(), valid_LDC, dummy_D.data(), valid_LDD, dummy_U.data(), valid_LDU, dummy_X.data(), dummy_Y.data(), valid_LDY, 0);
    EXPECT_EQ(info_result, -7);

    // Test invalid LDB (Column-Major)
    info_result = slicot_tf01md(valid_N, valid_M, valid_P, valid_NY, dummy_A.data(), valid_LDA, dummy_B.data(), 0, dummy_C.data(), valid_LDC, dummy_D.data(), valid_LDD, dummy_U.data(), valid_LDU, dummy_X.data(), dummy_Y.data(), valid_LDY, 0);
    EXPECT_EQ(info_result, -8);

    // Test NULL C (P > 0, N > 0)
    info_result = slicot_tf01md(valid_N, valid_M, valid_P, valid_NY, dummy_A.data(), valid_LDA, dummy_B.data(), valid_LDB, nullptr, valid_LDC, dummy_D.data(), valid_LDD, dummy_U.data(), valid_LDU, dummy_X.data(), dummy_Y.data(), valid_LDY, 0);
    EXPECT_EQ(info_result, -9);

    // Test invalid LDC (Column-Major)
    info_result = slicot_tf01md(valid_N, valid_M, valid_P, valid_NY, dummy_A.data(), valid_LDA, dummy_B.data(), valid_LDB, dummy_C.data(), 0, dummy_D.data(), valid_LDD, dummy_U.data(), valid_LDU, dummy_X.data(), dummy_Y.data(), valid_LDY, 0);
    EXPECT_EQ(info_result, -10);

    // Test NULL D (P > 0, M > 0)
    info_result = slicot_tf01md(valid_N, valid_M, valid_P, valid_NY, dummy_A.data(), valid_LDA, dummy_B.data(), valid_LDB, dummy_C.data(), valid_LDC, nullptr, valid_LDD, dummy_U.data(), valid_LDU, dummy_X.data(), dummy_Y.data(), valid_LDY, 0);
    EXPECT_EQ(info_result, -11);

    // Test invalid LDD (Column-Major)
    info_result = slicot_tf01md(valid_N, valid_M, valid_P, valid_NY, dummy_A.data(), valid_LDA, dummy_B.data(), valid_LDB, dummy_C.data(), valid_LDC, dummy_D.data(), 0, dummy_U.data(), valid_LDU, dummy_X.data(), dummy_Y.data(), valid_LDY, 0);
    EXPECT_EQ(info_result, -12);

    // Test NULL U (M > 0, NY > 0)
    info_result = slicot_tf01md(valid_N, valid_M, valid_P, valid_NY, dummy_A.data(), valid_LDA, dummy_B.data(), valid_LDB, dummy_C.data(), valid_LDC, dummy_D.data(), valid_LDD, nullptr, valid_LDU, dummy_X.data(), dummy_Y.data(), valid_LDY, 0);
    EXPECT_EQ(info_result, -13);

    // Test invalid LDU (Column-Major)
    info_result = slicot_tf01md(valid_N, valid_M, valid_P, valid_NY, dummy_A.data(), valid_LDA, dummy_B.data(), valid_LDB, dummy_C.data(), valid_LDC, dummy_D.data(), valid_LDD, dummy_U.data(), 0, dummy_X.data(), dummy_Y.data(), valid_LDY, 0);
    EXPECT_EQ(info_result, -14);

    // Test NULL X (N > 0)
    info_result = slicot_tf01md(valid_N, valid_M, valid_P, valid_NY, dummy_A.data(), valid_LDA, dummy_B.data(), valid_LDB, dummy_C.data(), valid_LDC, dummy_D.data(), valid_LDD, dummy_U.data(), valid_LDU, nullptr, dummy_Y.data(), valid_LDY, 0);
    EXPECT_EQ(info_result, -15);

    // Test NULL Y (P > 0, NY > 0)
    info_result = slicot_tf01md(valid_N, valid_M, valid_P, valid_NY, dummy_A.data(), valid_LDA, dummy_B.data(), valid_LDB, dummy_C.data(), valid_LDC, dummy_D.data(), valid_LDD, dummy_U.data(), valid_LDU, dummy_X.data(), nullptr, valid_LDY, 0);
    EXPECT_EQ(info_result, -16);

    // Test invalid LDY (Column-Major)
    info_result = slicot_tf01md(valid_N, valid_M, valid_P, valid_NY, dummy_A.data(), valid_LDA, dummy_B.data(), valid_LDB, dummy_C.data(), valid_LDC, dummy_D.data(), valid_LDD, dummy_U.data(), valid_LDU, dummy_X.data(), dummy_Y.data(), 0, 0);
    EXPECT_EQ(info_result, -17);

    // Test invalid LDA (Row-Major) with N=2, LDA=1 (< N cols)
    info_result = slicot_tf01md(2, 1, 1, 1, dummy_A.data(), 1, dummy_B.data(), 2, dummy_C.data(), 2, dummy_D.data(), 1, dummy_U.data(), 1, dummy_X.data(), dummy_Y.data(), 1, 1);
    EXPECT_EQ(info_result, -6);

    // Test invalid LDU (Row-Major) with M=1, NY=2, LDU=1 (< NY cols)
    info_result = slicot_tf01md(1, 1, 1, 2, dummy_A.data(), 1, dummy_B.data(), 1, dummy_C.data(), 1, dummy_D.data(), 1, dummy_U.data(), 1, dummy_X.data(), dummy_Y.data(), 2, 1);
    EXPECT_EQ(info_result, -14);

     // Test invalid LDY (Row-Major) with P=1, NY=2, LDY=1 (< NY cols)
    info_result = slicot_tf01md(1, 1, 1, 2, dummy_A.data(), 1, dummy_B.data(), 1, dummy_C.data(), 1, dummy_D.data(), 1, dummy_U.data(), 2, dummy_X.data(), dummy_Y.data(), 1, 1);
    EXPECT_EQ(info_result, -17);

}
