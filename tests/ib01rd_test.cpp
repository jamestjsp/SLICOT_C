#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm> // For std::max

#include "ib01rd.h"       // Include the wrapper header
#include "slicot_utils.h" // For transpose functions, MAX/MIN
// #include "test_utils.h" // Not used as data is embedded
// #include "test_config.h"  // Not used as data is embedded

// --- Base Test Fixture for IB01RD ---
class IB01RDTest : public ::testing::Test {
protected:
    // Parameters for the synthetic test case
    char JOB = 'N'; // D is not zero
    int N = 2;
    int M = 1;
    int L = 1;
    int NSMP = 3;
    double TOL = 0.0; // Use machine precision

    // Input matrices (column-major by default for _col versions)
    std::vector<double> A_col;
    std::vector<double> B_col;
    std::vector<double> C_col;
    std::vector<double> D_col;
    std::vector<double> U_input_col;
    std::vector<double> Y_input_col;

    // Expected output
    std::vector<double> X0_expected_col;
    int expected_info = 0;
    int expected_iwarn = 0;

    // Result variables
    std::vector<double> X0_result_col;
    int info_result = -999;
    int iwarn_result = -999;
    std::vector<double> dwork_info_result;


    // Leading dimensions (calculated in SetUp)
    // For column-major C arrays (passed to wrapper when row_major = 0)
    int LDA_c = 1, LDB_c = 1, LDC_c = 1, LDD_c = 1, LDU_c = 1, LDY_c = 1;
    // For row-major C arrays (passed to wrapper when row_major = 1)
    int LDA_r = 1, LDB_r = 1, LDC_r = 1, LDD_r = 1, LDU_r = 1, LDY_r = 1;


    // Verification tolerance
    double check_tol = 1e-7; // Adjusted for typical numerical precision

    void SetUpBase() {
        // Define synthetic data based on thought process
        // A = [[0.5, 0.1], [0.0, 0.8]] (Schur form)
        A_col = {0.5, 0.0, 0.1, 0.8};
        // B = [[1.0], [0.5]]
        B_col = {1.0, 0.5};
        // C = [[1.0, 0.0]]
        C_col = {1.0, 0.0}; // L x N = 1 x 2
        // D = [[0.2]]
        D_col = {0.2};
        // U = [[0.1], [0.2], [0.3]] (NSMP x M = 3 x 1)
        U_input_col = {0.1, 0.2, 0.3};
        // Y = [[1.02], [0.84], [0.825]] (NSMP x L = 3 x 1)
        Y_input_col = {1.02, 0.84, 0.825};
        // Expected X0_actual = [1.0, 2.0]
        X0_expected_col = {1.0, 2.0};

        X0_result_col.resize(N > 0 ? N : 0);
        dwork_info_result.resize(2); // For optimal LDWORK and RCOND

        // Calculate Column-Major C Leading Dimensions
        LDA_c = MAX(1, N);
        LDB_c = (N > 0 && M > 0) ? N : 1; // As per doc: LDB >= N if N>0,M>0
        LDC_c = MAX(1, L);
        LDD_c = (M > 0 && JOB == 'N') ? L : 1; // As per doc: LDD >= L if M>0,JOB='N'
        LDU_c = (M > 0) ? MAX(1, NSMP) : 1;
        LDY_c = MAX(1, NSMP);

        // Calculate Row-Major C Leading Dimensions (number of columns)
        LDA_r = MAX(1, N);    // A is N x N
        LDB_r = (N > 0 && M > 0) ? M : 1;    // B is N x M
        LDC_r = (L > 0 && N > 0) ? N : 1;    // C is L x N
        LDD_r = (L > 0 && M > 0 && JOB == 'N') ? M : 1;    // D is L x M
        LDU_r = (M > 0) ? M : 1;    // U is NSMP x M
        LDY_r = (L > 0) ? L : 1;    // Y is NSMP x L
    }
};

// --- Column-Major Test Fixture ---
class IB01RDTestColMajor : public IB01RDTest {
protected:
    void SetUp() override {
        SetUpBase();
    }
};

// --- Row-Major Test Fixture ---
class IB01RDTestRowMajor : public IB01RDTest {
protected:
    std::vector<double> A_rm, B_rm, C_rm, D_rm, U_input_rm, Y_input_rm;

    void SetUp() override {
        SetUpBase();
        // Transpose inputs for row-major tests
        if (N > 0) { A_rm.resize(A_col.size()); slicot_transpose_to_c(A_col.data(), A_rm.data(), N, N, sizeof(double));}
        if (N > 0 && M > 0) { B_rm.resize(B_col.size()); slicot_transpose_to_c(B_col.data(), B_rm.data(), N, M, sizeof(double));}
        if (L > 0 && N > 0) { C_rm.resize(C_col.size()); slicot_transpose_to_c(C_col.data(), C_rm.data(), L, N, sizeof(double));}
        if (L > 0 && M > 0 && JOB == 'N') { D_rm.resize(D_col.size()); slicot_transpose_to_c(D_col.data(), D_rm.data(), L, M, sizeof(double));}
        if (NSMP > 0 && M > 0) { U_input_rm.resize(U_input_col.size()); slicot_transpose_to_c(U_input_col.data(), U_input_rm.data(), NSMP, M, sizeof(double));}
        if (NSMP > 0 && L > 0) { Y_input_rm.resize(Y_input_col.size()); slicot_transpose_to_c(Y_input_col.data(), Y_input_rm.data(), NSMP, L, sizeof(double));}
    }
};


// --- Test Cases ---

TEST_F(IB01RDTestColMajor, SyntheticExample) {
    info_result = slicot_ib01rd(JOB, N, M, L, NSMP,
                                A_col.data(), LDA_c,
                                B_col.data(), LDB_c,
                                C_col.data(), LDC_c,
                                D_col.data(), LDD_c,
                                U_input_col.data(), LDU_c,
                                Y_input_col.data(), LDY_c,
                                X0_result_col.data(),
                                TOL, &iwarn_result, dwork_info_result.data(), 0 /* row_major */);

    ASSERT_EQ(info_result, expected_info);
    ASSERT_EQ(iwarn_result, expected_iwarn);

    if (info_result == 0) {
        for (int i = 0; i < N; ++i) {
            EXPECT_NEAR(X0_result_col[i], X0_expected_col[i], check_tol) << "Mismatch at X0[" << i << "]";
        }
        // Optionally check dwork_info_result
        // EXPECT_GT(dwork_info_result[0], 0); // Optimal LDWORK
        // EXPECT_GT(dwork_info_result[1], 0); // RCOND
    }
}

TEST_F(IB01RDTestRowMajor, SyntheticExample) {
    info_result = slicot_ib01rd(JOB, N, M, L, NSMP,
                                A_rm.empty() ? nullptr : A_rm.data(), LDA_r,
                                B_rm.empty() ? nullptr : B_rm.data(), LDB_r,
                                C_rm.empty() ? nullptr : C_rm.data(), LDC_r,
                                D_rm.empty() ? nullptr : D_rm.data(), LDD_r,
                                U_input_rm.empty() ? nullptr : U_input_rm.data(), LDU_r,
                                Y_input_rm.empty() ? nullptr : Y_input_rm.data(), LDY_r,
                                X0_result_col.data(), // X0 is vector, result always "column-major" like
                                TOL, &iwarn_result, dwork_info_result.data(), 1 /* row_major */);

    ASSERT_EQ(info_result, expected_info);
    ASSERT_EQ(iwarn_result, expected_iwarn);

    if (info_result == 0) {
        for (int i = 0; i < N; ++i) {
            EXPECT_NEAR(X0_result_col[i], X0_expected_col[i], check_tol) << "Mismatch at X0[" << i << "]";
        }
    }
}

TEST_F(IB01RDTestColMajor, ParameterValidation) {
    std::vector<double> dummy_x0(N > 0 ? N : 1);
    // Test invalid JOB
    info_result = slicot_ib01rd('X', N, M, L, NSMP, A_col.data(), LDA_c, B_col.data(), LDB_c, C_col.data(), LDC_c, D_col.data(), LDD_c, U_input_col.data(), LDU_c, Y_input_col.data(), LDY_c, dummy_x0.data(), TOL, &iwarn_result, nullptr, 0);
    EXPECT_EQ(info_result, -1);

    // Test invalid N
    info_result = slicot_ib01rd(JOB, -1, M, L, NSMP, A_col.data(), LDA_c, B_col.data(), LDB_c, C_col.data(), LDC_c, D_col.data(), LDD_c, U_input_col.data(), LDU_c, Y_input_col.data(), LDY_c, dummy_x0.data(), TOL, &iwarn_result, nullptr, 0);
    EXPECT_EQ(info_result, -2);

    // Test invalid M
    info_result = slicot_ib01rd(JOB, N, -1, L, NSMP, A_col.data(), LDA_c, B_col.data(), LDB_c, C_col.data(), LDC_c, D_col.data(), LDD_c, U_input_col.data(), LDU_c, Y_input_col.data(), LDY_c, dummy_x0.data(), TOL, &iwarn_result, nullptr, 0);
    EXPECT_EQ(info_result, -3);

    // Test invalid L
    info_result = slicot_ib01rd(JOB, N, M, 0, NSMP, A_col.data(), LDA_c, B_col.data(), LDB_c, C_col.data(), LDC_c, D_col.data(), LDD_c, U_input_col.data(), LDU_c, Y_input_col.data(), LDY_c, dummy_x0.data(), TOL, &iwarn_result, nullptr, 0);
    EXPECT_EQ(info_result, -4);

    // Test invalid NSMP
    if (N > 0) { // Only if N > 0, otherwise NSMP < N is NSMP < 0
      info_result = slicot_ib01rd(JOB, N, M, L, N - 1, A_col.data(), LDA_c, B_col.data(), LDB_c, C_col.data(), LDC_c, D_col.data(), LDD_c, U_input_col.data(), LDU_c, Y_input_col.data(), LDY_c, dummy_x0.data(), TOL, &iwarn_result, nullptr, 0);
      EXPECT_EQ(info_result, -5);
    }

    // Test NULL A if N > 0
    if (N > 0) {
        info_result = slicot_ib01rd(JOB, N, M, L, NSMP, nullptr, LDA_c, B_col.data(), LDB_c, C_col.data(), LDC_c, D_col.data(), LDD_c, U_input_col.data(), LDU_c, Y_input_col.data(), LDY_c, dummy_x0.data(), TOL, &iwarn_result, nullptr, 0);
        EXPECT_EQ(info_result, -6);
    }
    // Test invalid LDA if N > 0
    if (N > 0) {
        info_result = slicot_ib01rd(JOB, N, M, L, NSMP, A_col.data(), 0, B_col.data(), LDB_c, C_col.data(), LDC_c, D_col.data(), LDD_c, U_input_col.data(), LDU_c, Y_input_col.data(), LDY_c, dummy_x0.data(), TOL, &iwarn_result, nullptr, 0);
        EXPECT_EQ(info_result, -7);
    }
     // Test invalid TOL
    info_result = slicot_ib01rd(JOB, N, M, L, NSMP, A_col.data(), LDA_c, B_col.data(), LDB_c, C_col.data(), LDC_c, D_col.data(), LDD_c, U_input_col.data(), LDU_c, Y_input_col.data(), LDY_c, dummy_x0.data(), 1.1, &iwarn_result, nullptr, 0);
    EXPECT_EQ(info_result, -19);

}

TEST_F(IB01RDTestColMajor, ZeroN) {
    // Override fixture defaults for this N=0 test
    N = 0;
    M = 0;    // If M=0, B, D, U are not referenced or are zero-column matrices
    L = 1;    // L > 0 is required by SLICOT. Fixture default is 1.
    NSMP = 0; // NSMP >= N (0 >= 0). If NSMP=0, U and Y are zero-row matrices.
    JOB = 'N';// D is not referenced if M=0. Fixture default is 'N'.

    // LDs must be >= 1.
    LDA_c = 1; LDB_c = 1; LDC_c = 1; LDD_c = 1; LDU_c = 1; LDY_c = 1;

    double* x0_null_ptr = nullptr; // x0 is size N. If N=0, pass nullptr.

    dwork_info_result.assign(2, 0.0); // Initialize to known state
    iwarn_result = -999; // Initialize to known state, wrapper should set it to 0 if info=0

    // For N=0, M=0, NSMP=0, all matrix data pointers can be nullptr.
    // The wrapper is responsible for correctly handling these based on dimensions
    // and passing appropriate (often NULL) pointers to Fortran.
    info_result = slicot_ib01rd(JOB, N, M, L, NSMP,
                                nullptr, LDA_c, // A_ptr will be NULL in wrapper
                                nullptr, LDB_c, // B_ptr will be NULL
                                nullptr, LDC_c, // C_ptr will be NULL
                                nullptr, LDD_c, // D_ptr will be NULL (M=0 makes D not referenced for -12 error)
                                nullptr, LDU_c, // U_ptr will be NULL
                                nullptr, LDY_c, // Y_ptr will be NULL
                                x0_null_ptr,
                                TOL, &iwarn_result, dwork_info_result.data(), 0 /* row_major */);

    // The Fortran routine IB01RD should handle N=0 as a valid case,
    // likely performing no computation and returning INFO = 0.
    EXPECT_EQ(info_result, 0);
    // If INFO = 0, IWARN should be 0 unless a specific warning for N=0 is documented.
    // The wrapper sets *iwarn = local_iwarn. local_iwarn is 0 if Fortran call is successful
    // and Fortran's IWARN is 0.
    EXPECT_EQ(iwarn_result, 0);
}

TEST_F(IB01RDTestColMajor, ZeroM_JobZ) {
    M = 0; JOB = 'Z'; // D is not referenced
    // N, L, NSMP use fixture defaults (N=2, L=1, NSMP=3)
    // B, D, U can be NULL or point to zero-size data because M=0 or JOB='Z'
    LDB_c = 1; LDD_c = 1; LDU_c = 1; // Ensure LDs are valid
    X0_result_col.assign(N, 0.0); // Reinitialize for this test
    dwork_info_result.assign(2, 0.0);
    iwarn_result = -999;


    // Since M=0, B and U are effectively not used.
    // Since JOB='Z', D is not used.
    // The original Y_input_col was generated with M=1, D!=0.
    // For this test to be numerically meaningful for X0, Y_input_col would need
    // to correspond to a system with M=0.
    // Here, we primarily test that the call is valid and returns INFO=0.
    // A more rigorous test would compute an appropriate Y_input_col and X0_expected_col.

    info_result = slicot_ib01rd(JOB, N, M, L, NSMP,
                                A_col.data(), LDA_c,
                                nullptr, LDB_c, // B (M=0)
                                C_col.data(), LDC_c, // C
                                nullptr, LDD_c, // D (JOB='Z')
                                nullptr, LDU_c, // U (M=0)
                                Y_input_col.data(), LDY_c, // Y
                                X0_result_col.data(),
                                TOL, &iwarn_result, dwork_info_result.data(), 0);
    EXPECT_EQ(info_result, 0); // Expect success
    EXPECT_EQ(iwarn_result, 0); // Expect no warnings
}


TEST_F(IB01RDTestColMajor, ZeroM_JobN) {
    M = 0; JOB = 'N'; // D is specified as non-zero, but M=0 means Du(k) term is zero.
                     // D matrix itself is not referenced by Fortran if M=0.
    LDB_c = 1; LDD_c = 1; LDU_c = 1; // Ensure LDs are valid
    X0_result_col.assign(N, 0.0);
    dwork_info_result.assign(2, 0.0);
    iwarn_result = -999;

    // Similar to ZeroM_JobZ, this tests call validity.
    // Y_input_col and X0_expected_col would need to be specific to M=0.

    info_result = slicot_ib01rd(JOB, N, M, L, NSMP,
                                A_col.data(), LDA_c,
                                nullptr, LDB_c, // B (M=0)
                                C_col.data(), LDC_c, // C
                                nullptr, LDD_c, // D (not referenced by Fortran if M=0, wrapper validation passes if M=0)
                                nullptr, LDU_c, // U (M=0)
                                Y_input_col.data(), LDY_c, // Y
                                X0_result_col.data(),
                                TOL, &iwarn_result, dwork_info_result.data(), 0);
    EXPECT_EQ(info_result, 0);
    EXPECT_EQ(iwarn_result, 0);
}
