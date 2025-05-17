#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm> // For std::max

#include "sb10fd.h"       // Include the wrapper header
#include "slicot_utils.h" // For transpose functions if needed for setup/verification
// test_utils.h might be needed for load_test_data_from_csv, but not for this test
// #include "test_utils.h"
// test_config.h might be needed for TEST_DATA_DIR, but not for this test
// #include "test_config.h"

// --- Test Fixture Base ---
class SB10FDTest : public ::testing::Test {
protected:
    // Test parameters (from Python script / SB10FD.html example)
    int N = 6;
    int M = 5;
    int NP = 5;
    int NCON = 2;
    int NMEAS = 2;
    double GAMMA = 15.0;
    double TOL = 1e-8;

    // Verification tolerance
    double check_tol = 1e-5; // Adjusted for typical floating point comparisons

    // Input data vectors (will be initialized in SetUp)
    std::vector<double> A_data, B_data, C_data, D_data;

    // Output data vectors
    std::vector<double> AK_out, BK_out, CK_out, DK_out;
    std::vector<double> RCOND_out;

    // Expected results (from Python script output / SB10FD.html example)
    std::vector<double> AK_expected, BK_expected, CK_expected, DK_expected;
    std::vector<double> RCOND_expected;
    int expected_info = 0;

    // Leading dimensions (will be calculated in SetUp)
    int LDA, LDB, LDC, LDD;
    int LDAK, LDBK, LDCK, LDDK;

    // SetUp method: Initialize inputs, size outputs
    void SetUpBase() {
        // Input Matrices (Row-Major by default for std::vector initialization)
        A_data = {
            -1.0,  0.0,  4.0,  5.0, -3.0, -2.0,
            -2.0,  4.0, -7.0, -2.0,  0.0,  3.0,
            -6.0,  9.0, -5.0,  0.0,  2.0, -1.0,
            -8.0,  4.0,  7.0, -1.0, -3.0,  0.0,
             2.0,  5.0,  8.0, -9.0,  1.0, -4.0,
             3.0, -5.0,  8.0,  0.0,  2.0, -6.0
        };
        B_data = {
            -3.0, -4.0, -2.0,  1.0,  0.0,
             2.0,  0.0,  1.0, -5.0,  2.0,
            -5.0, -7.0,  0.0,  7.0, -2.0,
             4.0, -6.0,  1.0,  1.0, -2.0,
            -3.0,  9.0, -8.0,  0.0,  5.0,
             1.0, -2.0,  3.0, -6.0, -2.0
        };
        C_data = {
             1.0, -1.0,  2.0, -4.0,  0.0, -3.0,
            -3.0,  0.0,  5.0, -1.0,  1.0,  1.0,
            -7.0,  5.0,  0.0, -8.0,  2.0, -2.0,
             9.0, -3.0,  4.0,  0.0,  3.0,  7.0,
             0.0,  1.0, -2.0,  1.0, -6.0, -2.0
        };
        D_data = {
             1.0, -2.0, -3.0,  0.0,  0.0,
             0.0,  4.0,  0.0,  1.0,  0.0,
             5.0, -3.0, -4.0,  0.0,  1.0,
             0.0,  1.0,  0.0,  1.0, -3.0,
             0.0,  0.0,  1.0,  7.0,  1.0
        };

        // Expected Output Matrices (Row-Major for direct comparison after wrapper call)
        AK_expected = {
            -2.8042858,   14.73666134,  4.66577725,  8.15956427,  0.08483134,  2.52900744,
             4.66088313,   3.2756311,  -3.57536944, -2.89406341,  0.23932927,  8.29203923,
           -15.31267921,  23.55924187, -7.12293316,  2.75993492,  5.97754202, -2.02848042,
           -22.06914509,  16.47575801, 12.55229806,-16.3601873,   4.4300018,  -3.31678073,
            30.67892741,  -3.90264734, -1.3868182,   26.23574568, -8.82674599, 10.48603009,
            -5.74292649,   0.05765125, 10.82158014,-11.22754753,  1.50735199,-10.72440147
        };
        BK_expected = {
            -0.15814028, -0.07925238,
            -0.92372355, -0.57178942,
             0.79837902,  0.66270873,
             0.11445791,  0.14964322,
            -0.67426686, -0.23764975,
             0.01961329, -0.75984082
        };
        CK_expected = {
            -0.24795063, -0.17131818, -0.08799277,  0.15342581,  0.50164174, -0.07299765,
             2.88097336, -0.3658349,   1.30071651,  0.39448353,  1.2244281,   2.56900597
        };
        DK_expected = {
             0.05538274,  0.13335773,
            -0.31948765,  0.03333943
        };
        RCOND_expected = {1.00000000e+00, 1.00000000e+00, 1.12409352e-02, 8.04921216e-04};

        // Size output arrays
        AK_out.resize(N * N);
        BK_out.resize(N * NMEAS);
        CK_out.resize(NCON * N);
        DK_out.resize(NCON * NMEAS);
        RCOND_out.resize(4);
    }
};

// --- Column-Major Test Fixture ---
class SB10FDTestColMajor : public SB10FDTest {
protected:
    std::vector<double> A_cm, B_cm, C_cm, D_cm; // Column-major copies for input

    void SetUp() override {
        SetUpBase();
        // Calculate Leading Dimensions (Column Major Fortran style)
        LDA = std::max(1, N);
        LDB = std::max(1, N);
        LDC = std::max(1, NP);
        LDD = std::max(1, NP);
        LDAK = std::max(1, N);
        LDBK = std::max(1, N);
        LDCK = std::max(1, NCON);
        LDDK = std::max(1, NCON);

        // Convert row-major A_data, B_data, etc. to column-major for the test call
        if (N > 0) {
            A_cm.resize(N * N);
            slicot_transpose_to_fortran_with_ld(A_data.data(), A_cm.data(), N, N, N, LDA, sizeof(double));
        }
        if (N > 0 && M > 0) {
            B_cm.resize(N * M);
            slicot_transpose_to_fortran_with_ld(B_data.data(), B_cm.data(), N, M, M, LDB, sizeof(double));
        }
        if (NP > 0 && N > 0) {
            C_cm.resize(NP * N);
            slicot_transpose_to_fortran_with_ld(C_data.data(), C_cm.data(), NP, N, N, LDC, sizeof(double));
        }
        if (NP > 0 && M > 0) {
            D_cm.resize(NP * M);
            slicot_transpose_to_fortran_with_ld(D_data.data(), D_cm.data(), NP, M, M, LDD, sizeof(double));
        }
    }
};

// --- Row-Major Test Fixture ---
class SB10FDTestRowMajor : public SB10FDTest {
protected:
    void SetUp() override {
        SetUpBase();
        // Calculate Leading Dimensions (Row Major C style - number of columns)
        LDA = std::max(1, N);    // N cols for A (N x N)
        LDB = std::max(1, M);    // M cols for B (N x M)
        LDC = std::max(1, N);    // N cols for C (NP x N)
        LDD = std::max(1, M);    // M cols for D (NP x M)
        LDAK = std::max(1, N);   // N cols for AK (N x N)
        LDBK = std::max(1, NMEAS); // NMEAS cols for BK (N x NMEAS)
        LDCK = std::max(1, N);   // N cols for CK (NCON x N)
        LDDK = std::max(1, NMEAS); // NMEAS cols for DK (NCON x NMEAS)
        // Input data (A_data, B_data, etc.) is already in row-major format
    }
};


// --- Test Cases ---

TEST_F(SB10FDTestColMajor, DocExampleColMajor) {
    int info_result = slicot_sb10fd(N, M, NP, NCON, NMEAS, GAMMA,
                                   A_cm.empty() ? nullptr : A_cm.data(), LDA,
                                   B_cm.empty() ? nullptr : B_cm.data(), LDB,
                                   C_cm.empty() ? nullptr : C_cm.data(), LDC,
                                   D_cm.empty() ? nullptr : D_cm.data(), LDD,
                                   AK_out.data(), LDAK,
                                   BK_out.data(), LDBK,
                                   CK_out.data(), LDCK,
                                   DK_out.data(), LDDK,
                                   RCOND_out.data(), TOL,
                                   0 /* row_major = false */);

    ASSERT_EQ(info_result, expected_info);

    // Convert expected results to column-major for comparison as wrapper output is col-major here
    std::vector<double> AK_expected_cm(N*N), BK_expected_cm(N*NMEAS);
    std::vector<double> CK_expected_cm(NCON*N), DK_expected_cm(NCON*NMEAS);

    if (N > 0) slicot_transpose_to_fortran_with_ld(AK_expected.data(), AK_expected_cm.data(), N, N, N, LDAK, sizeof(double));
    if (N > 0 && NMEAS > 0) slicot_transpose_to_fortran_with_ld(BK_expected.data(), BK_expected_cm.data(), N, NMEAS, NMEAS, LDBK, sizeof(double));
    if (NCON > 0 && N > 0) slicot_transpose_to_fortran_with_ld(CK_expected.data(), CK_expected_cm.data(), NCON, N, N, LDCK, sizeof(double));
    if (NCON > 0 && NMEAS > 0) slicot_transpose_to_fortran_with_ld(DK_expected.data(), DK_expected_cm.data(), NCON, NMEAS, NMEAS, LDDK, sizeof(double));

    for (size_t i = 0; i < AK_expected_cm.size(); ++i) EXPECT_NEAR(AK_out[i], AK_expected_cm[i], check_tol) << "AK mismatch at index " << i;
    for (size_t i = 0; i < BK_expected_cm.size(); ++i) EXPECT_NEAR(BK_out[i], BK_expected_cm[i], check_tol) << "BK mismatch at index " << i;
    for (size_t i = 0; i < CK_expected_cm.size(); ++i) EXPECT_NEAR(CK_out[i], CK_expected_cm[i], check_tol) << "CK mismatch at index " << i;
    for (size_t i = 0; i < DK_expected_cm.size(); ++i) EXPECT_NEAR(DK_out[i], DK_expected_cm[i], check_tol) << "DK mismatch at index " << i;
    for (size_t i = 0; i < RCOND_expected.size(); ++i) EXPECT_NEAR(RCOND_out[i], RCOND_expected[i], check_tol) << "RCOND mismatch at index " << i;
}

TEST_F(SB10FDTestRowMajor, DocExampleRowMajor) {
    int info_result = slicot_sb10fd(N, M, NP, NCON, NMEAS, GAMMA,
                                   A_data.empty() ? nullptr : A_data.data(), LDA,
                                   B_data.empty() ? nullptr : B_data.data(), LDB,
                                   C_data.empty() ? nullptr : C_data.data(), LDC,
                                   D_data.empty() ? nullptr : D_data.data(), LDD,
                                   AK_out.data(), LDAK,
                                   BK_out.data(), LDBK,
                                   CK_out.data(), LDCK,
                                   DK_out.data(), LDDK,
                                   RCOND_out.data(), TOL,
                                   1 /* row_major = true */);

    ASSERT_EQ(info_result, expected_info);

    // Output AK_out, BK_out etc. are now in row-major, compare directly with _expected
    for (size_t i = 0; i < AK_expected.size(); ++i) EXPECT_NEAR(AK_out[i], AK_expected[i], check_tol) << "AK mismatch at index " << i;
    for (size_t i = 0; i < BK_expected.size(); ++i) EXPECT_NEAR(BK_out[i], BK_expected[i], check_tol) << "BK mismatch at index " << i;
    for (size_t i = 0; i < CK_expected.size(); ++i) EXPECT_NEAR(CK_out[i], CK_expected[i], check_tol) << "CK mismatch at index " << i;
    for (size_t i = 0; i < DK_expected.size(); ++i) EXPECT_NEAR(DK_out[i], DK_expected[i], check_tol) << "DK mismatch at index " << i;
    for (size_t i = 0; i < RCOND_expected.size(); ++i) EXPECT_NEAR(RCOND_out[i], RCOND_expected[i], check_tol) << "RCOND mismatch at index " << i;
}

TEST_F(SB10FDTestColMajor, ParameterValidation) {
    // Test invalid N
    int info_result = slicot_sb10fd(-1, M, NP, NCON, NMEAS, GAMMA, nullptr, 1, nullptr, 1, nullptr, 1, nullptr, 1, nullptr, 1, nullptr, 1, nullptr, 1, nullptr, 1, RCOND_out.data(), TOL, 0);
    EXPECT_EQ(info_result, -1); // Invalid N

    // Test invalid NCON
    info_result = slicot_sb10fd(N, M, NP, M + 1, NMEAS, GAMMA, A_cm.data(), LDA, B_cm.data(), LDB, C_cm.data(), LDC, D_cm.data(), LDD, AK_out.data(), LDAK, BK_out.data(), LDBK, CK_out.data(), LDCK, DK_out.data(), LDDK, RCOND_out.data(), TOL, 0);
    EXPECT_EQ(info_result, -4); // Invalid NCON

    // Test invalid LDA (col-major)
    if (N > 0) {
        info_result = slicot_sb10fd(N, M, NP, NCON, NMEAS, GAMMA, A_cm.data(), 0, B_cm.data(), LDB, C_cm.data(), LDC, D_cm.data(), LDD, AK_out.data(), LDAK, BK_out.data(), LDBK, CK_out.data(), LDCK, DK_out.data(), LDDK, RCOND_out.data(), TOL, 0);
        EXPECT_EQ(info_result, -8); // Invalid LDA
    }
     // Test invalid NP - NMEAS < NCON
    info_result = slicot_sb10fd(N, M, NP, NCON, NP, GAMMA, A_cm.data(), LDA, B_cm.data(), LDB, C_cm.data(), LDC, D_cm.data(), LDD, AK_out.data(), LDAK, BK_out.data(), LDBK, CK_out.data(), LDCK, DK_out.data(), LDDK, RCOND_out.data(), TOL, 0);
    if (NP - NP < NCON) { // if 0 < NCON
         EXPECT_EQ(info_result, -4); // Or -5, depending on interpretation. Wrapper uses -4 for NCON related.
    }


    // Test NULL pointer for required array A (N > 0)
    if (N > 0) {
        info_result = slicot_sb10fd(N, M, NP, NCON, NMEAS, GAMMA, nullptr, LDA, B_cm.data(), LDB, C_cm.data(), LDC, D_cm.data(), LDD, AK_out.data(), LDAK, BK_out.data(), LDBK, CK_out.data(), LDCK, DK_out.data(), LDDK, RCOND_out.data(), TOL, 0);
        EXPECT_EQ(info_result, -7); // A is NULL
    }
}

TEST_F(SB10FDTestColMajor, ZeroDimensionN) {
    int n_zero = 0;
    int m_zd = 2, np_zd = 2, ncon_zd = 1, nmeas_zd = 1; // M, NP, NCON, NMEAS must be valid
    // Ensure constraints: np_zd - nmeas_zd >= ncon_zd (2-1 >= 1 -> 1 >= 1 TRUE)
    //                     m_zd - ncon_zd >= nmeas_zd (2-1 >= 1 -> 1 >= 1 TRUE)

    std::vector<double> ak_zd, bk_zd, ck_zd, dk_zd;
    std::vector<double> rcond_zd(4);

    // Sizes for N=0 outputs
    // AK: 0x0, BK: 0xNMEAS, CK: NCONx0, DK: NCONxNMEAS
    if (ncon_zd > 0 && nmeas_zd > 0) dk_zd.resize(ncon_zd * nmeas_zd);


    int ld_one = 1; // Min leading dim for 0-size arrays
    int ld_ncon = std::max(1, ncon_zd);
    int ld_nmeas = std::max(1, nmeas_zd);


    int info_result = slicot_sb10fd(n_zero, m_zd, np_zd, ncon_zd, nmeas_zd, GAMMA,
                                   nullptr, ld_one, /* A */
                                   nullptr, ld_one, /* B */
                                   nullptr, std::max(1,np_zd), /* C */
                                   nullptr, std::max(1,np_zd), /* D */
                                   nullptr, ld_one, /* AK */
                                   nullptr, ld_one, /* BK */
                                   nullptr, ld_ncon, /* CK */
                                   dk_zd.empty() ? nullptr : dk_zd.data(), ld_ncon, /* DK */
                                   rcond_zd.data(), TOL,
                                   0 /* row_major = false */);

    // For N=0, the routine might still do some calculations for Dk if NCON, NMEAS > 0.
    // Depending on the Fortran code's behavior for N=0, info could be 0 or an error
    // if it expects N > 0 for certain paths.
    // The docs say N >= 0. The example has N=6.
    // If N=0, X and Y Riccati equations are trivial.
    // Let's expect info = 0 if N=0 is handled gracefully.
    // The problem is that the controller state dimension is N. If N=0, AK, BK, CK are empty or trivial.
    // DK = D22 (part of original D related to C2, B2) in some Hinf schemes.
    // The specific SB10FD formulas for N=0 would need checking in the original paper/Fortran.
    // Assuming it should run without error for N=0.
    EXPECT_EQ(info_result, 0);

    // If NCON=0 or NMEAS=0, DK is also empty.
    // If NCON > 0 and NMEAS > 0, DK (NCON x NMEAS) might be computed.
    // For the given test values (ncon_zd=1, nmeas_zd=1), DK is 1x1.
    // We don't have expected values for DK for N=0 from the example.
    // This test mainly checks if it runs without crashing for N=0.
}

