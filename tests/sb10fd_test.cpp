#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm> // For std::max
#include <iostream> // For debugging output

#include "sb10fd.h"       // Include the wrapper header
#include "slicot_utils.h" // For transpose functions if needed for setup/verification

// --- Test Fixture Base ---
class SB10FDTest : public ::testing::Test {
protected:
    // Test parameters (from Python script / SB10FD.html example)
    int N_param = 6; // Renamed to avoid conflict with potential N macro
    int M_param = 5;
    int NP_param = 5;
    int NCON_param = 2;
    int NMEAS_param = 2;
    double GAMMA_param = 15.0;
    double TOL_param = 1e-8;

    // Verification tolerance
    double check_tol = 1e-5; 

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

    void SetUpBase() {
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

        // Resize output arrays, ensuring a minimum size of 1 for N=0 cases to avoid .data() on empty vector
        AK_out.resize(std::max(1, N_param * N_param));
        BK_out.resize(std::max(1, N_param * NMEAS_param));
        CK_out.resize(std::max(1, NCON_param * N_param));
        DK_out.resize(std::max(1, NCON_param * NMEAS_param));
        RCOND_out.resize(4); // RCOND is always size 4
    }
};

class SB10FDTestColMajor : public SB10FDTest {
protected:
    std::vector<double> A_cm, B_cm, C_cm, D_cm; 

    void SetUp() override {
        SetUpBase();
        // Fortran-style Leading Dimensions (number of rows)
        LDA = std::max(1, N_param);      // A is N x N
        LDB = std::max(1, N_param);      // B is N x M
        LDC = std::max(1, NP_param);     // C is NP x N
        LDD = std::max(1, NP_param);     // D is NP x M
        LDAK = std::max(1, N_param);     // AK is N x N
        LDBK = std::max(1, N_param);     // BK is N x NMEAS
        LDCK = std::max(1, NCON_param);  // CK is NCON x N
        LDDK = std::max(1, NCON_param);  // DK is NCON x NMEAS

        // Transpose row-major data to column-major for Fortran call
        if (N_param > 0 && N_param > 0) { // Check both dimensions for A (NxN)
            A_cm.resize(N_param * N_param);
            // For A (NxN), C row-major LDA is N_param cols
            slicot_transpose_to_fortran_with_ld(A_data.data(), A_cm.data(), N_param, N_param, N_param, LDA, sizeof(double));
        }
        if (N_param > 0 && M_param > 0) { // Check both dimensions for B (NxM)
            B_cm.resize(N_param * M_param);
            // For B (NxM), C row-major LDB is M_param cols
            slicot_transpose_to_fortran_with_ld(B_data.data(), B_cm.data(), N_param, M_param, M_param, LDB, sizeof(double));
        }
        if (NP_param > 0 && N_param > 0) { // Check both dimensions for C (NPxN)
            C_cm.resize(NP_param * N_param);
            // For C (NPxN), C row-major LDC is N_param cols
            slicot_transpose_to_fortran_with_ld(C_data.data(), C_cm.data(), NP_param, N_param, N_param, LDC, sizeof(double));
        }
        if (NP_param > 0 && M_param > 0) { // Check both dimensions for D (NPxM)
            D_cm.resize(NP_param * M_param);
            // For D (NPxM), C row-major LDD is M_param cols
            slicot_transpose_to_fortran_with_ld(D_data.data(), D_cm.data(), NP_param, M_param, M_param, LDD, sizeof(double));
        }
    }
};

class SB10FDTestRowMajor : public SB10FDTest {
protected:
    void SetUp() override {
        SetUpBase();
        // C Row-Major Leading Dimensions (number of columns)
        LDA = std::max(1, N_param);       // A (NxN) has N_param columns
        LDB = std::max(1, M_param);       // B (NxM) has M_param columns
        LDC = std::max(1, N_param);       // C (NPxN) has N_param columns
        LDD = std::max(1, M_param);       // D (NPxM) has M_param columns
        LDAK = std::max(1, N_param);      // AK (NxN) has N_param columns
        LDBK = std::max(1, NMEAS_param);  // BK (NxNMEAS) has NMEAS_param columns
        LDCK = std::max(1, N_param);      // CK (NCONxN) has N_param columns
        LDDK = std::max(1, NMEAS_param);  // DK (NCONxNMEAS) has NMEAS_param columns
    }
};


TEST_F(SB10FDTestColMajor, DocExampleColMajor) {
    int info_result = slicot_sb10fd(N_param, M_param, NP_param, NCON_param, NMEAS_param, GAMMA_param,
                                   A_cm.empty() ? nullptr : A_cm.data(), LDA,
                                   B_cm.empty() ? nullptr : B_cm.data(), LDB,
                                   C_cm.empty() ? nullptr : C_cm.data(), LDC,
                                   D_cm.empty() ? nullptr : D_cm.data(), LDD,
                                   AK_out.data(), LDAK,
                                   BK_out.data(), LDBK,
                                   CK_out.data(), LDCK,
                                   DK_out.data(), LDDK,
                                   RCOND_out.data(), TOL_param,
                                   0 /* row_major = false */);

    ASSERT_EQ(info_result, expected_info);

    // Expected data is row-major, so convert to column-major for comparison
    // AK_out, etc from wrapper (when row_major=0) are already column-major.
    std::vector<double> AK_expected_cm, BK_expected_cm, CK_expected_cm, DK_expected_cm;

    if (N_param > 0 && N_param > 0) {
        AK_expected_cm.resize(N_param * N_param);
        slicot_transpose_to_fortran_with_ld(AK_expected.data(), AK_expected_cm.data(), N_param, N_param, N_param, LDAK, sizeof(double));
        for (size_t i = 0; i < AK_expected_cm.size(); ++i) EXPECT_NEAR(AK_out[i], AK_expected_cm[i], check_tol) << "AK mismatch at index " << i;
    }
    if (N_param > 0 && NMEAS_param > 0) {
        BK_expected_cm.resize(N_param * NMEAS_param);
        slicot_transpose_to_fortran_with_ld(BK_expected.data(), BK_expected_cm.data(), N_param, NMEAS_param, NMEAS_param, LDBK, sizeof(double));
        for (size_t i = 0; i < BK_expected_cm.size(); ++i) EXPECT_NEAR(BK_out[i], BK_expected_cm[i], check_tol) << "BK mismatch at index " << i;
    }
    if (NCON_param > 0 && N_param > 0) {
        CK_expected_cm.resize(NCON_param * N_param);
        slicot_transpose_to_fortran_with_ld(CK_expected.data(), CK_expected_cm.data(), NCON_param, N_param, N_param, LDCK, sizeof(double));
        for (size_t i = 0; i < CK_expected_cm.size(); ++i) EXPECT_NEAR(CK_out[i], CK_expected_cm[i], check_tol) << "CK mismatch at index " << i;
    }
    if (NCON_param > 0 && NMEAS_param > 0) {
        DK_expected_cm.resize(NCON_param * NMEAS_param);
        slicot_transpose_to_fortran_with_ld(DK_expected.data(), DK_expected_cm.data(), NCON_param, NMEAS_param, NMEAS_param, LDDK, sizeof(double));
        for (size_t i = 0; i < DK_expected_cm.size(); ++i) EXPECT_NEAR(DK_out[i], DK_expected_cm[i], check_tol) << "DK mismatch at index " << i;
    }
    for (size_t i = 0; i < RCOND_expected.size(); ++i) EXPECT_NEAR(RCOND_out[i], RCOND_expected[i], check_tol) << "RCOND mismatch at index " << i;
}

TEST_F(SB10FDTestRowMajor, DocExampleRowMajor) {
    int info_result = slicot_sb10fd(N_param, M_param, NP_param, NCON_param, NMEAS_param, GAMMA_param,
                                   A_data.empty() ? nullptr : A_data.data(), LDA,
                                   B_data.empty() ? nullptr : B_data.data(), LDB,
                                   C_data.empty() ? nullptr : C_data.data(), LDC,
                                   D_data.empty() ? nullptr : D_data.data(), LDD,
                                   AK_out.data(), LDAK,
                                   BK_out.data(), LDBK,
                                   CK_out.data(), LDCK,
                                   DK_out.data(), LDDK,
                                   RCOND_out.data(), TOL_param,
                                   1 /* row_major = true */);

    ASSERT_EQ(info_result, expected_info);

    // Output AK_out etc. are now in row-major, compare directly with _expected (which is also row-major)
    if (N_param > 0 && N_param > 0) for (size_t i = 0; i < AK_expected.size(); ++i) EXPECT_NEAR(AK_out[i], AK_expected[i], check_tol) << "AK mismatch at index " << i;
    if (N_param > 0 && NMEAS_param > 0) for (size_t i = 0; i < BK_expected.size(); ++i) EXPECT_NEAR(BK_out[i], BK_expected[i], check_tol) << "BK mismatch at index " << i;
    if (NCON_param > 0 && N_param > 0) for (size_t i = 0; i < CK_expected.size(); ++i) EXPECT_NEAR(CK_out[i], CK_expected[i], check_tol) << "CK mismatch at index " << i;
    if (NCON_param > 0 && NMEAS_param > 0) for (size_t i = 0; i < DK_expected.size(); ++i) EXPECT_NEAR(DK_out[i], DK_expected[i], check_tol) << "DK mismatch at index " << i;
    for (size_t i = 0; i < RCOND_expected.size(); ++i) EXPECT_NEAR(RCOND_out[i], RCOND_expected[i], check_tol) << "RCOND mismatch at index " << i;
}

TEST_F(SB10FDTestColMajor, ParameterValidation) {
    std::vector<double> dummy_rcond(4); // RCOND is always size 4
    int info_result = slicot_sb10fd(-1, M_param, NP_param, NCON_param, NMEAS_param, GAMMA_param, nullptr, 1, nullptr, 1, nullptr, 1, nullptr, 1, nullptr, 1, nullptr, 1, nullptr, 1, nullptr, 1, dummy_rcond.data(), TOL_param, 0);
    EXPECT_EQ(info_result, -1);

    // Test NCON > M
    info_result = slicot_sb10fd(N_param, M_param, NP_param, M_param + 1, NMEAS_param, GAMMA_param, 
                                A_cm.empty() ? nullptr : A_cm.data(), LDA, B_cm.empty() ? nullptr : B_cm.data(), LDB, 
                                C_cm.empty() ? nullptr : C_cm.data(), LDC, D_cm.empty() ? nullptr : D_cm.data(), LDD, 
                                AK_out.data(), LDAK, BK_out.data(), LDBK, CK_out.data(), LDCK, DK_out.data(), LDDK, 
                                dummy_rcond.data(), TOL_param, 0);
    EXPECT_EQ(info_result, -4); 

    if (N_param > 0) { // LDA check only if N > 0
        info_result = slicot_sb10fd(N_param, M_param, NP_param, NCON_param, NMEAS_param, GAMMA_param, 
                                    A_cm.empty() ? nullptr : A_cm.data(), 0, B_cm.empty() ? nullptr : B_cm.data(), LDB, 
                                    C_cm.empty() ? nullptr : C_cm.data(), LDC, D_cm.empty() ? nullptr : D_cm.data(), LDD, 
                                    AK_out.data(), LDAK, BK_out.data(), LDBK, CK_out.data(), LDCK, DK_out.data(), LDDK, 
                                    dummy_rcond.data(), TOL_param, 0);
        EXPECT_EQ(info_result, -8); 
    }
    
    // Test NP - NMEAS < NCON. Let NCON_test = NP - NMEAS + 1.
    int ncon_test_invalid = NP_param - NMEAS_param + 1;
    if (ncon_test_invalid <= M_param && ncon_test_invalid >=0) { // Ensure NCON_test is still valid w.r.t M and non-negative
       info_result = slicot_sb10fd(N_param, M_param, NP_param, ncon_test_invalid, NMEAS_param, GAMMA_param, 
                                   A_cm.empty() ? nullptr : A_cm.data(), LDA, B_cm.empty() ? nullptr : B_cm.data(), LDB, 
                                   C_cm.empty() ? nullptr : C_cm.data(), LDC, D_cm.empty() ? nullptr : D_cm.data(), LDD, 
                                   AK_out.data(), LDAK, BK_out.data(), LDBK, CK_out.data(), LDCK, DK_out.data(), LDDK, 
                                   dummy_rcond.data(), TOL_param, 0);
       EXPECT_EQ(info_result, -4); // This specific constraint violation maps to -4 in the wrapper
    }

    if (N_param > 0) { // A is required if N > 0
        info_result = slicot_sb10fd(N_param, M_param, NP_param, NCON_param, NMEAS_param, GAMMA_param, 
                                    nullptr, LDA, B_cm.empty() ? nullptr : B_cm.data(), LDB, 
                                    C_cm.empty() ? nullptr : C_cm.data(), LDC, D_cm.empty() ? nullptr : D_cm.data(), LDD, 
                                    AK_out.data(), LDAK, BK_out.data(), LDBK, CK_out.data(), LDCK, DK_out.data(), LDDK, 
                                    dummy_rcond.data(), TOL_param, 0);
        EXPECT_EQ(info_result, -7); 
    }
}

TEST_F(SB10FDTestColMajor, ZeroDimensionN) {
    int n_zero = 0;
    int m_zd = 2, np_zd = 2, ncon_zd = 1, nmeas_zd = 1;
    double gamma_zd = 1.0, tol_zd = 1e-8;

    // D is NPxM = 2x2. DK is NCONxNMEAS = 1x1.
    // Initialize D with non-zero values to satisfy D12, D21 rank conditions.
    // D = [D(0,0) D(0,1); D(1,0) D(1,1)] (matlab notation)
    // Stored column-major: d_zd_data = {D(0,0), D(1,0), D(0,1), D(1,1)}
    // D11 = D(0,0), D12 = D(0,1), D21 = D(1,0), D22 = D(1,1)
    // To make D12 (D(0,1)) and D21 (D(1,0)) non-zero:
    std::vector<double> d_zd_data = {1.0, 1.0, 1.0, 1.0}; // D(0,0)=1, D(1,0)=1 (D21), D(0,1)=1 (D12), D(1,1)=1
    
    std::vector<double> dk_zd_out(std::max(1, ncon_zd * nmeas_zd)); 
    std::vector<double> rcond_zd_out(4);

    int lda_zd = 1; 
    int ldb_zd = 1; 
    int ldc_zd = std::max(1, np_zd); 
    int ldd_zd = std::max(1, np_zd); // For D (np_zd x m_zd), Fortran LDD = np_zd
    int ldak_zd = 1;
    int ldbk_zd = 1;
    int ldck_zd = std::max(1, ncon_zd);
    int lddk_zd = std::max(1, ncon_zd); // For DK (ncon_zd x nmeas_zd), Fortran LDDK = ncon_zd


    int info_result = slicot_sb10fd(n_zero, m_zd, np_zd, ncon_zd, nmeas_zd, gamma_zd,
                                   nullptr, lda_zd,       
                                   nullptr, ldb_zd,       
                                   nullptr, ldc_zd,       
                                   d_zd_data.data(), ldd_zd, 
                                   nullptr, ldak_zd,      
                                   nullptr, ldbk_zd,      
                                   nullptr, ldck_zd,      
                                   dk_zd_out.data(), lddk_zd, 
                                   rcond_zd_out.data(), tol_zd,
                                   0 /* row_major = false */);
    EXPECT_EQ(info_result, 0); 
}

