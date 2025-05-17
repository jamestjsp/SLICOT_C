#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm> // For std::max
#include <iostream>

#include "sb10hd.h"       // Include the wrapper header
#include "slicot_utils.h" // For transpose functions

// --- Test Fixture Base ---
class SB10HDTest : public ::testing::Test {
protected:
    int N_param = 6, M_param = 5, NP_param = 5, NCON_param = 2, NMEAS_param = 2;
    double TOL_param = 1e-8; // From example: 0.00000001

    double check_tol = 5.0e-5; // Verification tolerance ADJUSTED from 1e-5

    std::vector<double> A_data, B_data, C_data, D_data;
    std::vector<double> AK_out, BK_out, CK_out, DK_out;
    std::vector<double> RCOND_out;

    std::vector<double> AK_expected, BK_expected, CK_expected, DK_expected;
    std::vector<double> RCOND_expected;
    int expected_info = 0;

    int LDA, LDB, LDC, LDD;
    int LDAK, LDBK, LDCK, LDDK;

    // Helper to set D11 block to zero for a row-major matrix
    void set_d11_zero_row_major(std::vector<double>& d_matrix, int np, int m, int ncon, int nmeas) {
        int m1 = m - ncon; // Columns in D11, D21
        int np1 = np - nmeas; // Rows in D11, D12
        if (np1 > 0 && m1 > 0) { // D11 exists
            for (int i = 0; i < np1; ++i) {
                for (int j = 0; j < m1; ++j) {
                    d_matrix[i * m + j] = 0.0; // D is NPxM, row-major ld is M
                }
            }
        }
    }
     // Helper to set D11 block to zero for a col-major matrix
    void set_d11_zero_col_major(std::vector<double>& d_matrix_cm, int np, int m, int ncon, int nmeas) {
        int m1 = m - ncon; 
        int np1 = np - nmeas;
        if (np1 > 0 && m1 > 0) { 
            for (int j = 0; j < m1; ++j) { // Iterate columns of D11
                for (int i = 0; i < np1; ++i) { // Iterate rows of D11
                    d_matrix_cm[i + j * np] = 0.0; // D_cm is NPxM, col-major ld is NP
                }
            }
        }
    }


    void SetUpBase() {
        // Data from SB10HD.html example
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
        // D from example. D11 must be zero for SB10HD.
        // D is NPxM (5x5). NCON=2, NMEAS=2.
        // M1 = M-NCON = 5-2 = 3. NP1 = NP-NMEAS = 5-2 = 3.
        // So D11 is the top-left 3x3 block of D.
        D_data = {
             0.0,  0.0,  0.0, -4.0, -1.0, // Row 1 (D11, D12)
             0.0,  0.0,  0.0,  1.0,  0.0, // Row 2 (D11, D12)
             0.0,  0.0,  0.0,  0.0,  1.0, // Row 3 (D11, D12)
             3.0,  1.0,  0.0,  1.0, -3.0, // Row 4 (D21, D22)
            -2.0,  0.0,  1.0,  7.0,  1.0  // Row 5 (D21, D22)
        };
        // Ensure D11 is zero, as per SB10HD requirement
        set_d11_zero_row_major(D_data, NP_param, M_param, NCON_param, NMEAS_param);


        AK_expected = { // These values are taken from SB10HD.html, rounded to 4 decimal places there.
             88.0015, -145.7298,  -46.2424,   82.2168,  -45.2996,  -31.1407,
             25.7489,  -31.4642,  -12.4198,    9.4625,   -3.5182,    2.7056,
             54.3008, -102.4013,  -41.4968,   50.8412,  -20.1286,  -26.7191,
            108.1006, -198.0785,  -45.4333,   70.3962,  -25.8591,  -37.2741,
           -115.8900,  226.1843,   47.2549,  -47.8435,  -12.5004,   34.7474,
             59.0362, -101.8471,  -20.1052,   36.7834,  -16.1063,  -26.4309
        };
        BK_expected = {
            3.7345,    3.4758,
           -0.3020,    0.6530,
            3.4735,    4.0499,
            4.3198,    7.2755,
           -3.9424,  -10.5942,
            2.1784,    2.5048
        };
        CK_expected = {
           -2.3346,    3.2556,    0.7150,   -0.9724,    0.6962,    0.4074,
            7.6899,   -8.4558,   -2.9642,    7.0365,   -4.2844,    0.1390
        };
        DK_expected = { 
            0.0000,    0.0000,
            0.0000,    0.0000
        };
        RCOND_expected = {0.23570e0, 0.26726e0, 0.22747e-1, 0.21130e-2};

        AK_out.resize(std::max(1, N_param * N_param));
        BK_out.resize(std::max(1, N_param * NMEAS_param));
        CK_out.resize(std::max(1, NCON_param * N_param));
        DK_out.resize(std::max(1, NCON_param * NMEAS_param));
        RCOND_out.resize(4);
    }
};

class SB10HDTestColMajor : public SB10HDTest {
protected:
    std::vector<double> A_cm, B_cm, C_cm, D_cm; 
    void SetUp() override {
        SetUpBase();
        LDA = std::max(1, N_param); LDB = std::max(1, N_param); LDC = std::max(1, NP_param); LDD = std::max(1, NP_param);
        LDAK = std::max(1, N_param); LDBK = std::max(1, N_param); LDCK = std::max(1, NCON_param); LDDK = std::max(1, NCON_param);

        if (N_param > 0) { A_cm.resize(N_param*N_param); slicot_transpose_to_fortran_with_ld(A_data.data(), A_cm.data(), N_param, N_param, N_param, LDA, sizeof(double));}
        if (N_param > 0 && M_param > 0) { B_cm.resize(N_param*M_param); slicot_transpose_to_fortran_with_ld(B_data.data(), B_cm.data(), N_param, M_param, M_param, LDB, sizeof(double));}
        if (NP_param > 0 && N_param > 0) { C_cm.resize(NP_param*N_param); slicot_transpose_to_fortran_with_ld(C_data.data(), C_cm.data(), NP_param, N_param, N_param, LDC, sizeof(double));}
        if (NP_param > 0 && M_param > 0) { D_cm.resize(NP_param*M_param); slicot_transpose_to_fortran_with_ld(D_data.data(), D_cm.data(), NP_param, M_param, M_param, LDD, sizeof(double));}
        // D_data (row-major) already has D11 zeroed in SetUpBase.
        // So, D_cm (column-major) will also have its D11 block as zero after transposition.
        // The set_d11_zero_col_major call is thus not strictly needed here if D_data is correctly prepared.
        // set_d11_zero_col_major(D_cm, NP_param, M_param, NCON_param, NMEAS_param); 
    }
};

class SB10HDTestRowMajor : public SB10HDTest {
protected:
    void SetUp() override {
        SetUpBase(); // D_data already has D11 zeroed
        LDA = std::max(1, N_param); LDB = std::max(1, M_param); LDC = std::max(1, N_param); LDD = std::max(1, M_param);
        LDAK = std::max(1, N_param); LDBK = std::max(1, NMEAS_param); LDCK = std::max(1, N_param); LDDK = std::max(1, NMEAS_param);
    }
};

TEST_F(SB10HDTestColMajor, DocExampleColMajor) {
    int info_result = slicot_sb10hd(N_param, M_param, NP_param, NCON_param, NMEAS_param,
                                   A_cm.empty() ? nullptr : A_cm.data(), LDA,
                                   B_cm.empty() ? nullptr : B_cm.data(), LDB,
                                   C_cm.empty() ? nullptr : C_cm.data(), LDC,
                                   D_cm.empty() ? nullptr : D_cm.data(), LDD,
                                   AK_out.data(), LDAK, BK_out.data(), LDBK,
                                   CK_out.data(), LDCK, DK_out.data(), LDDK,
                                   RCOND_out.data(), TOL_param,
                                   0 /* row_major = false */);
    ASSERT_EQ(info_result, expected_info);
    std::vector<double> AK_exp_cm, BK_exp_cm, CK_exp_cm, DK_exp_cm;
    if(!AK_expected.empty()){ AK_exp_cm.resize(AK_expected.size()); slicot_transpose_to_fortran_with_ld(AK_expected.data(), AK_exp_cm.data(), N_param, N_param, N_param, LDAK, sizeof(double)); for(size_t i=0; i<AK_out.size(); ++i) EXPECT_NEAR(AK_out[i], AK_exp_cm[i], check_tol) << "AK["<<i<<"]";}
    if(!BK_expected.empty()){ BK_exp_cm.resize(BK_expected.size()); slicot_transpose_to_fortran_with_ld(BK_expected.data(), BK_exp_cm.data(), N_param, NMEAS_param, NMEAS_param, LDBK, sizeof(double)); for(size_t i=0; i<BK_out.size(); ++i) EXPECT_NEAR(BK_out[i], BK_exp_cm[i], check_tol) << "BK["<<i<<"]";}
    if(!CK_expected.empty()){ CK_exp_cm.resize(CK_expected.size()); slicot_transpose_to_fortran_with_ld(CK_expected.data(), CK_exp_cm.data(), NCON_param, N_param, N_param, LDCK, sizeof(double)); for(size_t i=0; i<CK_out.size(); ++i) EXPECT_NEAR(CK_out[i], CK_exp_cm[i], check_tol) << "CK["<<i<<"]";}
    if(!DK_expected.empty()){ DK_exp_cm.resize(DK_expected.size()); slicot_transpose_to_fortran_with_ld(DK_expected.data(), DK_exp_cm.data(), NCON_param, NMEAS_param, NMEAS_param, LDDK, sizeof(double)); for(size_t i=0; i<DK_out.size(); ++i) EXPECT_NEAR(DK_out[i], DK_exp_cm[i], check_tol) << "DK["<<i<<"]";}
    for(size_t i=0; i<RCOND_out.size(); ++i) EXPECT_NEAR(RCOND_out[i], RCOND_expected[i], check_tol) << "RCOND["<<i<<"]";
}

TEST_F(SB10HDTestRowMajor, DocExampleRowMajor) {
    int info_result = slicot_sb10hd(N_param, M_param, NP_param, NCON_param, NMEAS_param,
                                   A_data.empty() ? nullptr : A_data.data(), LDA,
                                   B_data.empty() ? nullptr : B_data.data(), LDB,
                                   C_data.empty() ? nullptr : C_data.data(), LDC,
                                   D_data.empty() ? nullptr : D_data.data(), LDD,
                                   AK_out.data(), LDAK, BK_out.data(), LDBK,
                                   CK_out.data(), LDCK, DK_out.data(), LDDK,
                                   RCOND_out.data(), TOL_param,
                                   1 /* row_major = true */);
    ASSERT_EQ(info_result, expected_info);
    if(!AK_expected.empty()) for(size_t i=0; i<AK_out.size(); ++i) EXPECT_NEAR(AK_out[i], AK_expected[i], check_tol) << "AK["<<i<<"]";
    if(!BK_expected.empty()) for(size_t i=0; i<BK_out.size(); ++i) EXPECT_NEAR(BK_out[i], BK_expected[i], check_tol) << "BK["<<i<<"]";
    if(!CK_expected.empty()) for(size_t i=0; i<CK_out.size(); ++i) EXPECT_NEAR(CK_out[i], CK_expected[i], check_tol) << "CK["<<i<<"]";
    if(!DK_expected.empty()) for(size_t i=0; i<DK_out.size(); ++i) EXPECT_NEAR(DK_out[i], DK_expected[i], check_tol) << "DK["<<i<<"]";
    for(size_t i=0; i<RCOND_out.size(); ++i) EXPECT_NEAR(RCOND_out[i], RCOND_expected[i], check_tol) << "RCOND["<<i<<"]";
}

TEST_F(SB10HDTestColMajor, ParameterValidation) {
    std::vector<double> dummy_rcond(4);
    int info_result = slicot_sb10hd(-1, M_param, NP_param, NCON_param, NMEAS_param, nullptr,1,nullptr,1,nullptr,1,nullptr,1,nullptr,1,nullptr,1,nullptr,1,nullptr,1,dummy_rcond.data(),TOL_param,0);
    EXPECT_EQ(info_result, -1); // Invalid N
    info_result = slicot_sb10hd(N_param, M_param, NP_param, M_param + 1, NMEAS_param, A_cm.data(),LDA,B_cm.data(),LDB,C_cm.data(),LDC,D_cm.data(),LDD,AK_out.data(),LDAK,BK_out.data(),LDBK,CK_out.data(),LDCK,DK_out.data(),LDDK,dummy_rcond.data(),TOL_param,0);
    EXPECT_EQ(info_result, -4); // NCON > M
    if (N_param > 0) {
      info_result = slicot_sb10hd(N_param,M_param,NP_param,NCON_param,NMEAS_param, A_cm.data(),0,B_cm.data(),LDB,C_cm.data(),LDC,D_cm.data(),LDD,AK_out.data(),LDAK,BK_out.data(),LDBK,CK_out.data(),LDCK,DK_out.data(),LDDK,dummy_rcond.data(),TOL_param,0);
      EXPECT_EQ(info_result, -7); // Invalid LDA
    }
}

TEST_F(SB10HDTestColMajor, ZeroDimensionN) {
    int n_zero = 0;
    int m_zd = 2, np_zd = 2, ncon_zd = 1, nmeas_zd = 1; // M, NP, NCON, NMEAS
    double tol_zd = 1e-8;

    // D is NPxM = 2x2. DK is NCONxNMEAS = 1x1.
    // D11 must be zero. M1 = M-NCON = 1, NP1 = NP-NMEAS = 1. So D11 is D(0,0).
    // D12 (NP1xNCON = 1x1) is D(0,1)
    // D21 (NMEASxM1 = 1x1) is D(1,0)
    // D22 (NMEASxNCON = 1x1) is D(1,1)
    // For column major: D = {D(0,0), D(1,0), D(0,1), D(1,1)}
    // To make D12 and D21 non-zero (full rank):
    std::vector<double> d_zd_data = {0.0, 1.0, 1.0, 1.0}; // D(0,0)=0 (D11), D(1,0)=1 (D21), D(0,1)=1 (D12), D(1,1)=1 (D22)
    
    std::vector<double> dk_zd_out(std::max(1,ncon_zd * nmeas_zd));
    std::vector<double> rcond_zd_out(4);

    int ldd_zd = std::max(1, np_zd); 
    int lddk_zd = std::max(1, ncon_zd);

    int info_result = slicot_sb10hd(n_zero, m_zd, np_zd, ncon_zd, nmeas_zd,
                                   nullptr, 1, nullptr, 1, nullptr, std::max(1,np_zd),
                                   d_zd_data.data(), ldd_zd,
                                   nullptr, 1, nullptr, 1,
                                   nullptr, std::max(1,ncon_zd), dk_zd_out.data(), lddk_zd,
                                   rcond_zd_out.data(), tol_zd,
                                   0 /* row_major = false */);
    EXPECT_EQ(info_result, 0);
}

