#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <string>
#include <stdexcept>
#include <algorithm>

#include "sb10dd.h"
#include "slicot_utils.h"
#include "test_config.h"

class SB10DDTest : public ::testing::Test {
protected:
    int N, M, NP, NCON, NMEAS;
    double GAMMA, TOL;

    std::vector<double> A_in, B_in, C_in, D_in;
    std::vector<double> AK_expected, BK_expected, CK_expected, DK_expected;
    std::vector<double> RCOND_expected;

    std::vector<double> AK_out, BK_out, CK_out, DK_out;
    std::vector<double> X_out, Z_out;
    std::vector<double> RCOND_out;

    int LDA, LDB, LDC, LDD, LDAK, LDBK, LDCK, LDDK, LDX, LDZ;

    int INFO_expected = 0;
    int INFO_result = -999;
    double check_tol = 1e-4;

    // Method to initialize data directly in the fixture
    void InitializeData(bool row_major_c_storage) {
        N = 6; M = 5; NP = 5; NCON = 2; NMEAS = 2;
        GAMMA = 111.294; TOL = 0.00000001;

        // Data from SB10DD.dat / SB10DD.res, stored in column-major order first
        std::vector<double> A_colmajor = {
            -0.7, -0.6, -0.5, -0.7,  0.0,  0.5,
             0.0,  0.2,  0.7,  0.0,  0.3, -0.8,
             0.3, -0.4, -0.1,  0.0,  0.6,  0.0,
             0.0, -0.3,  0.0, -0.5, -0.9,  0.0,
            -0.5,  0.0,  0.0, -1.0,  0.1,  0.2,
            -0.1,  0.0, -0.8,  0.0, -0.4, -0.9
        };
        std::vector<double> B_colmajor = {
            -1.0,  1.0, -3.0,  1.0,  0.0,  1.0,
            -2.0,  0.0, -4.0, -2.0,  1.0,  0.0,
            -2.0,  1.0,  0.0,  1.0, -2.0,  3.0,
             1.0, -2.0,  2.0,  0.0,  0.0, -1.0,
             0.0,  1.0, -2.0, -1.0,  3.0, -2.0
        };
        std::vector<double> C_colmajor = {
             1.0, -3.0,  0.0,  1.0,  0.0,
            -1.0,  0.0,  2.0, -3.0,  1.0,
             2.0,  1.0,  0.0,  0.0, -2.0,
            -2.0, -1.0, -4.0,  0.0,  1.0,
             0.0,  1.0,  0.0,  3.0,  0.0,
            -3.0,  0.0, -2.0,  1.0, -2.0
        };
        std::vector<double> D_colmajor = {
             1.0,  0.0,  2.0,  0.0,  0.0,
            -1.0,  1.0, -1.0,  1.0,  0.0,
            -2.0,  0.0, -3.0,  0.0,  1.0,
             0.0,  1.0,  0.0,  1.0,  2.0,
             0.0,  0.0,  1.0, -1.0,  1.0
        };
        std::vector<double> AK_expected_colmajor = {
            -18.0030,  18.8203, -26.5994, -21.4163,  -0.8911,  -5.3286,
             52.0376, -57.6244,  77.9693,  62.1719,   4.2787,  16.1955,
             26.0831, -29.0938,  39.0368,  30.7507,   2.3286,   8.4824,
             -0.4271,   0.5870,  -1.4020,  -0.9201,  -0.2424,  -0.2489,
            -40.9022,  45.3309, -60.1129, -48.6221,  -3.0376, -12.2348,
             18.0857, -19.8644,  26.6910,  21.8351,   1.2169,   5.1590
        };
        std::vector<double> BK_expected_colmajor = {
             16.9788, -18.9215,  25.2046,  20.1122,   1.4104,   5.3181,
             14.1648, -15.6726,  21.2848,  16.8322,   1.2040,   4.5149
        };
        std::vector<double> CK_expected_colmajor = {
            -9.1941,  3.6490,
            27.5165,-10.6194,
            13.7364, -5.2772,
            -0.3639,  0.2432,
           -21.5983,  8.1108,
             9.6025, -3.6293
        };
        std::vector<double> DK_expected_colmajor = {
            9.0317, -3.4006,
            7.5348, -2.8219
        };
        RCOND_expected = {0.24960, 0.98548, 0.99186, 0.0000063733, 0.48625, 0.029430, 0.0056942, 0.012470};

        if (row_major_c_storage) {
            A_in.resize(A_colmajor.size()); slicot_transpose_to_c(A_colmajor.data(), A_in.data(), N, N, sizeof(double));
            B_in.resize(B_colmajor.size()); slicot_transpose_to_c(B_colmajor.data(), B_in.data(), N, M, sizeof(double));
            C_in.resize(C_colmajor.size()); slicot_transpose_to_c(C_colmajor.data(), C_in.data(), NP, N, sizeof(double));
            D_in.resize(D_colmajor.size()); slicot_transpose_to_c(D_colmajor.data(), D_in.data(), NP, M, sizeof(double));
            AK_expected.resize(AK_expected_colmajor.size()); slicot_transpose_to_c(AK_expected_colmajor.data(), AK_expected.data(), N, N, sizeof(double));
            BK_expected.resize(BK_expected_colmajor.size()); slicot_transpose_to_c(BK_expected_colmajor.data(), BK_expected.data(), N, NMEAS, sizeof(double));
            CK_expected.resize(CK_expected_colmajor.size()); slicot_transpose_to_c(CK_expected_colmajor.data(), CK_expected.data(), NCON, N, sizeof(double));
            DK_expected.resize(DK_expected_colmajor.size()); slicot_transpose_to_c(DK_expected_colmajor.data(), DK_expected.data(), NCON, NMEAS, sizeof(double));

            LDA = std::max(1,N); LDB = std::max(1,M); LDC = std::max(1,N); LDD = std::max(1,M);
            LDAK = std::max(1,N); LDBK = std::max(1,NMEAS); LDCK = std::max(1,N); LDDK = std::max(1,NMEAS);
            LDX = std::max(1,N); LDZ = std::max(1,N);
        } else {
            A_in = A_colmajor;
            B_in = B_colmajor;
            C_in = C_colmajor;
            D_in = D_colmajor;
            AK_expected = AK_expected_colmajor;
            BK_expected = BK_expected_colmajor;
            CK_expected = CK_expected_colmajor;
            DK_expected = DK_expected_colmajor;

            LDA = std::max(1,N); LDB = std::max(1,N); LDC = std::max(1,NP); LDD = std::max(1,NP);
            LDAK = std::max(1,N); LDBK = std::max(1,N); LDCK = std::max(1,NCON); LDDK = std::max(1,NCON);
            LDX = std::max(1,N); LDZ = std::max(1,N);
        }

        AK_out.resize((size_t)LDAK * N);
        BK_out.resize((size_t)LDBK * NMEAS);
        CK_out.resize((size_t)LDCK * N);
        DK_out.resize((size_t)LDDK * NMEAS);
        X_out.resize((size_t)LDX * N);
        Z_out.resize((size_t)LDZ * N);
        RCOND_out.resize(8);
    }
};

class SB10DDTestColMajor : public SB10DDTest {
protected:
    void SetUp() override { InitializeData(false); }
};

class SB10DDTestRowMajor : public SB10DDTest {
protected:
    void SetUp() override { InitializeData(true); }
};

TEST_F(SB10DDTestColMajor, DocExample) {
    INFO_result = slicot_sb10dd(N, M, NP, NCON, NMEAS, GAMMA,
                  A_in.data(), LDA, B_in.data(), LDB, C_in.data(), LDC, D_in.data(), LDD,
                  AK_out.data(), LDAK, BK_out.data(), LDBK, CK_out.data(), LDCK, DK_out.data(), LDDK,
                  X_out.data(), LDX, Z_out.data(), LDZ, RCOND_out.data(), TOL,
                  0 /* row_major = false */);

    ASSERT_EQ(INFO_result, INFO_expected);
    if (INFO_result == 0) {
        for(size_t i=0; i < AK_expected.size(); ++i) EXPECT_NEAR(AK_out[i], AK_expected[i], check_tol) << "AK mismatch at index " << i;
        for(size_t i=0; i < BK_expected.size(); ++i) EXPECT_NEAR(BK_out[i], BK_expected[i], check_tol) << "BK mismatch at index " << i;
        for(size_t i=0; i < CK_expected.size(); ++i) EXPECT_NEAR(CK_out[i], CK_expected[i], check_tol) << "CK mismatch at index " << i;
        for(size_t i=0; i < DK_expected.size(); ++i) EXPECT_NEAR(DK_out[i], DK_expected[i], check_tol) << "DK mismatch at index " << i;
        for(size_t i=0; i < RCOND_expected.size(); ++i) EXPECT_NEAR(RCOND_out[i], RCOND_expected[i], check_tol) << "RCOND mismatch at index " << i;
    }
}

TEST_F(SB10DDTestRowMajor, DocExample) {
     INFO_result = slicot_sb10dd(N, M, NP, NCON, NMEAS, GAMMA,
                  A_in.data(), LDA, B_in.data(), LDB, C_in.data(), LDC, D_in.data(), LDD,
                  AK_out.data(), LDAK, BK_out.data(), LDBK, CK_out.data(), LDCK, DK_out.data(), LDDK,
                  X_out.data(), LDX, Z_out.data(), LDZ, RCOND_out.data(), TOL,
                  1 /* row_major = true */);

    ASSERT_EQ(INFO_result, INFO_expected);
    if (INFO_result == 0) {
        // When row_major = 1, the wrapper slicot_sb10dd is expected to return outputs
        // (AK_out, BK_out, CK_out, DK_out) in row-major format.
        // The expected arrays (AK_expected, etc.) are already in row-major format
        // as prepared by InitializeData(true).
        // Therefore, direct comparison is needed.

        for(size_t i=0; i < AK_expected.size(); ++i) EXPECT_NEAR(AK_out[i], AK_expected[i], check_tol) << "AK_rm mismatch at index " << i;
        for(size_t i=0; i < BK_expected.size(); ++i) EXPECT_NEAR(BK_out[i], BK_expected[i], check_tol) << "BK_rm mismatch at index " << i;
        for(size_t i=0; i < CK_expected.size(); ++i) EXPECT_NEAR(CK_out[i], CK_expected[i], check_tol) << "CK_rm mismatch at index " << i;
        for(size_t i=0; i < DK_expected.size(); ++i) EXPECT_NEAR(DK_out[i], DK_expected[i], check_tol) << "DK_rm mismatch at index " << i;
        for(size_t i=0; i < RCOND_expected.size(); ++i) EXPECT_NEAR(RCOND_out[i], RCOND_expected[i], check_tol) << "RCOND mismatch at index " << i;
    }
}

TEST_F(SB10DDTestColMajor, ParameterValidation) {
    // Test invalid N
    INFO_result = slicot_sb10dd(-1, M, NP, NCON, NMEAS, GAMMA, A_in.data(), LDA, B_in.data(), LDB, C_in.data(), LDC, D_in.data(), LDD, AK_out.data(), LDAK, BK_out.data(), LDBK, CK_out.data(), LDCK, DK_out.data(), LDDK, X_out.data(), LDX, Z_out.data(), LDZ, RCOND_out.data(), TOL, 0);
    EXPECT_EQ(INFO_result, -1);

    // Test invalid M
    INFO_result = slicot_sb10dd(N, -1, NP, NCON, NMEAS, GAMMA, A_in.data(), LDA, B_in.data(), LDB, C_in.data(), LDC, D_in.data(), LDD, AK_out.data(), LDAK, BK_out.data(), LDBK, CK_out.data(), LDCK, DK_out.data(), LDDK, X_out.data(), LDX, Z_out.data(), LDZ, RCOND_out.data(), TOL, 0);
    EXPECT_EQ(INFO_result, -2);

    // Test invalid NP
    INFO_result = slicot_sb10dd(N, M, -1, NCON, NMEAS, GAMMA, A_in.data(), LDA, B_in.data(), LDB, C_in.data(), LDC, D_in.data(), LDD, AK_out.data(), LDAK, BK_out.data(), LDBK, CK_out.data(), LDCK, DK_out.data(), LDDK, X_out.data(), LDX, Z_out.data(), LDZ, RCOND_out.data(), TOL, 0);
    EXPECT_EQ(INFO_result, -3);

    // Test invalid NCON
    INFO_result = slicot_sb10dd(N, M, NP, -1, NMEAS, GAMMA, A_in.data(), LDA, B_in.data(), LDB, C_in.data(), LDC, D_in.data(), LDD, AK_out.data(), LDAK, BK_out.data(), LDBK, CK_out.data(), LDCK, DK_out.data(), LDDK, X_out.data(), LDX, Z_out.data(), LDZ, RCOND_out.data(), TOL, 0);
    EXPECT_EQ(INFO_result, -4);
    INFO_result = slicot_sb10dd(N, M, NP, M + 1, NMEAS, GAMMA, A_in.data(), LDA, B_in.data(), LDB, C_in.data(), LDC, D_in.data(), LDD, AK_out.data(), LDAK, BK_out.data(), LDBK, CK_out.data(), LDCK, DK_out.data(), LDDK, X_out.data(), LDX, Z_out.data(), LDZ, RCOND_out.data(), TOL, 0);
    EXPECT_EQ(INFO_result, -4);

    // Test invalid NMEAS
    INFO_result = slicot_sb10dd(N, M, NP, NCON, -1, GAMMA, A_in.data(), LDA, B_in.data(), LDB, C_in.data(), LDC, D_in.data(), LDD, AK_out.data(), LDAK, BK_out.data(), LDBK, CK_out.data(), LDCK, DK_out.data(), LDDK, X_out.data(), LDX, Z_out.data(), LDZ, RCOND_out.data(), TOL, 0);
    EXPECT_EQ(INFO_result, -5);
    INFO_result = slicot_sb10dd(N, M, NP, NCON, NP + 1, GAMMA, A_in.data(), LDA, B_in.data(), LDB, C_in.data(), LDC, D_in.data(), LDD, AK_out.data(), LDAK, BK_out.data(), LDBK, CK_out.data(), LDCK, DK_out.data(), LDDK, X_out.data(), LDX, Z_out.data(), LDZ, RCOND_out.data(), TOL, 0);
    EXPECT_EQ(INFO_result, -5);

    // Test invalid GAMMA
    INFO_result = slicot_sb10dd(N, M, NP, NCON, NMEAS, -1.0, A_in.data(), LDA, B_in.data(), LDB, C_in.data(), LDC, D_in.data(), LDD, AK_out.data(), LDAK, BK_out.data(), LDBK, CK_out.data(), LDCK, DK_out.data(), LDDK, X_out.data(), LDX, Z_out.data(), LDZ, RCOND_out.data(), TOL, 0);
    EXPECT_EQ(INFO_result, -6);

    // Test invalid LDA
    if (N > 0) {
        INFO_result = slicot_sb10dd(N, M, NP, NCON, NMEAS, GAMMA, A_in.data(), 0, B_in.data(), LDB, C_in.data(), LDC, D_in.data(), LDD, AK_out.data(), LDAK, BK_out.data(), LDBK, CK_out.data(), LDCK, DK_out.data(), LDDK, X_out.data(), LDX, Z_out.data(), LDZ, RCOND_out.data(), TOL, 0);
        EXPECT_EQ(INFO_result, -8);
    }
}

TEST_F(SB10DDTestColMajor, ZeroDimensionN) {
    N = 0;
    // Re-initialize with N=0 for this specific test
    InitializeData(false); // Or true, doesn't matter much for N=0 as arrays will be empty or NULL
    N = 0; // Override N after InitializeData

    // Update LDs for N=0. Most data pointers will be NULL.
    LDA = std::max(1,N); LDB = std::max(1,N); // LDC, LDD are not N-dependent in the same way
    LDAK = std::max(1,N); LDBK = std::max(1,N); // LDCK, LDDK are not N-dependent
    LDX = std::max(1,N); LDZ = std::max(1,N);

    // Resize output vectors to be empty or appropriately sized for N=0
    AK_out.assign((size_t)LDAK * N, 0.0);      // size 0
    BK_out.assign((size_t)LDBK * NMEAS, 0.0); // size 0 if N=0 (LDBK becomes 1, N=0)
    CK_out.assign((size_t)LDCK * N, 0.0);      // size 0
    X_out.assign((size_t)LDX * N, 0.0);        // size 0
    Z_out.assign((size_t)LDZ * N, 0.0);        // size 0
    // DK_out and RCOND_out are not directly N-dependent in their primary dimension for allocation
    // DK_out is NCON x NMEAS. RCOND is fixed size 8.
    // D_in is NP x M

    INFO_expected = 0; 
    INFO_result = slicot_sb10dd(N, M, NP, NCON, NMEAS, GAMMA,
                  nullptr, LDA, nullptr, LDB, nullptr, LDC, D_in.data(), LDD,
                  nullptr, LDAK, nullptr, LDBK, nullptr, LDCK, DK_out.data(), LDDK,
                  nullptr, LDX, nullptr, LDZ, RCOND_out.data(), TOL,
                  0);
    EXPECT_EQ(INFO_result, INFO_expected);
}

