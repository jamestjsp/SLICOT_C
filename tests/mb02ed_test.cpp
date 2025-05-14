#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include "mb02ed.h"
#include "slicot_utils.h"

class MB02EDTestColMajor : public ::testing::Test {
protected:
    // Example parameters from documentation
    int K = 3;
    int N = 3;
    int NRHS = 2;
    char TYPET = 'C';
    int LDT = 9;   // Fortran: LDT >= N*K = 9
    int LDB = 9;   // Fortran: LDB >= N*K = 9
    int row_major = 0;
    double tol = 1e-4;

    // Input arrays (column-major)
    std::vector<double> T; // (N*K) x K
    std::vector<double> B; // (N*K) x NRHS

    // Expected solution (from doc output, column-major)
    std::vector<double> X_expected;

    void SetUp() override {
        // T: (N*K) x K = 9 x 3, column-major
        T = {
            3.0000, 1.0000, 0.2000, 0.1000, 0.2000, 0.0500, 0.1000, 0.0400, 0.0100, // col 1
            1.0000, 4.0000, 0.4000, 0.1000, 0.0400, 0.2000, 0.0300, 0.0200, 0.0300, // col 2
            0.2000, 0.4000, 5.0000, 0.2000, 0.0300, 0.1000, 0.1000, 0.2000, 0.0200  // col 3
        };
        // B: (N*K) x NRHS = 9 x 2, column-major
        B = {
            1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, // col 1
            2.0000, 2.0000, 2.0000, 2.0000, 2.0000, 2.0000, 2.0000, 2.0000, 2.0000  // col 2
        };
        // X_expected: (N*K) x NRHS, column-major
        X_expected = {
            0.2408, 0.1558, 0.1534, 0.2302, 0.1467, 0.1537, 0.2349, 0.1498, 0.1653, // col 1
            0.4816, 0.3116, 0.3068, 0.4603, 0.2934, 0.3075, 0.4698, 0.2995, 0.3307  // col 2
        };
    }
};

TEST_F(MB02EDTestColMajor, DocExample_TypetC) {
    std::vector<double> T_in = T; // Copy, as T is in/out
    std::vector<double> B_in = B; // Copy, as B is in/out

    int info = slicot_mb02ed(TYPET, K, N, NRHS, T_in.data(), LDT, B_in.data(), LDB, row_major);

    ASSERT_EQ(info, 0);

    // Check B_in (solution X) matches expected
    for (size_t i = 0; i < X_expected.size(); ++i) {
        EXPECT_NEAR(B_in[i], X_expected[i], tol) << "Mismatch at index " << i;
    }
}

class MB02EDTestRowMajor : public MB02EDTestColMajor {
protected:
    int row_major_rm = 1;
    int LDT_rm = 3; // For row-major, LDT is number of columns (K)
    int LDB_rm = 2; // For row-major, LDB is number of columns (NRHS)
    std::vector<double> T_rm; // (N*K) rows x K cols, row-major
    std::vector<double> B_rm; // (N*K) rows x NRHS cols, row-major
    std::vector<double> X_expected_rm;

    void SetUp() override {
        MB02EDTestColMajor::SetUp();
        // Convert T and B to row-major
        T_rm.resize(T.size());
        B_rm.resize(B.size());
        slicot_transpose_to_c(T.data(), T_rm.data(), N*K, K, sizeof(double));
        slicot_transpose_to_c(B.data(), B_rm.data(), N*K, NRHS, sizeof(double));
        // Convert expected X to row-major for comparison
        X_expected_rm.resize(X_expected.size());
        slicot_transpose_to_c(X_expected.data(), X_expected_rm.data(), N*K, NRHS, sizeof(double));
    }
};

TEST_F(MB02EDTestRowMajor, DocExample_TypetC_RowMajor) {
    std::vector<double> T_in = T_rm;
    std::vector<double> B_in = B_rm;

    int info = slicot_mb02ed(TYPET, K, N, NRHS, T_in.data(), LDT_rm, B_in.data(), LDB_rm, row_major_rm);

    ASSERT_EQ(info, 0);

    // Check B_in (solution X) matches expected (row-major)
    for (size_t i = 0; i < X_expected_rm.size(); ++i) {
        EXPECT_NEAR(B_in[i], X_expected_rm[i], tol) << "Mismatch at index " << i;
    }
}

TEST_F(MB02EDTestColMajor, ParameterValidation) {
    std::vector<double> dummy(1, 0.0);

    // Invalid K
    int info = slicot_mb02ed(TYPET, -1, N, NRHS, dummy.data(), LDT, dummy.data(), LDB, 0);
    EXPECT_EQ(info, -2);

    // Invalid N
    info = slicot_mb02ed(TYPET, K, -1, NRHS, dummy.data(), LDT, dummy.data(), LDB, 0);
    EXPECT_EQ(info, -3);

    // Invalid NRHS
    info = slicot_mb02ed(TYPET, K, N, -1, dummy.data(), LDT, dummy.data(), LDB, 0);
    EXPECT_EQ(info, -4);

    // Invalid TYPET
    info = slicot_mb02ed('Z', K, N, NRHS, dummy.data(), LDT, dummy.data(), LDB, 0);
    EXPECT_EQ(info, -1);

    // Invalid LDT (col-major)
    info = slicot_mb02ed(TYPET, K, N, NRHS, dummy.data(), 0, dummy.data(), LDB, 0);
    EXPECT_EQ(info, -6);

    // Invalid LDB (col-major)
    info = slicot_mb02ed(TYPET, K, N, NRHS, dummy.data(), LDT, dummy.data(), 0, 0);
    EXPECT_EQ(info, -8);

    // Invalid LDT (row-major)
    info = slicot_mb02ed(TYPET, K, N, NRHS, dummy.data(), 0, dummy.data(), LDB, 1);
    EXPECT_EQ(info, -6);

    // Invalid LDB (row-major)
    info = slicot_mb02ed(TYPET, K, N, NRHS, dummy.data(), LDT, dummy.data(), 0, 1);
    EXPECT_EQ(info, -8);
}
