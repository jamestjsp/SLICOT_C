#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm> // For std::max
#include <stdexcept> // For std::runtime_error

#include "mb05nd.h"
#include "slicot_utils.h"

// --- Column-Major Test Fixture ---
class MB05NDTestColMajor : public ::testing::Test {
protected:
    int N = 5;                // Order of the matrix A
    double DELTA = 0.1;       // Scalar delta
    double TOL = 1e-4;        // Tolerance for Pade approximation
    int INFO = -999;          // Result from slicot_mb05nd call

    // Input matrix A (column-major format)
    std::vector<double> A;

    // Output matrices EX and EXINT (column-major format)
    std::vector<double> EX;
    std::vector<double> EXINT;

    // Expected results for EX and EXINT
    std::vector<double> EX_expected;
    std::vector<double> EXINT_expected;

    // Leading dimensions
    int LDA = 0;
    int LDEX = 0;
    int LDEXIN = 0;

    // Verification tolerance
    double check_tol = 1e-4;

    void SetUp() override {
        // Set leading dimensions
        LDA = std::max(1, N);
        LDEX = std::max(1, N);
        LDEXIN = std::max(1, N);

        // Initialize matrix A (column-major format)
        A = {
            5.0, 1.0, 2.0, 1.0, 2.0,
            4.0, 6.0, 0.0, 3.0, 5.0,
            3.0, 0.0, 7.0, 1.0, 7.0,
            2.0, 4.0, 6.0, 8.0, 1.0,
            1.0, 3.0, 5.0, 7.0, 9.0
        };

        // Resize output matrices
        EX.resize((size_t)LDEX * N);
        EXINT.resize((size_t)LDEXIN * N);

        // Expected results for EX and EXINT (from MB05ND documentation example)
        EX_expected = {
            1.8391, 0.3359, 0.6335, 0.4804, 0.7105,
            0.9476, 2.2262, 0.6776, 1.1561, 1.4244,
            0.7920, 0.4013, 2.6933, 0.9110, 1.8835,
            0.8216, 1.0078, 1.6155, 2.7461, 1.0966,
            0.7811, 1.0957, 1.8502, 2.0854, 3.4134
        };

        EXINT_expected = {
            0.1347, 0.0114, 0.0218, 0.0152, 0.0240,
            0.0352, 0.1477, 0.0178, 0.0385, 0.0503,
            0.0284, 0.0104, 0.1624, 0.0267, 0.0679,
            0.0272, 0.0369, 0.0580, 0.1660, 0.0317,
            0.0231, 0.0368, 0.0619, 0.0732, 0.1863
        };
    }
};

// --- Row-Major Test Fixture ---
class MB05NDTestRowMajor : public MB05NDTestColMajor {
protected:
    std::vector<double> A_rm;
    std::vector<double> EX_rm;
    std::vector<double> EXINT_rm;

    std::vector<double> EX_expected_rm;
    std::vector<double> EXINT_expected_rm;

    void SetUp() override {
        MB05NDTestColMajor::SetUp();

        // Resize row-major matrices
        A_rm.resize((size_t)N * LDA);
        EX_rm.resize((size_t)N * LDEX);
        EXINT_rm.resize((size_t)N * LDEXIN);

        // Transpose column-major inputs to row-major
        slicot_transpose_to_c(A.data(), A_rm.data(), N, N, sizeof(double));

        // Transpose expected results to row-major
        EX_expected_rm.resize((size_t)N * N);
        EXINT_expected_rm.resize((size_t)N * N);
        slicot_transpose_to_c(EX_expected.data(), EX_expected_rm.data(), N, N, sizeof(double));
        slicot_transpose_to_c(EXINT_expected.data(), EXINT_expected_rm.data(), N, N, sizeof(double));
    }
};

// --- Test Cases ---

// Test: Documentation Example (Column-Major)
TEST_F(MB05NDTestColMajor, DocExample) {
    INFO = slicot_mb05nd(N, DELTA, A.data(), LDA, EX.data(), LDEX, EXINT.data(), LDEXIN, TOL, 0 /* column-major */);

    ASSERT_EQ(INFO, 0) << "slicot_mb05nd call failed with INFO = " << INFO;

    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            EXPECT_NEAR(EX[j * LDEX + i], EX_expected[j * N + i], check_tol)
                << "EX mismatch at (" << i << "," << j << ")";
            EXPECT_NEAR(EXINT[j * LDEXIN + i], EXINT_expected[j * N + i], check_tol)
                << "EXINT mismatch at (" << i << "," << j << ")";
        }
    }
}

// Test: Documentation Example (Row-Major)
TEST_F(MB05NDTestRowMajor, DocExample) {
    INFO = slicot_mb05nd(N, DELTA, A_rm.data(), LDA, EX_rm.data(), LDEX, EXINT_rm.data(), LDEXIN, TOL, 1 /* row-major */);

    ASSERT_EQ(INFO, 0) << "slicot_mb05nd call failed with INFO = " << INFO;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            EXPECT_NEAR(EX_rm[i * LDEX + j], EX_expected_rm[i * N + j], check_tol)
                << "EX_rm mismatch at (" << i << "," << j << ")";
            EXPECT_NEAR(EXINT_rm[i * LDEXIN + j], EXINT_expected_rm[i * N + j], check_tol)
                << "EXINT_rm mismatch at (" << i << "," << j << ")";
        }
    }
}

// Test: Zero Dimensions
TEST_F(MB05NDTestColMajor, ZeroDimensions) {
    int n_zero = 0;
    int lda_z = 1, ldex_z = 1, ldexin_z = 1;
    double delta_z = 0.1, tol_z = 1e-4;

    INFO = slicot_mb05nd(n_zero, delta_z, nullptr, lda_z, nullptr, ldex_z, nullptr, ldexin_z, tol_z, 0 /* column-major */);

    EXPECT_EQ(INFO, 0) << "slicot_mb05nd failed for zero dimensions with INFO = " << INFO;
}

// Test: Parameter Validation
TEST_F(MB05NDTestColMajor, ParameterValidation) {
    // Test invalid N
    INFO = slicot_mb05nd(-1, DELTA, A.data(), LDA, EX.data(), LDEX, EXINT.data(), LDEXIN, TOL, 0);
    EXPECT_EQ(INFO, -1);

    // Test invalid LDA
    INFO = slicot_mb05nd(N, DELTA, A.data(), 0, EX.data(), LDEX, EXINT.data(), LDEXIN, TOL, 0);
    EXPECT_EQ(INFO, -4);

    // Test invalid LDEX
    INFO = slicot_mb05nd(N, DELTA, A.data(), LDA, EX.data(), 0, EXINT.data(), LDEXIN, TOL, 0);
    EXPECT_EQ(INFO, -6);

    // Test invalid LDEXIN
    INFO = slicot_mb05nd(N, DELTA, A.data(), LDA, EX.data(), LDEX, EXINT.data(), 0, TOL, 0);
    EXPECT_EQ(INFO, -8);
}
