#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm> // For std::max
#include <stdexcept> // For std::runtime_error

#include "mc01td.h"

// --- Test Fixture ---
class MC01TDTest : public ::testing::Test {
protected:
    char DICO_in = 'C';    // Continuous-time system
    int DP_in = 3;         // Degree of the polynomial
    int INFO_out = -999;   // Result from slicot_mc01td call

    // Input polynomial coefficients
    std::vector<double> P_in;

    // Outputs
    int STABLE_out = -1;
    int NZ_out = -1;
    int IWARN_out = -1;

    // Expected results
    int EXPECTED_STABLE = 1; // Stable
    int EXPECTED_NZ = 0;     // No unstable zeros
    int EXPECTED_IWARN = 0;  // No warnings
    int EXPECTED_INFO = 0;   // Successful execution

    void SetUp() override {
        // Initialize polynomial coefficients (P(x) = x^3 + 2x^2 + 3x + 4)
        P_in = {4.0, 3.0, 2.0, 1.0};
    }
};

// --- Test Cases ---

// Test: Continuous-time stability check
TEST_F(MC01TDTest, ContinuousTimeStability) {
    INFO_out = slicot_mc01td(DICO_in, DP_in, P_in.data(), &STABLE_out, &NZ_out, &IWARN_out);

    ASSERT_EQ(INFO_out, EXPECTED_INFO) << "slicot_mc01td call failed with INFO = " << INFO_out;
    EXPECT_EQ(STABLE_out, EXPECTED_STABLE) << "STABLE mismatch";
    EXPECT_EQ(NZ_out, EXPECTED_NZ) << "NZ mismatch";
    EXPECT_EQ(IWARN_out, EXPECTED_IWARN) << "IWARN mismatch";
}

// Test: Discrete-time stability check
TEST_F(MC01TDTest, DiscreteTimeStability) {
    DICO_in = 'D'; // Discrete-time system
    INFO_out = slicot_mc01td(DICO_in, DP_in, P_in.data(), &STABLE_out, &NZ_out, &IWARN_out);

    // Recalculate expected results based on the roots of P(x) = x^3 + 2x^2 + 3x + 4
    EXPECTED_STABLE = 0; // The polynomial is unstable in the discrete-time case
    EXPECTED_NZ = 3;     // All 3 roots are outside the unit circle

    ASSERT_EQ(INFO_out, EXPECTED_INFO) << "slicot_mc01td call failed with INFO = " << INFO_out;
    EXPECT_EQ(STABLE_out, EXPECTED_STABLE) << "STABLE mismatch";
    EXPECT_EQ(NZ_out, EXPECTED_NZ) << "NZ mismatch";
    EXPECT_EQ(IWARN_out, EXPECTED_IWARN) << "IWARN mismatch";
}

// Test: Zero-degree polynomial
TEST_F(MC01TDTest, ZeroDegreePolynomial) {
    DP_in = 0; // Zero-degree polynomial
    P_in = {1.0}; // Constant polynomial P(x) = 1

    INFO_out = slicot_mc01td(DICO_in, DP_in, P_in.data(), &STABLE_out, &NZ_out, &IWARN_out);

    ASSERT_EQ(INFO_out, EXPECTED_INFO) << "slicot_mc01td call failed with INFO = " << INFO_out;
    EXPECT_EQ(STABLE_out, EXPECTED_STABLE) << "STABLE mismatch";
    EXPECT_EQ(NZ_out, EXPECTED_NZ) << "NZ mismatch";
    EXPECT_EQ(IWARN_out, EXPECTED_IWARN) << "IWARN mismatch";
}

// Test: Invalid degree
TEST_F(MC01TDTest, InvalidDegree) {
    DP_in = -1; // Invalid degree
    INFO_out = slicot_mc01td(DICO_in, DP_in, P_in.data(), &STABLE_out, &NZ_out, &IWARN_out);

    EXPECT_EQ(INFO_out, -2) << "Expected INFO = -2 for invalid degree";
}

// Test: Invalid DICO
TEST_F(MC01TDTest, InvalidDICO) {
    DICO_in = 'X'; // Invalid DICO
    INFO_out = slicot_mc01td(DICO_in, DP_in, P_in.data(), &STABLE_out, &NZ_out, &IWARN_out);

    EXPECT_EQ(INFO_out, -1) << "Expected INFO = -1 for invalid DICO";
}
