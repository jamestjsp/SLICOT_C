#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <string>
#include <stdexcept> // For std::runtime_error
#include <algorithm> // For std::max
#include <iostream>  // For std::cerr
#include <filesystem> // Required for path joining (C++17)

// Include headers for both wrappers and utilities
#include "ib01ad.h"
#include "ib01bd.h"
#include "slicot_utils.h"
#include "test_utils.h" // For load_test_data_from_csv
#include "test_config.h" // Include the new test configuration

// Using ib01ad.csv as it contains the U,Y inputs needed to generate R for IB01BD
const std::string DATA_FILE_PATH = TEST_DATA_DIR "ib01ad.csv";


// --- Column-Major Test Fixture for IB01BD ---
class IB01BDTestColMajor : public ::testing::Test {
protected:
    // === Parameters for IB01AD (used in SetUp) ===
    int ad_NOBR = 15;
    int ad_M = 1; // Corresponds to {"U"}
    int ad_L = 1; // Corresponds to {"Y"}
    int ad_NSMP = 1000; // Default/Expected number of samples for IB01AD
    double ad_RCOND = 0.0;
    double ad_TOL = -1.0;
    char ad_METH = 'M'; // MOESP
    char ad_ALG = 'C';  // Cholesky
    // Using JOBD = 'M' as required by IB01BD METH='M', JOB='A'/'B'/'D'
    char ad_JOBD = 'M';
    char ad_BATCH = 'O';
    char ad_CONCT = 'N';
    char ad_CTRL = 'N';
    // Column names for IB01AD loader
    std::vector<std::string> ad_input_columns = {"U"}; // Use actual header name from ib01ad.csv
    std::vector<std::string> ad_output_columns = {"Y"};// Use actual header name from ib01ad.csv

    // === Parameters for IB01BD (based on example) ===
    char bd_METH = 'M'; // MOESP (consistent with ad_METH for this example)
    // **DEBUG: Changed JOB to 'C' and JOBCK to 'N' to test A/C calculation only**
    char bd_JOB = 'C';  // Compute A, C only
    char bd_JOBCK = 'N';// Do not compute Covariances/Kalman Gain
    // NOBR, N, M, L will be set based on IB01AD results/params
    int bd_NSMPL = 0; // Total samples, will be set from ad_NSMP (not used if JOBCK='N')
    double bd_TOL = 0.0; // Use default tolerance for IB01BD rank estimation

    // === Data Vectors ===
    // Inputs for IB01AD
    std::vector<double> U_ad;
    std::vector<double> Y_ad;
    // Outputs from IB01AD (Inputs for IB01BD)
    std::vector<double> R_computed; // Computed by IB01AD (Col-Major)
    std::vector<double> SV_computed; // Computed by IB01AD
    int N_computed = 0; // System order computed by IB01AD
    int LDR_computed = 0; // Leading dimension of R_computed
    // Outputs from IB01BD
    std::vector<double> A_out, C_out; // B, D, Q, RY, S, K not computed with JOB='C', JOBCK='N'
    // std::vector<double> B_out, D_out; // Removed for this debug step
    // std::vector<double> Q_out, RY_out, S_out, K_out; // Removed for this debug step

    // === Expected Results for IB01BD (Column-Major) ===
    int expected_info_bd = 0;
    int expected_iwarn_bd = 0; // Example shows no warning
    // Expected matrices from example output (COLUMN-MAJOR format)
    std::vector<double> A_expected = { 0.8924, -0.0837, 0.0052, 0.0055, 0.3887, 0.6186, 0.1307, 0.0734, 0.1285, -0.6273, 0.6685, -0.2148, 0.1716, -0.4582, -0.6755, 0.4788 };
    std::vector<double> C_expected = { -0.4442, 0.6663, 0.3961, 0.4102 };
    // B, D, Q, Ry, S, K expected values removed for this debug step

    // === Other Variables ===
    double check_tol = 1e-4; // Tolerance for checking matrix elements
    // Leading dimensions (calculated in SetUp)
    int LDA = 1, LDC = 1; // LDB, LDD etc. removed for this debug step

    // SetUp method: Load data, run IB01AD, size outputs for IB01BD
    void SetUp() override {
        // --- 1. Load Data for IB01AD ---
        int samples_loaded = 0;
        ASSERT_EQ(ad_input_columns.size(), ad_M) << "IB01AD Fixture M != input_columns size";
        ASSERT_EQ(ad_output_columns.size(), ad_L) << "IB01AD Fixture L != output_columns size";

        try {
            bool success = load_test_data_from_csv(
                DATA_FILE_PATH, ad_input_columns, ad_output_columns,
                U_ad, Y_ad, samples_loaded
            );
            ASSERT_TRUE(success) << "CSV loading reported failure for " << DATA_FILE_PATH;
            ASSERT_GT(samples_loaded, 0) << "No data samples loaded from CSV: " << DATA_FILE_PATH;
            ad_NSMP = samples_loaded;
            bd_NSMPL = samples_loaded;

        } catch (const std::runtime_error& e) {
            FAIL() << "CSV data loading failed: " << e.what();
        } catch (...) {
            FAIL() << "Caught unknown exception during CSV data loading.";
        }

        // --- 2. Run IB01AD to get R and N ---
        int n_ad_out = 0;
        int iwarn_ad = 0;
        int info_ad = -999;
        LDR_computed = std::max(1, 2 * (ad_M + ad_L) * ad_NOBR);
        if (ad_METH == 'M' && ad_JOBD == 'M') {
             LDR_computed = std::max(LDR_computed, 3 * ad_M * ad_NOBR);
        }
        R_computed.resize((size_t)LDR_computed * 2 * (ad_M + ad_L) * ad_NOBR);
        SV_computed.resize((size_t)ad_L * ad_NOBR);
        info_ad = slicot_ib01ad(ad_METH, ad_ALG, ad_JOBD, ad_BATCH, ad_CONCT, ad_CTRL,
                                ad_NOBR, ad_M, ad_L, ad_NSMP,
                                (ad_M > 0 ? U_ad.data() : nullptr), ad_NSMP,
                                (ad_L > 0 ? Y_ad.data() : nullptr), ad_NSMP,
                                &n_ad_out, R_computed.data(), LDR_computed, SV_computed.data(),
                                ad_RCOND, ad_TOL, &iwarn_ad, 0 /* col-major */);
        ASSERT_EQ(info_ad, 0) << "slicot_ib01ad failed in SetUp with info = " << info_ad;
        N_computed = n_ad_out;
        ASSERT_GT(N_computed, 0) << "IB01AD returned non-positive system order N.";
        ASSERT_EQ(N_computed, 4) << "IB01AD computed N (" << N_computed << ") differs from expected (4).";

        // --- 3. Size Output Arrays for IB01BD (A and C only) ---
        int n = N_computed;
        int l = ad_L;
        if (bd_JOB == 'A' || bd_JOB == 'C') { // Condition is true for bd_JOB='C'
             A_out.resize((size_t)n * n);
             C_out.resize((size_t)l * n);
        }

        // --- 4. Calculate Leading Dimensions for IB01BD (A and C only) ---
        LDA = (bd_JOB == 'A' || bd_JOB == 'C') ? std::max(1, n) : 1;
        LDC = (bd_JOB == 'A' || bd_JOB == 'C') ? std::max(1, l) : 1;
        // LDB, LDD etc. not needed for this test configuration
    }
};

// --- Row-Major Test Fixture for IB01BD ---
class IB01BDTestRowMajor : public IB01BDTestColMajor {
protected:
    // Row-major vectors for outputs (for verification)
    std::vector<double> A_expected_rm, C_expected_rm;
    // B, D, Q etc expected removed for this debug step

    // Row-major leading dimensions (cols)
    int LDA_rm = 1, LDC_rm = 1;
    // LDB_rm, LDD_rm etc. removed for this debug step

    void SetUp() override {
        // Run base class SetUp first to load data and run IB01AD (with JOBD='M')
        IB01BDTestColMajor::SetUp();

        int n = N_computed;
        int l = ad_L;

        // --- Convert Expected Results to Row-Major for Comparison ---
        if (bd_JOB == 'A' || bd_JOB == 'C') { // True for bd_JOB='C'
            A_expected_rm.resize(A_expected.size()); if (!A_expected.empty()) slicot_transpose_to_c(A_expected.data(), A_expected_rm.data(), n, n, sizeof(double));
            C_expected_rm.resize(C_expected.size()); if (!C_expected.empty()) slicot_transpose_to_c(C_expected.data(), C_expected_rm.data(), l, n, sizeof(double));
        }

        // --- Calculate Row-Major Leading Dimensions (cols) ---
        LDA_rm = (bd_JOB == 'A' || bd_JOB == 'C') ? std::max(1, n) : 1;
        LDC_rm = (bd_JOB == 'A' || bd_JOB == 'C') ? std::max(1, n) : 1; // n cols for C
    }
};


// --- Test Cases ---

// Test: Documentation Example (Column-Major) - Modified for A/C only
TEST_F(IB01BDTestColMajor, DocExample_AC_Only) { // Renamed test slightly
    int iwarn_bd = -1;
    int n = N_computed; // Use N computed by IB01AD
    int m = ad_M;
    int l = ad_L;
    int info_result = -999; // Initialize local info

    // Call slicot_ib01bd with JOB='C', JOBCK='N'
    // Pass NULL for pointers not needed (B, D, Q, Ry, S, K) and LD=1
    info_result = slicot_ib01bd(bd_METH, bd_JOB, bd_JOBCK, ad_NOBR, n, m, l,
                                bd_NSMPL, R_computed.data(), LDR_computed,
                                A_out.data(), LDA, C_out.data(), LDC,
                                nullptr, 1, nullptr, 1, // B, D
                                nullptr, 1, nullptr, 1, nullptr, 1, // Q, Ry, S
                                nullptr, 1, // K
                                bd_TOL, &iwarn_bd, 0 /* row_major = false */);

    // Verify return code
    ASSERT_EQ(info_result, expected_info_bd) << "slicot_ib01bd (A/C only) failed with info = " << info_result;
    EXPECT_EQ(iwarn_bd, expected_iwarn_bd) << "slicot_ib01bd (A/C only) returned warning = " << iwarn_bd;

    // Verify output matrices A and C (Column-Major comparison)
    if (bd_JOB == 'A' || bd_JOB == 'C') { // True for bd_JOB='C'
        ASSERT_EQ(A_out.size(), A_expected.size()); for (size_t i = 0; i < A_expected.size(); ++i) EXPECT_NEAR(A_out[i], A_expected[i], check_tol) << "A mismatch at index " << i;
        ASSERT_EQ(C_out.size(), C_expected.size()); for (size_t i = 0; i < C_expected.size(); ++i) EXPECT_NEAR(C_out[i], C_expected[i], check_tol) << "C mismatch at index " << i;
    }
    // No need to check B, D, Q, Ry, S, K
}

// Test: Documentation Example (Row-Major) - Modified for A/C only
TEST_F(IB01BDTestRowMajor, DocExample_AC_Only) { // Renamed test slightly
    int iwarn_bd = -1;
    int n = N_computed; // Use N computed by IB01AD
    int m = ad_M;
    int l = ad_L;
    int info_result = -999; // Initialize local info

    // Call slicot_ib01bd with JOB='C', JOBCK='N'
    info_result = slicot_ib01bd(bd_METH, bd_JOB, bd_JOBCK, ad_NOBR, n, m, l,
                                bd_NSMPL, R_computed.data(), LDR_computed, // R is always col-major
                                A_out.data(), LDA_rm, C_out.data(), LDC_rm, // Pass output buffers and row-major LDs
                                nullptr, 1, nullptr, 1, // B, D
                                nullptr, 1, nullptr, 1, nullptr, 1, // Q, Ry, S
                                nullptr, 1, // K
                                bd_TOL, &iwarn_bd, 1 /* row_major = true */);

    // Verify return code
    ASSERT_EQ(info_result, expected_info_bd) << "slicot_ib01bd (A/C only) failed with info = " << info_result;
    EXPECT_EQ(iwarn_bd, expected_iwarn_bd) << "slicot_ib01bd (A/C only) returned warning = " << iwarn_bd;

    // Verify output matrices A and C (Row-Major comparison)
    if (bd_JOB == 'A' || bd_JOB == 'C') { // True for bd_JOB='C'
        ASSERT_EQ(A_out.size(), A_expected_rm.size()); for (size_t i = 0; i < A_expected_rm.size(); ++i) EXPECT_NEAR(A_out[i], A_expected_rm[i], check_tol) << "A row-major mismatch at index " << i;
        ASSERT_EQ(C_out.size(), C_expected_rm.size()); for (size_t i = 0; i < C_expected_rm.size(); ++i) EXPECT_NEAR(C_out[i], C_expected_rm[i], check_tol) << "C row-major mismatch at index " << i;
    }
     // No need to check B, D, Q, Ry, S, K
}


// Test: Parameter Validation
TEST_F(IB01BDTestColMajor, ParameterValidation) {
    // Use dummy inputs where possible, focus on validating wrapper checks
    int n = N_computed > 0 ? N_computed : 4; // Use computed N or default
    int m = ad_M;
    int l = ad_L;
    int nobr = ad_NOBR;
    int nsmpl = bd_NSMPL > 0 ? bd_NSMPL : 1000; // Use actual or dummy
    int ldr = LDR_computed > 0 ? LDR_computed : 2*(m+l)*nobr; // Use actual or min
    std::vector<double> dummy_r(ldr * 2*(m+l)*nobr, 0.0); // Need valid R size
    std::vector<double> dummy_a(1), dummy_c(1), dummy_b(1), dummy_d(1);
    std::vector<double> dummy_q(1), dummy_ry(1), dummy_s(1), dummy_k(1);
    int iwarn = 0;
    int info;

    // Test invalid METH
    info = slicot_ib01bd('X', bd_JOB, bd_JOBCK, nobr, n, m, l, nsmpl, dummy_r.data(), ldr,
                         dummy_a.data(), 1, dummy_c.data(), 1, dummy_b.data(), 1, dummy_d.data(), 1,
                         dummy_q.data(), 1, dummy_ry.data(), 1, dummy_s.data(), 1, dummy_k.data(), 1,
                         bd_TOL, &iwarn, 0);
    EXPECT_EQ(info, -1);

    // Test invalid JOB
    info = slicot_ib01bd(bd_METH, 'X', bd_JOBCK, nobr, n, m, l, nsmpl, dummy_r.data(), ldr,
                         dummy_a.data(), 1, dummy_c.data(), 1, dummy_b.data(), 1, dummy_d.data(), 1,
                         dummy_q.data(), 1, dummy_ry.data(), 1, dummy_s.data(), 1, dummy_k.data(), 1,
                         bd_TOL, &iwarn, 0);
    EXPECT_EQ(info, -2);

    // Test invalid JOBCK
     info = slicot_ib01bd(bd_METH, bd_JOB, 'X', nobr, n, m, l, nsmpl, dummy_r.data(), ldr,
                         dummy_a.data(), 1, dummy_c.data(), 1, dummy_b.data(), 1, dummy_d.data(), 1,
                         dummy_q.data(), 1, dummy_ry.data(), 1, dummy_s.data(), 1, dummy_k.data(), 1,
                         bd_TOL, &iwarn, 0);
    EXPECT_EQ(info, -3);

    // Test invalid NOBR (<= 1)
     info = slicot_ib01bd(bd_METH, bd_JOB, bd_JOBCK, 1, n, m, l, nsmpl, dummy_r.data(), ldr,
                         dummy_a.data(), 1, dummy_c.data(), 1, dummy_b.data(), 1, dummy_d.data(), 1,
                         dummy_q.data(), 1, dummy_ry.data(), 1, dummy_s.data(), 1, dummy_k.data(), 1,
                         bd_TOL, &iwarn, 0);
    EXPECT_EQ(info, -4);

    // Test invalid N (<= 0)
     info = slicot_ib01bd(bd_METH, bd_JOB, bd_JOBCK, nobr, 0, m, l, nsmpl, dummy_r.data(), ldr,
                         dummy_a.data(), 1, dummy_c.data(), 1, dummy_b.data(), 1, dummy_d.data(), 1,
                         dummy_q.data(), 1, dummy_ry.data(), 1, dummy_s.data(), 1, dummy_k.data(), 1,
                         bd_TOL, &iwarn, 0);
    EXPECT_EQ(info, -5);

     // Test invalid N (>= NOBR)
     info = slicot_ib01bd(bd_METH, bd_JOB, bd_JOBCK, nobr, nobr, m, l, nsmpl, dummy_r.data(), ldr,
                         dummy_a.data(), 1, dummy_c.data(), 1, dummy_b.data(), 1, dummy_d.data(), 1,
                         dummy_q.data(), 1, dummy_ry.data(), 1, dummy_s.data(), 1, dummy_k.data(), 1,
                         bd_TOL, &iwarn, 0);
    EXPECT_EQ(info, -5);

    // Test invalid LDR
     info = slicot_ib01bd(bd_METH, bd_JOB, bd_JOBCK, nobr, n, m, l, nsmpl, dummy_r.data(), 1, // LDR too small
                         dummy_a.data(), 1, dummy_c.data(), 1, dummy_b.data(), 1, dummy_d.data(), 1,
                         dummy_q.data(), 1, dummy_ry.data(), 1, dummy_s.data(), 1, dummy_k.data(), 1,
                         bd_TOL, &iwarn, 0);
    EXPECT_EQ(info, -10);

    // Add more checks for other LDs and NULL pointers based on JOB/JOBCK/METH...
}
