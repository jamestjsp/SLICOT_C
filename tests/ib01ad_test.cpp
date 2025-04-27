#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm> // For std::max
#include <iostream>  // For std::cerr
#include <string>    // For std::string

#include "ib01ad.h"       // Include the updated wrapper header
#include "slicot_utils.h" // For transpose functions
#include "test_utils.h"   // For data loading utilities (contains the updated loader)
#include "test_config.h"  // Include the new test configuration

// Path to test data file (ensure this path is correct relative to execution)
// ASSUMPTION: This CSV file has a header row with columns named "U" and "Y".
const std::string DATA_FILE_PATH = TEST_DATA_DIR "ib01ad.csv";

// --- Column-Major Test Fixture ---
class IB01ADTestColMajor : public ::testing::Test {
protected:
    // Parameters from IB01AD.dat example (NSMP might be updated by loader)
    int NOBR = 15;
    int M = 1; // Corresponds to {"U"}
    int L = 1; // Corresponds to {"Y"}
    int NSMP = 1000; // Default/Expected number of samples
    double RCOND = 0.0;
    double TOL = -1.0;
    char METH = 'M';
    char ALG = 'C';
    char JOBD = 'N';
    char BATCH = 'O';
    char CONCT = 'N';
    char CTRL = 'N';

    // Column names to load from CSV - CORRECTED
    std::vector<std::string> input_columns = {"U"}; // Use actual header name
    std::vector<std::string> output_columns = {"Y"};// Use actual header name

    double check_tol = 1e-4; // Tolerance for checking singular values

    // Input data (U and Y) - Column Major - Loaded from CSV
    std::vector<double> U;
    std::vector<double> Y;

    // Expected results from IB01AD.res
    int expected_N = 4;
    std::vector<double> SV_expected = {
        69.8841, 14.9963, 3.6675, 1.9677, 0.3000, 0.2078, 0.1651, 0.1373,
        0.1133,  0.1059,  0.0856, 0.0784, 0.0733, 0.0678, 0.0571
    };

    // Output arrays
    std::vector<double> R;
    std::vector<double> SV;

    // Result variables
    int n_result = -1;
    int iwarn_result = -1;

    // Calculated leading dimensions (based on potentially updated NSMP)
    int ldu = 1;
    int ldy = 1;
    int ldr = 1;

    void SetUp() override {
        // --- Load data using the updated CSV loader ---
        int samples_loaded = 0;
        ASSERT_EQ(input_columns.size(), M) << "Test parameter M doesn't match number of input_columns.";
        ASSERT_EQ(output_columns.size(), L) << "Test parameter L doesn't match number of output_columns.";

        try {
            // Ensure test_utils.h declares the function with the correct signature
            bool success = load_test_data_from_csv(
                DATA_FILE_PATH,
                input_columns,
                output_columns,
                U, // Output U vector
                Y, // Output Y vector
                samples_loaded // Output number of samples loaded
            );
            ASSERT_TRUE(success) << "CSV loading reported failure."; // Should throw on error, but check bool too
            ASSERT_GT(samples_loaded, 0) << "No data samples loaded from CSV: " << DATA_FILE_PATH;

            // Update NSMP based on actual data loaded
            if (samples_loaded != NSMP) {
                 std::cout << "Note: Loaded " << samples_loaded << " samples, updating NSMP from " << NSMP << "." << std::endl;
                 NSMP = samples_loaded;
            }

        } catch (const std::runtime_error& e) {
            // Fail the test immediately if loading fails
            FAIL() << "CSV data loading failed: " << e.what();
        } catch (...) {
            FAIL() << "Caught unknown exception during CSV data loading.";
        }

        // --- Array Sizing (using potentially updated NSMP) ---
        // Leading dimensions (Column Major)
        ldu = (M > 0) ? NSMP : 1;
        ldy = (L > 0) ? NSMP : 1; // L must be > 0 per docs
        int min_ldr_f = (METH == 'M' && JOBD == 'N') ? 2*(M+L)*NOBR :
                       ((METH == 'M' && JOBD == 'M') ? std::max(2*(M+L)*NOBR, 3*M*NOBR) : 2*(M+L)*NOBR);
        ldr = std::max(1, min_ldr_f);

        // R array size: ldr rows, 2*(M + L)*NOBR columns
        R.resize((size_t)ldr * 2 * (M + L) * NOBR);

        // SV array size: L*NOBR
        SV.resize((size_t)L * NOBR);
    }
};

// --- Row-Major Test Fixture ---
class IB01ADTestRowMajor : public IB01ADTestColMajor { // Inherit parameters and expected results
public:
    // Input data in row-major format
    std::vector<double> U_rm;
    std::vector<double> Y_rm;

    // Row-major leading dimensions
    int ldu_rm = 1;
    int ldy_rm = 1;
    int ldr_rm = 1;

    void SetUp() override {
        // --- Load column-major data first using the updated CSV loader ---
        int samples_loaded = 0;
        ASSERT_EQ(input_columns.size(), M) << "Test parameter M doesn't match number of input_columns.";
        ASSERT_EQ(output_columns.size(), L) << "Test parameter L doesn't match number of output_columns.";

        std::vector<double> U_col, Y_col; // Temporary column-major storage
        try {
             bool success = load_test_data_from_csv(
                DATA_FILE_PATH,
                input_columns, // Use corrected names {"U"}, {"Y"}
                output_columns,
                U_col, // Load into temporary col-major U
                Y_col, // Load into temporary col-major Y
                samples_loaded
            );
            ASSERT_TRUE(success) << "CSV loading reported failure.";
            ASSERT_GT(samples_loaded, 0) << "No data samples loaded from CSV: " << DATA_FILE_PATH;

            // Update NSMP based on actual data loaded
             if (samples_loaded != NSMP) {
                 std::cout << "Note: Loaded " << samples_loaded << " samples, updating NSMP from " << NSMP << "." << std::endl;
                 NSMP = samples_loaded; // Update base class NSMP as well
            }

        } catch (const std::runtime_error& e) {
            FAIL() << "CSV data loading failed: " << e.what();
        } catch (...) {
            FAIL() << "Caught unknown exception during CSV data loading.";
        }

        // --- Convert to row-major format ---
        U_rm.resize(U_col.size());
        Y_rm.resize(Y_col.size());
        // Use U_col and Y_col for transpose source
        if (M > 0 && U_col.size() > 0) slicot_transpose_to_c(U_col.data(), U_rm.data(), NSMP, M, sizeof(double));
        if (L > 0 && Y_col.size() > 0) slicot_transpose_to_c(Y_col.data(), Y_rm.data(), NSMP, L, sizeof(double));

        // --- Array Sizing (using potentially updated NSMP) ---
        // Leading dimensions (Row Major)
        ldu_rm = (M > 0) ? M : 1;
        ldy_rm = (L > 0) ? L : 1;
        ldr_rm = 2 * (M + L) * NOBR;

        // R array size (based on Fortran requirements)
        int min_ldr_f = (METH == 'M' && JOBD == 'N') ? 2*(M+L)*NOBR :
                       ((METH == 'M' && JOBD == 'M') ? std::max(2*(M+L)*NOBR, 3*M*NOBR) : 2*(M+L)*NOBR);
        int r_rows_f = std::max(1, min_ldr_f);
        int r_cols_f = 2 * (M + L) * NOBR;
        R.resize((size_t)r_rows_f * r_cols_f);

        // SV array size
        SV.resize((size_t)L * NOBR);
    }
};


// --- Test Cases ---

// Test using data from IB01AD.dat/.res (Column-Major)
TEST_F(IB01ADTestColMajor, DocExample) {
    int info;
    int n = 0; // Auto-detect order

    // Run Actual Computation - Workspace handled internally by wrapper
    info = slicot_ib01ad(
        METH, ALG, JOBD, BATCH, CONCT, CTRL,
        NOBR, M, L, NSMP,
        (M > 0 ? U.data() : nullptr), ldu, (L > 0 ? Y.data() : nullptr), ldy,
        &n, R.data(), ldr, SV.data(), RCOND,
        TOL, &iwarn_result, 0  // Column-major
    );

    // Save result
    n_result = n;

    // Verify results
    ASSERT_EQ(info, 0) << "slicot_ib01ad failed with info = " << info;
    EXPECT_EQ(iwarn_result, 0) << "slicot_ib01ad returned warning = " << iwarn_result;
    EXPECT_EQ(n_result, expected_N) << "System order N does not match expected value.";

    // Check singular values
    ASSERT_GE(SV.size(), SV_expected.size()) << "Output SV vector too small.";
    for (size_t i = 0; i < SV_expected.size(); ++i) { // Only check the expected number of values
        EXPECT_NEAR(SV[i], SV_expected[i], check_tol)
            << "Mismatch in singular value at index " << i;
    }
}

// Test using data from IB01AD.dat/.res (Row-Major)
TEST_F(IB01ADTestRowMajor, DocExample) {
    int info;
    int n = 0; // Auto-detect order

    // Run Actual Computation - Workspace handled internally by wrapper
    info = slicot_ib01ad(
        METH, ALG, JOBD, BATCH, CONCT, CTRL,
        NOBR, M, L, NSMP,
        (M > 0 ? U_rm.data() : nullptr), ldu_rm, (L > 0 ? Y_rm.data() : nullptr), ldy_rm,
        &n, R.data(), ldr_rm, SV.data(), RCOND, // Use ldr_rm for row-major call
        TOL, &iwarn_result, 1  // Row-major
    );

    // Save result
    n_result = n;

    // Verify results
    ASSERT_EQ(info, 0) << "slicot_ib01ad failed with info = " << info;
    EXPECT_EQ(iwarn_result, 0) << "slicot_ib01ad returned warning = " << iwarn_result;
    EXPECT_EQ(n_result, expected_N) << "System order N does not match expected value.";

    // Check singular values
    ASSERT_GE(SV.size(), SV_expected.size()) << "Output SV vector too small.";
     for (size_t i = 0; i < SV_expected.size(); ++i) { // Only check the expected number of values
        EXPECT_NEAR(SV[i], SV_expected[i], check_tol)
            << "Mismatch in singular value at index " << i;
    }
}


// Test parameter validation (using Column-Major fixture - doesn't rely on loaded data)
TEST_F(IB01ADTestColMajor, ParameterValidation) {
    // This test doesn't need real data, just checks parameter validation in the wrapper
    std::vector<double> dummy_u(1);
    std::vector<double> dummy_y(1);
    std::vector<double> dummy_r(1);
    std::vector<double> dummy_sv(1);
    int n_out = 0;
    int iwarn = 0;
    int info;
    int dummy_nsmp = 100;
    int dummy_ld = 100;
    int dummy_ldr = std::max(1, (METH == 'M' && JOBD == 'N') ? 2*(M+L)*NOBR : ((METH == 'M' && JOBD == 'M') ? std::max(2*(M+L)*NOBR, 3*M*NOBR) : 2*(M+L)*NOBR));

    // Test invalid METH parameter
    info = slicot_ib01ad('X', ALG, JOBD, BATCH, CONCT, CTRL, NOBR, M, L, dummy_nsmp,
                         dummy_u.data(), dummy_ld, dummy_y.data(), dummy_ld,
                         &n_out, dummy_r.data(), dummy_ldr, dummy_sv.data(), RCOND, TOL,
                         &iwarn, 0);
    EXPECT_EQ(info, -1) << "Invalid METH parameter not detected correctly (expected -1)";

    // Test invalid ALG parameter
    info = slicot_ib01ad(METH, 'X', JOBD, BATCH, CONCT, CTRL, NOBR, M, L, dummy_nsmp,
                         dummy_u.data(), dummy_ld, dummy_y.data(), dummy_ld,
                         &n_out, dummy_r.data(), dummy_ldr, dummy_sv.data(), RCOND, TOL,
                         &iwarn, 0);
    EXPECT_EQ(info, -2) << "Invalid ALG parameter not detected correctly (expected -2)";

    // Test invalid NOBR parameter (<= 0)
    info = slicot_ib01ad(METH, ALG, JOBD, BATCH, CONCT, CTRL, 0, M, L, dummy_nsmp,
                         dummy_u.data(), dummy_ld, dummy_y.data(), dummy_ld,
                         &n_out, dummy_r.data(), dummy_ldr, dummy_sv.data(), RCOND, TOL,
                         &iwarn, 0);
    EXPECT_EQ(info, -7) << "Invalid NOBR parameter not detected correctly (expected -7)";

    // Test invalid L parameter (<= 0)
    info = slicot_ib01ad(METH, ALG, JOBD, BATCH, CONCT, CTRL, NOBR, M, 0, dummy_nsmp,
                         dummy_u.data(), dummy_ld, nullptr, 1,
                         &n_out, dummy_r.data(), dummy_ldr, dummy_sv.data(), RCOND, TOL,
                         &iwarn, 0);
    EXPECT_EQ(info, -9) << "Invalid L parameter not detected correctly (expected -9)";

    // Test invalid NSMP parameter for non-sequential processing (BATCH='O')
    int nsmp_too_small = 2 * NOBR - 1;
    info = slicot_ib01ad(METH, ALG, JOBD, 'O', CONCT, CTRL, NOBR, M, L, nsmp_too_small,
                         dummy_u.data(), std::max(1,nsmp_too_small),
                         dummy_y.data(), std::max(1,nsmp_too_small),
                         &n_out, dummy_r.data(), dummy_ldr, dummy_sv.data(), RCOND, TOL,
                         &iwarn, 0);
    EXPECT_EQ(info, -10) << "Invalid NSMP parameter not detected correctly (expected -10)";
}
