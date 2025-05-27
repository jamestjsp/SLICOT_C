#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm> // For std::max
#include <numeric>   // For std::iota
#include <stdexcept> // For std::runtime_error
#include <iostream>  // For std::cout in tests

#include "ib03ad.h"       // Include the wrapper header
#include "slicot_utils.h" // For transpose functions, MAX/MIN
#include "test_utils.h"   // For load_test_data_from_csv
#include "test_config.h"  // For TEST_DATA_DIR

// --- Base Test Fixture for IB03AD ---
class IB03ADTest : public ::testing::Test {
protected:
    // Parameters matching IB03AD.html example
    char INIT = 'B';
    char ALG = 'D';
    char STOR = 'F';

    int NOBR = 10;
    int M_val = 1;
    int L_val = 1;
    int NSMP_val = 1024; // Will be updated by CSV loader if successful
    int N_val = 4;
    int NN_val = 12;

    int ITMAX1 = 500;
    int ITMAX2 = 1000;
    int NPRINT = 0;

    double TOL1 = 1e-5;
    double TOL2 = 1e-5;

    // Input data
    std::vector<double> U_data;
    std::vector<double> Y_data;
    std::vector<double> X_params; // Parameter vector

    int LX_val = 0; // Length of X_params, to be calculated

    // Results
    int info_result = -999;
    int iwarn_result = -999;
    int n_returned = -1; // To store N returned by the routine
    int lx_returned = -1; // To store LX returned by the routine
    std::vector<int> iwork_summary;
    std::vector<double> dwork_summary;


    // Leading dimensions (calculated in SetUp)
    int LDU_c = 1, LDY_c = 1; // For column-major C arrays
    int LDU_r = 1, LDY_r = 1; // For row-major C arrays

    // CSV configuration
    std::string csv_filename = "ib03ad.csv";
    std::vector<std::string> input_columns = {"U1"}; // Assuming M=1, header "U1"
    std::vector<std::string> output_columns = {"Y1"};// Assuming L=1, header "Y1"


    void SetUpBase(bool use_csv = true) { // Added flag to control CSV loading for specialized tests
        if (use_csv) {
            std::string full_csv_path = std::string(TEST_DATA_DIR) + "/" + csv_filename;
            int samples_loaded = 0;
            bool csv_load_success = false;
            try {
                 csv_load_success = load_test_data_from_csv(
                    full_csv_path, input_columns, output_columns,
                    U_data, Y_data, samples_loaded
                );
            } catch (const std::exception& e) {
                 FAIL() << "CSV data loading failed with exception: " << e.what() << " for file " << full_csv_path;
            }

            if (csv_load_success && samples_loaded > 0) {
                NSMP_val = samples_loaded;
                 // Verify that NSMP_val matches the example if CSV is used for the main example test
                if (csv_filename == "ib03ad.csv") { // Check if it's the specific example data
                    // The example data has 1024 samples.
                    // The load_test_data_from_csv might read fewer if the file is truncated or malformed.
                    // For this specific example, we expect 1024.
                    // Allow some flexibility if the CSV reader has issues with the very long file.
                    // However, for a correct run against example results, NSMP should be 1024.
                    // For now, we just update NSMP_val. If it's not 1024, numerical results might differ.
                }
            } else {
                // Fallback to smaller synthetic data if CSV loading fails or provides no samples,
                // or if use_csv is false.
                NSMP_val = 10; // Default small NSMP for non-CSV or failed CSV
                U_data.resize((size_t)NSMP_val * M_val);
                for (size_t i = 0; i < U_data.size(); ++i) {
                    U_data[i] = sin((double)i * 0.1) * 0.5;
                }
                Y_data.resize((size_t)NSMP_val * L_val);
                double current_x_state = 0.0;
                for (int t = 0; t < NSMP_val; ++t) {
                    double u_t = (M_val > 0) ? U_data[t * M_val] : 0.0;
                    double z_t = current_x_state + ((M_val > 0) ? (0.1 * u_t) : 0.0) ;
                    Y_data[t * L_val] = z_t + 0.05 * z_t * z_t + (rand() % 100 / 2000.0 - 0.025);
                    current_x_state = 0.8 * current_x_state + ((M_val > 0) ? (0.5 * u_t) : 0.0);
                }
                if (use_csv) { // Only fail if CSV was intended but failed
                     FAIL() << "CSV data loading failed or returned no samples from " << full_csv_path
                            << ". Ensure CSV file is correct and headers match: U1, Y1 (for M=1, L=1).";
                }
            }
        } else { // Not using CSV, setup for minimal tests
             NSMP_val = 10; // Default for specialized tests
             // Ensure M_val and L_val are set correctly by the specialized test's SetUp
             U_data.assign((size_t)NSMP_val * M_val, 0.1);
             Y_data.assign((size_t)NSMP_val * L_val, 0.1);
        }


        // Calculate LX based on N_val, M_val, L_val, NN_val
        // LX = (NN*(L+2) + 1)*L + N*(L+M+1) + L*M
        int current_L = L_val > 0 ? L_val : 0; // Use 0 if L_val is 0 for LX calc
        int current_M = M_val > 0 ? M_val : 0;
        int current_N = N_val > 0 ? N_val : 0;
        int current_NN = NN_val > 0 ? NN_val : 0;

        int bsn = current_NN * (current_L + 2) + 1;
        int num_nonlinear_params = bsn * current_L;
        int num_linear_params = current_N * (current_L + current_M + 1) + current_L * current_M;
        LX_val = num_nonlinear_params + num_linear_params;
        if (LX_val < 0) LX_val = 0;

        X_params.assign(LX_val > 0 ? LX_val : 1, 0.0); // Ensure min size 1 if LX_val is 0
                                                  // For INIT='B', X is output only. For other INIT, X might be input.

        iwork_summary.resize(3);
        dwork_summary.resize(10 + MAX_EXPECTED_RCONDS_IN_DWORK);

        // Calculate Column-Major C Leading Dimensions
        LDU_c = (M_val > 0) ? MAX(1, NSMP_val) : 1;
        LDY_c = (L_val > 0) ? MAX(1, NSMP_val) : 1;

        // Calculate Row-Major C Leading Dimensions (number of columns)
        LDU_r = (M_val > 0) ? M_val : 1;
        LDY_r = (L_val > 0) ? L_val : 1;
    }
    // Estimate for sizing dwork_summary to hold potential RCOND values
    static const int MAX_EXPECTED_RCONDS_IN_DWORK = 30;
};

// --- Column-Major Test Fixture ---
class IB03ADTestColMajor : public IB03ADTest {
protected:
    void SetUp() override {
        // Parameters for the main CSV example
        INIT = 'B'; ALG = 'D'; STOR = 'F';
        NOBR = 10; M_val = 1; L_val = 1; NSMP_val = 1024; // NSMP_val will be updated by CSV
        N_val = 4; NN_val = 12;
        ITMAX1 = 500; ITMAX2 = 1000; NPRINT = 0;
        TOL1 = 1e-5; TOL2 = 1e-5;
        csv_filename = "ib03ad.csv";
        input_columns = {"U1"};
        output_columns = {"Y1"};
        SetUpBase(true); // Use CSV for these main tests
    }
};

// --- Row-Major Test Fixture ---
class IB03ADTestRowMajor : public IB03ADTest {
protected:
    std::vector<double> U_data_rm;
    std::vector<double> Y_data_rm;

    void SetUp() override {
        // Parameters for the main CSV example
        INIT = 'B'; ALG = 'D'; STOR = 'F';
        NOBR = 10; M_val = 1; L_val = 1; NSMP_val = 1024; // NSMP_val will be updated by CSV
        N_val = 4; NN_val = 12;
        ITMAX1 = 500; ITMAX2 = 1000; NPRINT = 0;
        TOL1 = 1e-5; TOL2 = 1e-5;
        csv_filename = "ib03ad.csv";
        input_columns = {"U1"};
        output_columns = {"Y1"};
        SetUpBase(true); // Use CSV for these main tests

        if (M_val > 0 && NSMP_val > 0 && !U_data.empty()) {
            U_data_rm.resize(U_data.size());
            slicot_transpose_to_c_with_ld(U_data.data(), U_data_rm.data(), NSMP_val, M_val, LDU_c, LDU_r, sizeof(double));
        }
        if (L_val > 0 && NSMP_val > 0 && !Y_data.empty()) {
            Y_data_rm.resize(Y_data.size());
            slicot_transpose_to_c_with_ld(Y_data.data(), Y_data_rm.data(), NSMP_val, L_val, LDY_c, LDY_r, sizeof(double));
        }
    }
};

// --- Test Cases ---

TEST_F(IB03ADTestColMajor, CSVExample_InitB_AlgD_StorF) {
    n_returned = N_val;
    lx_returned = LX_val;

    // Ensure U_data and Y_data are populated (SetUpBase should handle this)
    ASSERT_FALSE(U_data.empty()) << "U_data is empty before calling slicot_ib03ad.";
    ASSERT_FALSE(Y_data.empty()) << "Y_data is empty before calling slicot_ib03ad.";
    ASSERT_GT(LX_val, 0) << "LX_val is not positive before calling slicot_ib03ad.";


    info_result = slicot_ib03ad(INIT, ALG, STOR, NOBR, M_val, L_val, NSMP_val,
                                &n_returned, NN_val, ITMAX1, ITMAX2, NPRINT,
                                U_data.data(), LDU_c,
                                Y_data.data(), LDY_c,
                                X_params.data(), &lx_returned,
                                TOL1, TOL2,
                                iwork_summary.data(), dwork_summary.data(),
                                &iwarn_result, 0 /* row_major */);

    EXPECT_EQ(info_result, 0) << "SLICOT routine IB03AD failed with CSV data.";
    if (info_result == 0) {
        EXPECT_EQ(n_returned, N_val);
        EXPECT_EQ(lx_returned, LX_val);
        EXPECT_EQ(iwarn_result, 0);
        std::cout << "CSVExample_InitB_AlgD_StorF (ColMajor): INFO = " << info_result
                  << ", IWARN = " << iwarn_result
                  << ", Final Error Norm (DWORK[1]): " << (dwork_summary.size() > 1 ? dwork_summary[1] : -1.0)
                  << ", Iterations (DWORK[2]): " << (dwork_summary.size() > 2 ? dwork_summary[2] : -1.0)
                  << std::endl;
    }
}

TEST_F(IB03ADTestRowMajor, CSVExample_InitB_AlgD_StorF) {
    n_returned = N_val;
    lx_returned = LX_val;

    ASSERT_FALSE(U_data_rm.empty()) << "U_data_rm is empty before calling slicot_ib03ad.";
    ASSERT_FALSE(Y_data_rm.empty()) << "Y_data_rm is empty before calling slicot_ib03ad.";
    ASSERT_GT(LX_val, 0) << "LX_val is not positive before calling slicot_ib03ad.";

    info_result = slicot_ib03ad(INIT, ALG, STOR, NOBR, M_val, L_val, NSMP_val,
                                &n_returned, NN_val, ITMAX1, ITMAX2, NPRINT,
                                U_data_rm.data(), LDU_r,
                                Y_data_rm.data(), LDY_r,
                                X_params.data(), &lx_returned,
                                TOL1, TOL2,
                                iwork_summary.data(), dwork_summary.data(),
                                &iwarn_result, 1 /* row_major */);

    EXPECT_EQ(info_result, 0) << "SLICOT routine IB03AD failed with row-major CSV data.";
    if (info_result == 0) {
        EXPECT_EQ(n_returned, N_val);
        EXPECT_EQ(lx_returned, LX_val);
        EXPECT_EQ(iwarn_result, 0);
         std::cout << "CSVExample_InitB_AlgD_StorF (RowMajor): INFO = " << info_result
                  << ", IWARN = " << iwarn_result
                  << ", Final Error Norm (DWORK[1]): " << (dwork_summary.size() > 1 ? dwork_summary[1] : -1.0)
                  << ", Iterations (DWORK[2]): " << (dwork_summary.size() > 2 ? dwork_summary[2] : -1.0)
                  << std::endl;
    }
}

class IB03ADTestParamValidation : public IB03ADTest {
protected:
    void SetUp() override {
        INIT = 'N'; 
        N_val = 1; M_val = 1; L_val = 1; NSMP_val = 2; NN_val = 1; NOBR = 2; // NOBR must be > N if INIT='L'/'B'
        ITMAX1 = 1; ITMAX2 = 1; // Minimal iterations for validation runs
        SetUpBase(false); 
    }
};


TEST_F(IB03ADTestParamValidation, ParameterValidation) {
    int n_temp = N_val; int lx_temp = LX_val;
    // Ensure X_params is sized correctly for this minimal N_val, M_val, etc.
    // LX = (NN*(L+2) + 1)*L + N*(L+M+1) + L*M
    // LX = (1*(1+2)+1)*1 + 1*(1+1+1) + 1*1 = 4*1 + 1*3 + 1 = 4+3+1 = 8
    ASSERT_EQ(LX_val, 8) << "LX_val calculation mismatch for param validation setup.";
    std::vector<double> x_temp(LX_val > 0 ? LX_val : 1, 0.0); // Use actual LX_val
    std::vector<double> u_temp( (size_t)NSMP_val * M_val, 0.0);
    std::vector<double> y_temp( (size_t)NSMP_val * L_val, 0.0);


    // Invalid INIT
    info_result = slicot_ib03ad('X', ALG, STOR, NOBR, M_val, L_val, NSMP_val, &n_temp, NN_val, ITMAX1, ITMAX2, NPRINT, u_temp.data(), LDU_c, y_temp.data(), LDY_c, x_temp.data(), &lx_temp, TOL1, TOL2, nullptr, nullptr, nullptr, 0);
    EXPECT_EQ(info_result, -1);
    n_temp = N_val; lx_temp = LX_val; // Reset for next test

    // Invalid ALG
    info_result = slicot_ib03ad(INIT, 'X', STOR, NOBR, M_val, L_val, NSMP_val, &n_temp, NN_val, ITMAX1, ITMAX2, NPRINT, u_temp.data(), LDU_c, y_temp.data(), LDY_c, x_temp.data(), &lx_temp, TOL1, TOL2, nullptr, nullptr, nullptr, 0);
    EXPECT_EQ(info_result, -2);
    n_temp = N_val; lx_temp = LX_val;

    // Invalid STOR with ALG='D'
    info_result = slicot_ib03ad(INIT, 'D', 'X', NOBR, M_val, L_val, NSMP_val, &n_temp, NN_val, ITMAX1, ITMAX2, NPRINT, u_temp.data(), LDU_c, y_temp.data(), LDY_c, x_temp.data(), &lx_temp, TOL1, TOL2, nullptr, nullptr, nullptr, 0);
    EXPECT_EQ(info_result, -3);
    n_temp = N_val; lx_temp = LX_val;
    
    // Invalid M
    info_result = slicot_ib03ad(INIT, ALG, STOR, NOBR, -1, L_val, NSMP_val, &n_temp, NN_val, ITMAX1, ITMAX2, NPRINT, u_temp.data(), LDU_c, y_temp.data(), LDY_c, x_temp.data(), &lx_temp, TOL1, TOL2, nullptr, nullptr, nullptr, 0);
    EXPECT_EQ(info_result, -5);
    n_temp = N_val; lx_temp = LX_val;
}

class IB03ADTestZeroDim : public IB03ADTest {
protected:
    void SetUp() override {
        // Parameters will be set per test case
        SetUpBase(false); // Do not load CSV initially
    }
};

TEST_F(IB03ADTestZeroDim, ZeroM_InitN) {
    INIT = 'N'; 
    M_val = 0;
    L_val = 1; 
    N_val = 0;
    NN_val = 1;
    NSMP_val = 5; 
    NOBR = 2; 
    ITMAX1 = 1; ITMAX2 = 1;

    // Recalculate LX for M=0
    int bsn_zero_m = NN_val * (L_val + 2) + 1; // 1*(1+2)+1 = 4
    int num_nl_zero_m = bsn_zero_m * L_val;   // 4*1 = 4
    int num_l_zero_m = N_val * (L_val + M_val + 1) + L_val * M_val; // 1*(1+0+1) + 1*0 = 2
    LX_val = num_nl_zero_m + num_l_zero_m; // 4 + 2 = 6
    ASSERT_GT(LX_val, 0);

    X_params.assign(LX_val, 0.1); 
    U_data.clear(); 
    Y_data.assign((size_t)NSMP_val * L_val, 0.2); 

    n_returned = N_val;
    lx_returned = LX_val;
    
    // For M=0, LDU passed to C wrapper should be >=1.
    // The Fortran routine expects LDU >= 1 if M=0.
    // Passing MAX(1, NSMP_val) for LDU to be safe, though 1 should suffice.
    int ldu_for_m_zero = MAX(1, NSMP_val);


    info_result = slicot_ib03ad(INIT, ALG, STOR, NOBR, M_val, L_val, NSMP_val,
                                &n_returned, NN_val, ITMAX1, ITMAX2, NPRINT,
                                nullptr, ldu_for_m_zero, 
                                Y_data.empty() ? nullptr : Y_data.data(), MAX(1,NSMP_val), 
                                X_params.data(), &lx_returned,
                                TOL1, TOL2,
                                iwork_summary.data(), dwork_summary.data(),
                                &iwarn_result, 0 /* row_major */);
    
    EXPECT_GE(info_result, 0) << "Expected non-negative info_result for M=0, INIT='N'";
     if (info_result == 0) {
        EXPECT_EQ(iwarn_result, 1); // Or other acceptable warning code
    }
}

TEST_F(IB03ADTestZeroDim, ZeroN_InitS) {
    INIT = 'S'; 
    N_val = 0;  
    M_val = 1;
    L_val = 1;
    NN_val = 1;
    NSMP_val = 5;
    NOBR = 1; 
    ITMAX1 = 1; ITMAX2 = 1;


    int bsn_zero_n = NN_val * (L_val + 2) + 1; // 1*(1+2)+1 = 4
    int num_nl_zero_n = bsn_zero_n * L_val;   // 4*1 = 4
    int num_l_zero_n = N_val * (L_val + M_val + 1) + L_val * M_val; // 0*(1+1+1)+1*1 = 1
    LX_val = num_nl_zero_n + num_l_zero_n; // 4 + 1 = 5
    ASSERT_GT(LX_val, 0);

    X_params.assign(LX_val, 0.1); 
    U_data.assign((size_t)NSMP_val * M_val, 0.3);
    Y_data.assign((size_t)NSMP_val * L_val, 0.4);


    n_returned = N_val;
    lx_returned = LX_val;

    info_result = slicot_ib03ad(INIT, ALG, STOR, NOBR, M_val, L_val, NSMP_val,
                                &n_returned, NN_val, ITMAX1, ITMAX2, NPRINT,
                                U_data.empty() ? nullptr : U_data.data(), MAX(1,NSMP_val), 
                                Y_data.empty() ? nullptr : Y_data.data(), MAX(1,NSMP_val), 
                                X_params.data(), &lx_returned,
                                TOL1, TOL2,
                                iwork_summary.data(), dwork_summary.data(),
                                &iwarn_result, 0 /* row_major */);

    EXPECT_GE(info_result, 0) << "Wrapper validation failed or Fortran routine indicated fatal error for N=0, INIT='S'. INFO: " << info_result;
    if (info_result == 0 && iwarn_result != 0) {
        // A warning is acceptable if INFO is 0.
        std::cout << "ZeroN_InitS: INFO=0, IWARN=" << iwarn_result 
                  << " (Numerical warning from Fortran, acceptable for wrapper test)" << std::endl;
    } else if (info_result == 0 && iwarn_result == 0) {
        // Also acceptable
    } else if (info_result != 0) {
        // This will be caught by EXPECT_GE if info_result < 0
    }
}

