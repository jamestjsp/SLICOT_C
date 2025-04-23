#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm> // For std::max

#include "ib01ad.h"
#include "slicot_utils.h" // For transpose functions
#include "test_utils.h"   // For data loading utilities

// Path to test data file
const std::string DATA_FILE_PATH = "data/ib01ad.csv";

// Column-major test fixture
class IB01ADTestColMajor : public ::testing::Test {
protected:
    // Parameters from IB01AD.dat example
    int NOBR = 15;
    int M = 1;          // Number of inputs (Fortran example reads U if M>0)
    int L = 1;          // Number of outputs
    int NSMP = 1000;    // Number of samples
    double RCOND = 0.0; // Use default tolerance for rank determination
    double TOL = -1.0;  // Use default tolerance for order detection (-1.0 means largest gap)
    char METH = 'M';    // MOESP algorithm
    char ALG = 'C';     // Cholesky algorithm
    char JOBD = 'N';    // Don't compute B/D using MOESP (from .dat file)
    char BATCH = 'O';   // One batch only
    char CONCT = 'N';   // No connection between blocks
    char CTRL = 'N';    // No user confirmation for order N
    
    double check_tol = 1e-4; // Tolerance for checking singular values

    // Input data (U and Y)
    std::vector<double> U; // Input data (M columns)
    std::vector<double> Y; // Output data (L columns) - from IB01AD.dat

    // Expected results from IB01AD.res
    int expected_N = 4; 
    std::vector<double> SV_expected = {
        69.8841, 14.9963, 3.6675, 1.9677, 0.3000, 0.2078, 0.1651, 0.1373,
        0.1133,  0.1059,  0.0856, 0.0784, 0.0733, 0.0678, 0.0571
    };

    // Output arrays
    std::vector<double> R;   // Triangular factor
    std::vector<double> SV;  // Singular values

    // Workspace arrays
    std::vector<int> IWORK;
    std::vector<double> DWORK;
    int LDWORK = 100000; // Use large workspace size to avoid query issues

    // Result variables
    int n_result = -1;
    int iwarn_result = -1;
    
    void SetUp() override {
        // --- Load data from CSV file ---
        bool loaded = false;
        try {
            loaded = load_test_data_from_csv(DATA_FILE_PATH, U, Y);
        } catch (const std::exception& e) {
            std::cerr << "Exception loading CSV: " << e.what() << std::endl;
            loaded = false;
        }
        
        // If loading fails, use hardcoded data as fallback
        if (!loaded || U.size() < NSMP || Y.size() < NSMP) {
            std::cerr << "Warning: Failed to load data from CSV, using hardcoded data" << std::endl;
            
            // Initialize U with hardcoded data from IB01AD.dat
            U.resize(NSMP);
            for (int i = 0; i < NSMP; i++) {
                U[i] = (i % 2 == 0) ? 6.41 : 3.41;
            }
            
            // Initialize Y with a simple pattern
            Y.resize(NSMP);
            for (int i = 0; i < NSMP; i++) {
                Y[i] = 4.75 + 1.25 * sin(i * 0.01);
            }
        }
        
        // --- Array Sizing ---
        // LDR calculation (Fortran leading dimension for R)
        int LDR = std::max(1, 2 * (M + L) * NOBR); // Fortran row count
        
        // R array size: LDR rows, 2*(M + L)*NOBR columns
        R.resize((size_t)LDR * 2 * (M + L) * NOBR);

        // SV array size: L*NOBR
        SV.resize((size_t)L * NOBR);

        // IWORK size
        int LIWORK = std::max(3, (M + L) * NOBR); // Use larger size for safety
        IWORK.resize(LIWORK);

        // Workspace array
        DWORK.resize(LDWORK);
    }
};

// Row-major test fixture
class IB01ADTestRowMajor : public ::testing::Test {
protected:
    // Parameters from IB01AD.dat example
    int NOBR = 15;
    int M = 1;          // Number of inputs (Fortran example reads U if M>0)
    int L = 1;          // Number of outputs
    int NSMP = 1000;    // Number of samples
    double RCOND = 0.0; // Use default tolerance for rank determination
    double TOL = -1.0;  // Use default tolerance for order detection (-1.0 means largest gap)
    char METH = 'M';    // MOESP algorithm
    char ALG = 'C';     // Cholesky algorithm
    char JOBD = 'N';    // Don't compute B/D using MOESP (from .dat file)
    char BATCH = 'O';   // One batch only
    char CONCT = 'N';   // No connection between blocks
    char CTRL = 'N';    // No user confirmation for order N
    
    double check_tol = 1e-4; // Tolerance for checking singular values

    // Input data (U and Y) in row-major format
    std::vector<double> U; // Input data in row-major format
    std::vector<double> Y; // Output data in row-major format

    // Expected results from IB01AD.res
    int expected_N = 4; 
    std::vector<double> SV_expected = {
        69.8841, 14.9963, 3.6675, 1.9677, 0.3000, 0.2078, 0.1651, 0.1373,
        0.1133,  0.1059,  0.0856, 0.0784, 0.0733, 0.0678, 0.0571
    };

    // Output arrays
    std::vector<double> R;   // Triangular factor
    std::vector<double> SV;  // Singular values

    // Workspace arrays
    std::vector<int> IWORK;
    std::vector<double> DWORK;
    int LDWORK = 100000; // Use large workspace size to avoid query issues

    // Result variables
    int n_result = -1;
    int iwarn_result = -1;
    
    void SetUp() override {
        // --- Load column-major data from CSV ---
        std::vector<double> U_col, Y_col;
        bool loaded = false;
        
        try {
            loaded = load_test_data_from_csv(DATA_FILE_PATH, U_col, Y_col);
        } catch (const std::exception& e) {
            std::cerr << "Exception loading CSV: " << e.what() << std::endl;
            loaded = false;
        }
        
        // If loading fails, use hardcoded data as fallback
        if (!loaded || U_col.size() < NSMP || Y_col.size() < NSMP) {
            std::cerr << "Warning: Failed to load data from CSV, using hardcoded data" << std::endl;
            
            // Initialize U with hardcoded data
            U_col.resize(NSMP);
            for (int i = 0; i < NSMP; i++) {
                U_col[i] = (i % 2 == 0) ? 6.41 : 3.41;
            }
            
            // Initialize Y with a simple pattern
            Y_col.resize(NSMP);
            for (int i = 0; i < NSMP; i++) {
                Y_col[i] = 4.75 + 1.25 * sin(i * 0.01);
            }
        }
        
        // --- Convert to row-major format ---
        // Allocate row-major arrays
        U.resize(U_col.size());
        Y.resize(Y_col.size());
        
        // Convert U and Y to row-major format
        slicot_transpose_to_c(U_col.data(), U.data(), NSMP, M, sizeof(double));
        slicot_transpose_to_c(Y_col.data(), Y.data(), NSMP, L, sizeof(double));
        
        // --- Array Sizing ---
        // LDR calculation for row-major
        int LDR_rm = 2 * (M + L) * NOBR;
        
        // R array size: 2*(M + L)*NOBR rows, LDR_rm columns (for row-major)
        R.resize((size_t)2 * (M + L) * NOBR * LDR_rm);

        // SV array size: L*NOBR
        SV.resize((size_t)L * NOBR);

        // IWORK size
        int LIWORK = std::max(3, (M + L) * NOBR);
        IWORK.resize(LIWORK);

        // Workspace array
        DWORK.resize(LDWORK);
    }
};

// Test using data from IB01AD.dat/.res (Column-Major)
TEST_F(IB01ADTestColMajor, DocExample) {
    // Run the actual computation with the large workspace
    int ldu = NSMP;  // Column-major leading dimension
    int ldy = NSMP;
    int ldr = 2 * (M + L) * NOBR;
    int n = 0;  // Auto-detect order
    
    int info = slicot_ib01ad(
        METH, ALG, JOBD, BATCH, CONCT, CTRL,
        NOBR, M, L, NSMP,
        U.data(), ldu, Y.data(), ldy,
        &n, R.data(), ldr, SV.data(), RCOND,
        TOL, IWORK.data(), DWORK.data(), LDWORK,
        &iwarn_result, 0  // Column-major
    );

    // Save result
    n_result = n;
    
    // Verify results
    ASSERT_EQ(info, 0) << "slicot_ib01ad failed with info = " << info;
    EXPECT_EQ(iwarn_result, 0) << "slicot_ib01ad returned warning = " << iwarn_result;
    EXPECT_EQ(n_result, expected_N) << "System order N does not match expected value.";

    // Check singular values
    ASSERT_EQ(SV.size(), SV_expected.size()) << "Output SV size mismatch.";
    for (size_t i = 0; i < SV.size(); ++i) {
        EXPECT_NEAR(SV[i], SV_expected[i], check_tol) 
            << "Mismatch in singular value at index " << i;
    }
}

// Test using data from IB01AD.dat/.res (Row-Major)
TEST_F(IB01ADTestRowMajor, DocExample) {
    // Run the actual computation with the large workspace
    int ldu = M;  // Row-major leading dimension
    int ldy = L;
    int ldr = 2 * (M + L) * NOBR;
    int n = 0;  // Auto-detect order
    
    int info = slicot_ib01ad(
        METH, ALG, JOBD, BATCH, CONCT, CTRL,
        NOBR, M, L, NSMP,
        U.data(), ldu, Y.data(), ldy,
        &n, R.data(), ldr, SV.data(), RCOND,
        TOL, IWORK.data(), DWORK.data(), LDWORK,
        &iwarn_result, 1  // Row-major
    );
    
    // Save result
    n_result = n;
    
    // Verify results
    ASSERT_EQ(info, 0) << "slicot_ib01ad failed with info = " << info;
    EXPECT_EQ(iwarn_result, 0) << "slicot_ib01ad returned warning = " << iwarn_result;
    EXPECT_EQ(n_result, expected_N) << "System order N does not match expected value.";

    // Check singular values
    ASSERT_EQ(SV.size(), SV_expected.size()) << "Output SV size mismatch.";
    for (size_t i = 0; i < SV.size(); ++i) {
        EXPECT_NEAR(SV[i], SV_expected[i], check_tol) 
            << "Mismatch in singular value at index " << i;
    }
}

// Test parameter validation (using Column-Major fixture)
TEST_F(IB01ADTestColMajor, ParameterValidation) {
    int n_out = 0;
    int iwarn = 0;
    int info;
    int dummy_ldwork = 10; 
    DWORK.resize(dummy_ldwork);
    int ldu = NSMP;
    int ldy = NSMP;
    int ldr = 2 * (M + L) * NOBR;
    
    // Test invalid METH parameter
    info = slicot_ib01ad('X', ALG, JOBD, BATCH, CONCT, CTRL,
                        NOBR, M, L, NSMP,
                        U.data(), ldu, Y.data(), ldy,
                        &n_out, R.data(), ldr, SV.data(), RCOND,
                        TOL, IWORK.data(), DWORK.data(), dummy_ldwork,
                        &iwarn, 0);  // Column-major
    EXPECT_LT(info, 0) << "Invalid METH parameter not detected";
    
    // Test invalid ALG parameter
    info = slicot_ib01ad(METH, 'X', JOBD, BATCH, CONCT, CTRL,
                        NOBR, M, L, NSMP,
                        U.data(), ldu, Y.data(), ldy,
                        &n_out, R.data(), ldr, SV.data(), RCOND,
                        TOL, IWORK.data(), DWORK.data(), dummy_ldwork,
                        &iwarn, 0);  // Column-major
    EXPECT_LT(info, 0) << "Invalid ALG parameter not detected";
    
    // Test invalid NOBR parameter
    info = slicot_ib01ad(METH, ALG, JOBD, BATCH, CONCT, CTRL,
                        0, M, L, NSMP,
                        U.data(), ldu, Y.data(), ldy,
                        &n_out, R.data(), ldr, SV.data(), RCOND,
                        TOL, IWORK.data(), DWORK.data(), dummy_ldwork,
                        &iwarn, 0);  // Column-major
    EXPECT_LT(info, 0) << "Invalid NOBR parameter not detected";
    
    // Test invalid L parameter
    info = slicot_ib01ad(METH, ALG, JOBD, BATCH, CONCT, CTRL,
                        NOBR, M, 0, NSMP,
                        U.data(), ldu, Y.data(), ldy,
                        &n_out, R.data(), ldr, SV.data(), RCOND,
                        TOL, IWORK.data(), DWORK.data(), dummy_ldwork,
                        &iwarn, 0);  // Column-major
    EXPECT_LT(info, 0) << "Invalid L parameter not detected";
    
    // Test invalid NSMP parameter for non-sequential processing
    info = slicot_ib01ad(METH, ALG, JOBD, 'O', CONCT, CTRL,
                        NOBR, M, L, 1, // NSMP too small
                        U.data(), ldu, Y.data(), ldy,
                        &n_out, R.data(), ldr, SV.data(), RCOND,
                        TOL, IWORK.data(), DWORK.data(), dummy_ldwork,
                        &iwarn, 0);  // Column-major
    EXPECT_LT(info, 0) << "Invalid NSMP parameter not detected";
}