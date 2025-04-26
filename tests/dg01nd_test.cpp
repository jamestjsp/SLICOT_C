#include <gtest/gtest.h>  
#include <vector>  
#include <cmath>  
#include <string>  
#include <stdexcept>
#include <algorithm>

#include "dg01nd.h"
#include "slicot_utils.h"
#include "test_utils.h"

// Path to the CSV test data file  
const std::string DATA_FILE_PATH = "data/dg01nd.csv";

// --- Column-Major Test Fixture ---  
class DG01NDTestColMajor : public ::testing::Test {  
protected:  
    // Test parameters based on example
    int N = 8; // Half the number of real samples
    char INDI = 'D'; // Direct Fourier transform
    
    // Number of samples - will be set in SetUp
    int NSMP = 0; 

    // Column names to load from CSV (**MUST match CSV header exactly**)  
    std::vector<std::string> input_columns = {"A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", 
                                             "A9", "A10", "A11", "A12", "A13", "A14", "A15", "A16"};  
    std::vector<std::string> output_columns = {};

    // Verification tolerance  
    double check_tol = 1e-4;

    // Input arrays for the function
    std::vector<double> XR; // Odd part of the input signal (size N)
    std::vector<double> XI; // Even part of the input signal (size N)
    
    // Raw input data from CSV (will contain all A values)
    std::vector<double> A; 
    // Dummy output vector (needed for CSV loader interface)
    std::vector<double> Y_dummy;

    // Expected results (from the program results section)  
    std::vector<double> XR_expected = {4.0420, -3.1322, 0.1862, -2.1312, 1.5059, 2.1927, -1.4462, -0.5757, -0.2202};
    std::vector<double> XI_expected = {0.0000, -0.2421, -1.4675, -1.1707, -1.3815, -0.1908, 2.0327, 1.4914, 0.0000};
    int expected_info = 0;  

    // Result variables
    int info_result = -999; 

    void SetUp() override {  
        // --- Load Sample Data ---  
        int samples_loaded = 0;  

        try {  
            bool success = load_test_data_from_csv(  
                DATA_FILE_PATH, input_columns, output_columns,  
                A, Y_dummy, samples_loaded  
            );  
            ASSERT_TRUE(success) << "CSV loading reported failure for " << DATA_FILE_PATH;  
            ASSERT_GT(samples_loaded, 0) << "No data samples loaded from CSV: " << DATA_FILE_PATH;  
            NSMP = samples_loaded;

        } catch (const std::runtime_error& e) {  
            FAIL() << "CSV data loading failed: " << e.what();  
        } catch (...) {  
            FAIL() << "Caught unknown exception during CSV data loading.";  
        }

        // Split A into XR (odd indices in Fortran) and XI (even indices in Fortran)
        XR.resize(N);
        XI.resize(N);
        
        for (int i = 0; i < N; i++) {
            XR[i] = A[2*i];      // A[0], A[2], ... - Fortran's odd indices A(1), A(3), ...
            XI[i] = A[2*i + 1];  // A[1], A[3], ... - Fortran's even indices A(2), A(4), ...
        }
        
        // Resize XR and XI to accommodate output which will be N+1 for INDI='D'
        // This is only needed for the direct transform
        if (INDI == 'D' || INDI == 'd') {
            XR.resize(N+1);
            XI.resize(N+1);
        }
    }  
};

// --- Row-Major Test Fixture --- (same as col-major since we're dealing with vectors)
class DG01NDTestRowMajor : public DG01NDTestColMajor {  
protected:  
    void SetUp() override {  
        // Use the col-major setup since we're just dealing with vectors
        DG01NDTestColMajor::SetUp();
    }  
};

// --- Test Cases ---

// Test: Documentation Example (Column-Major)  
TEST_F(DG01NDTestColMajor, DocExample) {  
    // Call C wrapper function
    info_result = slicot_dg01nd(INDI, N, XR.data(), XI.data(), 0 /* col_major */);

    // Verify return code  
    ASSERT_EQ(info_result, expected_info);

    // Verify output results against expected values
    ASSERT_EQ(XR.size(), XR_expected.size()) << "XR size mismatch";
    ASSERT_EQ(XI.size(), XI_expected.size()) << "XI size mismatch";
    
    for (size_t i = 0; i < XR_expected.size(); i++) {
        EXPECT_NEAR(XR[i], XR_expected[i], check_tol) << "XR mismatch at index " << i;
    }
    
    for (size_t i = 0; i < XI_expected.size(); i++) {
        EXPECT_NEAR(XI[i], XI_expected[i], check_tol) << "XI mismatch at index " << i;
    }
}

// Test: Documentation Example (Row-Major) - for vectors, should be the same as col-major
TEST_F(DG01NDTestRowMajor, DocExample) {  
    // Call C wrapper function with row_major=1
    info_result = slicot_dg01nd(INDI, N, XR.data(), XI.data(), 1 /* row_major */);

    // Verify return code  
    ASSERT_EQ(info_result, expected_info);

    // Verify output results against expected values
    ASSERT_EQ(XR.size(), XR_expected.size()) << "XR size mismatch";
    ASSERT_EQ(XI.size(), XI_expected.size()) << "XI size mismatch";
    
    for (size_t i = 0; i < XR_expected.size(); i++) {
        EXPECT_NEAR(XR[i], XR_expected[i], check_tol) << "XR mismatch at index " << i;
    }
    
    for (size_t i = 0; i < XI_expected.size(); i++) {
        EXPECT_NEAR(XI[i], XI_expected[i], check_tol) << "XI mismatch at index " << i;
    }
}

// Test: Parameter Validation
TEST_F(DG01NDTestColMajor, ParameterValidation) {  
    std::vector<double> dummy_xr(2), dummy_xi(2);
    
    // Test invalid N (too small)
    info_result = slicot_dg01nd(INDI, 1, dummy_xr.data(), dummy_xi.data(), 0);  
    EXPECT_EQ(info_result, -2);

    // Test N not a power of 2
    info_result = slicot_dg01nd(INDI, 3, dummy_xr.data(), dummy_xi.data(), 0);  
    EXPECT_EQ(info_result, -2);

    // Test invalid INDI
    info_result = slicot_dg01nd('X', N, dummy_xr.data(), dummy_xi.data(), 0);  
    EXPECT_EQ(info_result, -1);

    // Test NULL pointers
    info_result = slicot_dg01nd(INDI, N, nullptr, dummy_xi.data(), 0);  
    EXPECT_EQ(info_result, -3);
    
    info_result = slicot_dg01nd(INDI, N, dummy_xr.data(), nullptr, 0);  
    EXPECT_EQ(info_result, -4);
}

// Test inverse transform
TEST_F(DG01NDTestColMajor, InverseTransform) {
    // First do a direct transform
    info_result = slicot_dg01nd(INDI, N, XR.data(), XI.data(), 0 /* col_major */);
    ASSERT_EQ(info_result, expected_info);
    
    // Save the transform result
    std::vector<double> xr_transform(XR.begin(), XR.end());
    std::vector<double> xi_transform(XI.begin(), XI.end());
    
    // Now do an inverse transform on the result
    char indi_inverse = 'I'; // Switch to inverse
    info_result = slicot_dg01nd(indi_inverse, N, xr_transform.data(), xi_transform.data(), 0);
    ASSERT_EQ(info_result, expected_info);
    
    // The result should be proportional to the original signal
    // In DG01ND, the scaling factor is 2N
    double scale_factor = 2.0 * N;
    
    // Original data is in A
    // For INDI='I', output is size N and contains odd/even parts
    for (int i = 0; i < N; i++) {
        EXPECT_NEAR(xr_transform[i] / scale_factor, A[2*i], check_tol) 
            << "Inverse XR mismatch at index " << i;
        EXPECT_NEAR(xi_transform[i] / scale_factor, A[2*i+1], check_tol) 
            << "Inverse XI mismatch at index " << i;
    }
}
