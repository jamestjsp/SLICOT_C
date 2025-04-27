#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>

#include "de01od.h"
#include "test_utils.h"
#include "test_config.h" // Include the CMake-generated configuration

// Use the TEST_DATA_DIR macro defined in test_config.h
const std::string DATA_FILE_PATH = TEST_DATA_DIR "de01od.csv";

// --- Column-Major Test Fixture ---
class DE01ODTestColMajor : public ::testing::Test {
protected:
    // Test parameters
    int N = 8;           // Number of samples
    char CONV = 'C';     // Convolution
    int NSMP = 0;        // Will be determined by loader in SetUp
    
    // Column names to load from CSV (MUST match CSV header exactly)
    std::vector<std::string> input_columns = {"A", "B"};
    std::vector<std::string> output_columns = {"A_expected"};
    
    // Verification tolerance
    double check_tol = 1e-4;
    
    // Input/Output data vectors
    std::vector<double> A;  // First signal (input/output)
    std::vector<double> B;  // Second signal (input)
    std::vector<double> A_expected;  // Expected output
    
    // Results
    int info_result = -999;  // Initialize to indicate not run
    
    // SetUp method: Load data, initialize inputs, size outputs
    void SetUp() override {
        // Load test data from CSV
        int samples_loaded = 0;
        
        try {
            bool success = load_test_data_from_csv(
                DATA_FILE_PATH, input_columns, output_columns,
                A, A_expected, samples_loaded
            );
            ASSERT_TRUE(success) << "CSV loading reported failure for " << DATA_FILE_PATH;
            ASSERT_GT(samples_loaded, 0) << "No data samples loaded from CSV: " << DATA_FILE_PATH;
            
            // Update NSMP based on actual data loaded
            NSMP = samples_loaded;
            
            // Extract column B from A (the loader puts all input columns in A)
            B.resize(NSMP);
            for (int i = 0; i < NSMP; i++) {
                B[i] = A[NSMP + i];  // B starts after A in the loaded data
            }
            
            // Resize A to contain only the first signal data
            A.resize(NSMP);
            
        } catch (const std::runtime_error& e) {
            FAIL() << "CSV data loading failed: " << e.what();
        } catch (...) {
            FAIL() << "Caught unknown exception during CSV data loading.";
        }
    }
};

// --- Row-Major Test Fixture ---
class DE01ODTestRowMajor : public DE01ODTestColMajor {
protected:
    // Row-major versions of vectors
    std::vector<double> A_rm;
    std::vector<double> B_rm;
    
    void SetUp() override {
        // Load column-major data using parent SetUp
        DE01ODTestColMajor::SetUp();
        
        // Convert column-major A and B to row-major (for this simple case, just copy)
        A_rm = A;
        B_rm = B;
    }
};

// --- Test Cases ---

// Test: Documentation Example (Column-Major)
TEST_F(DE01ODTestColMajor, DocExample) {
    // Make a copy of B since it will be overwritten
    std::vector<double> B_copy = B;
    
    // Call C wrapper function
    info_result = slicot_de01od(CONV, N, A.data(), B_copy.data(), 0);
    
    // Verify return code
    ASSERT_EQ(info_result, 0);
    
    // Verify output against expected
    ASSERT_EQ(A.size(), A_expected.size());
    for (size_t i = 0; i < A_expected.size(); ++i) {
        EXPECT_NEAR(A[i], A_expected[i], check_tol) << "Mismatch at index " << i;
    }
}

// Test: Documentation Example (Row-Major)
TEST_F(DE01ODTestRowMajor, DocExample) {
    // Make a copy of B_rm since it will be overwritten
    std::vector<double> B_rm_copy = B_rm;
    
    // Call C wrapper function with row_major=1
    info_result = slicot_de01od(CONV, N, A_rm.data(), B_rm_copy.data(), 1);
    
    // Verify return code
    ASSERT_EQ(info_result, 0);
    
    // Verify output against expected
    ASSERT_EQ(A_rm.size(), A_expected.size());
    for (size_t i = 0; i < A_expected.size(); ++i) {
        EXPECT_NEAR(A_rm[i], A_expected[i], check_tol) << "Mismatch at index " << i;
    }
}

// Test: Deconvolution (Column-Major)
TEST_F(DE01ODTestColMajor, Deconvolution) {
    // For deconvolution, we'll start with the result of convolution
    // and try to recover the original signal
    
    // Make a copy of vectors for deconvolution test
    std::vector<double> A_conv = A_expected;  // Start with the convolution result
    std::vector<double> B_copy = B;          // Second signal stays the same
    
    // Call C wrapper function with CONV = 'D' (deconvolution)
    info_result = slicot_de01od('D', N, A_conv.data(), B_copy.data(), 0);
    
    // Verify return code
    ASSERT_EQ(info_result, 0);
    
    // Verify output against original A (tolerance may need to be larger for numerical reasons)
    double deconv_tol = 1e-2;
    ASSERT_EQ(A_conv.size(), A.size());
    for (size_t i = 0; i < A.size(); ++i) {
        EXPECT_NEAR(A_conv[i], A[i], deconv_tol) << "Deconvolution mismatch at index " << i;
    }
}

// Test: Parameter Validation
TEST_F(DE01ODTestColMajor, ParameterValidation) {
    // Create dummy variables for validation tests
    std::vector<double> dummy_a(N, 0.0);
    std::vector<double> dummy_b(N, 0.0);
    
    // Test invalid N (not a power of 2)
    info_result = slicot_de01od(CONV, 3, dummy_a.data(), dummy_b.data(), 0);
    EXPECT_EQ(info_result, -2);
    
    // Test N < 2
    info_result = slicot_de01od(CONV, 1, dummy_a.data(), dummy_b.data(), 0);
    EXPECT_EQ(info_result, -2);
    
    // Test invalid CONV character
    info_result = slicot_de01od('X', N, dummy_a.data(), dummy_b.data(), 0);
    EXPECT_EQ(info_result, -1);
    
    // Test NULL pointer for A
    info_result = slicot_de01od(CONV, N, nullptr, dummy_b.data(), 0);
    EXPECT_EQ(info_result, -3);
    
    // Test NULL pointer for B
    info_result = slicot_de01od(CONV, N, dummy_a.data(), nullptr, 0);
    EXPECT_EQ(info_result, -4);
}
