#include <gtest/gtest.h>  
#include <vector>  
#include <cmath>  
#include <string>  
#include <stdexcept>
#include <algorithm>

#include "dk01md.h"
#include "slicot_utils.h"
#include "test_utils.h"

// Path to the CSV test data file  
const std::string DATA_FILE_PATH = "data/dk01md.csv";

// --- Column-Major Test Fixture ---  
class DK01MDTestColMajor : public ::testing::Test {  
protected:  
    // Test parameters
    int N = 8; // Length of vector A (number of samples)
    char TYPE = 'M'; // Hamming window
    int NSMP = 0;

    // Column names to load from CSV (**MUST match CSV header exactly**)  
    std::vector<std::string> input_columns = {"X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8"};
    std::vector<std::string> output_columns = {}; // No output columns needed

    // Verification tolerance  
    double check_tol = 1e-4;

    // Input/Output data vectors (column-major)  
    std::vector<double> A; // Input signal that will be modified
    std::vector<double> Y_dummy; // Dummy output vector for CSV loader
    
    // Expected results - from program results section
    std::vector<double> A_expected = {0.3262, 0.8326, -0.6591, 0.4286, -0.0754, 0.0820, 0.0661, -0.0262};
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

            // Transform the data from row format to column format if needed
            // With this CSV format, we might get a single row with N columns
            // We need to reshape it to get N elements in a single vector
            if (A.size() != N) {
                // If we got multiple columns, flatten them into a single vector
                std::vector<double> tmp_A = A;
                A.clear();
                A.resize(N);
                int copied = 0;
                
                // Copy as much data as we can from the loaded array
                for (size_t i = 0; i < tmp_A.size() && copied < N; ++i) {
                    A[copied++] = tmp_A[i];
                }
                
                // If we don't have enough data, fill with zeros
                while (copied < N) {
                    A[copied++] = 0.0;
                }
            }

        } catch (const std::runtime_error& e) {  
            FAIL() << "CSV data loading failed: " << e.what();  
        } catch (...) {  
            FAIL() << "Caught unknown exception during CSV data loading.";  
        }
    }  
};

// --- Row-Major Test Fixture --- (same as column-major for vectors)
class DK01MDTestRowMajor : public DK01MDTestColMajor {  
protected:  
    void SetUp() override {  
        // Use the col-major setup since we're dealing with vectors
        DK01MDTestColMajor::SetUp();
    }  
};

// --- Test Cases ---

// Test: Documentation Example (Column-Major)  
TEST_F(DK01MDTestColMajor, DocExample) {  
    // Call C wrapper function
    info_result = slicot_dk01md(TYPE, N, A.data(), 0 /* col_major */);

    // Verify return code  
    ASSERT_EQ(info_result, expected_info);

    // Verify output results against expected values
    ASSERT_EQ(A.size(), A_expected.size());
    for (size_t i = 0; i < A_expected.size(); i++) {
        EXPECT_NEAR(A[i], A_expected[i], check_tol) << "A mismatch at index " << i;
    }
}

// Test: Documentation Example (Row-Major) - for vectors, should be the same as col-major
TEST_F(DK01MDTestRowMajor, DocExample) {  
    // Call C wrapper function with row_major=1
    info_result = slicot_dk01md(TYPE, N, A.data(), 1 /* row_major */);

    // Verify return code  
    ASSERT_EQ(info_result, expected_info);

    // Verify output results against expected values
    ASSERT_EQ(A.size(), A_expected.size());
    for (size_t i = 0; i < A_expected.size(); i++) {
        EXPECT_NEAR(A[i], A_expected[i], check_tol) << "A mismatch at index " << i;
    }
}

// Test: Parameter Validation
TEST_F(DK01MDTestColMajor, ParameterValidation) {  
    std::vector<double> dummy_a(5);
    
    // Test invalid N (negative)
    info_result = slicot_dk01md(TYPE, -1, dummy_a.data(), 0);  
    EXPECT_EQ(info_result, -2);

    // Test invalid N (zero)
    info_result = slicot_dk01md(TYPE, 0, dummy_a.data(), 0);  
    EXPECT_EQ(info_result, -2);

    // Test NULL pointer
    info_result = slicot_dk01md(TYPE, N, nullptr, 0);  
    EXPECT_EQ(info_result, -3);
    
    // Test invalid TYPE
    info_result = slicot_dk01md('X', N, dummy_a.data(), 0);  
    EXPECT_EQ(info_result, -1);
}

// Test different window types
TEST_F(DK01MDTestColMajor, WindowTypes) {  
    // Test Hann window
    std::vector<double> A_hann = A;
    info_result = slicot_dk01md('N', N, A_hann.data(), 0);  
    EXPECT_EQ(info_result, 0);
    
    // Test Quadratic window
    std::vector<double> A_quad = A;
    info_result = slicot_dk01md('Q', N, A_quad.data(), 0);  
    EXPECT_EQ(info_result, 0);
}
