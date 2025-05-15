#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm> // For std::min, std::max
#include <iostream>
#include <iomanip>   // For std::setprecision

#include "mb05md.h"       // Header for the C wrapper
#include "slicot_utils.h" // For MAX, MIN, slicot_transpose_... etc.

// Helper function to get 2D array index for column-major storage
// Fortran: A(row, col) (1-indexed)
// C: a_flat[ col_idx * lda_c + row_idx ] (0-indexed)
size_t get_2d_col_major_idx(int r_idx, int c_idx, int lda_c) {
    // Using C-style cast as a workaround for static_cast issues
    return (size_t)c_idx * lda_c + (size_t)r_idx;
}

// Helper function to get 2D array index for row-major storage
// C: a_flat[ r_idx * ldc_c + c_idx ] (0-indexed) (ldc_c is number of columns)
size_t get_2d_row_major_idx(int r_idx, int c_idx, int num_cols_c) {
    // Using C-style cast as a workaround for static_cast issues
    return (size_t)r_idx * num_cols_c + (size_t)c_idx;
}


// --- Test Fixture ---
class MB05MDTest : public ::testing::Test {
protected:
    char BALANC_default = 'N';
    int N_default = 4;
    double DELTA_default = 1.0;
    double check_tol = 1e-4; // Based on example output precision

    // Variables for tests
    char BALANC;
    int N;
    double DELTA;
    std::vector<double> A_data, A_expected_data;
    std::vector<double> V_data, V_expected_data;
    std::vector<double> Y_data, Y_expected_data;
    std::vector<double> VALR_data, VALR_expected_data;
    std::vector<double> VALI_data, VALI_expected_data;
    int LDA, LDV, LDY;

    int info_result = -999; 

    // Initializes dimensions and resizes data vectors
    void InitializeData(int n_val, char balanc_val = 'N', double delta_val = 1.0, bool use_example_data = false) {
        N = n_val;
        BALANC = balanc_val;
        DELTA = delta_val;

        LDA = std::max(1, N);
        LDV = std::max(1, N);
        LDY = std::max(1, N);
        
        size_t matrix_elements = (size_t)N * N; // Using C-style cast for size_t conversion

        if (N > 0) {
            A_data.assign(matrix_elements, 0.0);
            V_data.assign(matrix_elements, 0.0);
            Y_data.assign(matrix_elements, 0.0);
            VALR_data.resize(N);
            VALI_data.resize(N);

            if (use_example_data && N == 4) {
                // Input A from MB05MD.html example (column-major order)
                A_data = {
                    0.5,  0.0,  2.3, -2.6,
                    0.0,  0.5, -1.4, -0.7,
                    2.3, -1.4,  0.5,  0.0,
                   -2.6, -0.7,  0.0,  0.5
                };

                // Expected output exp(A*delta) (column-major)
                A_expected_data = {
                   26.8551,  -3.2824,  18.7409, -19.4430,
                   -3.2824,   4.3474,  -5.1848,   0.2700,
                   18.7409,  -5.1848,  15.6012, -11.7228,
                  -19.4430,   0.2700, -11.7228,  15.6012
                };
                // Expected eigenvalues
                VALR_expected_data = {-3.0, 4.0, -1.0, 2.0};
                VALI_expected_data = {0.0, 0.0, 0.0, 0.0};
                // Expected eigenvectors V (column-major)
                V_expected_data = {
                  -0.7,  0.1,  0.5, -0.5,
                   0.7, -0.1,  0.5, -0.5,
                   0.1,  0.7,  0.5,  0.5,
                  -0.1, -0.7,  0.5,  0.5
                };
                 // Expected Y (column-major)
                Y_expected_data = {
                  -0.0349,  0.0050,  0.0249, -0.0249,
                  38.2187, -5.4598, 27.2991,-27.2991,
                   0.0368,  0.2575,  0.1839,  0.1839,
                  -0.7389, -5.1723,  3.6945,  3.6945
                };
            }
        } else { // N = 0
            A_data.clear(); V_data.clear(); Y_data.clear();
            VALR_data.clear(); VALI_data.clear();
            A_expected_data.clear(); V_expected_data.clear(); Y_expected_data.clear();
            VALR_expected_data.clear(); VALI_expected_data.clear();
        }
    }
};

// --- Test Cases ---

TEST_F(MB05MDTest, ParameterValidation) {
    InitializeData(2); // Base valid dimensions for most checks
    
    // BALANC invalid
    info_result = slicot_mb05md('X', N, DELTA, A_data.data(), LDA, V_data.data(), LDV, Y_data.data(), LDY, VALR_data.data(), VALI_data.data(), 0);
    EXPECT_EQ(info_result, -1);
    // N invalid
    info_result = slicot_mb05md(BALANC, -1, DELTA, nullptr, 1, nullptr, 1, nullptr, 1, nullptr, nullptr, 0);
    EXPECT_EQ(info_result, -2);
    
    // A is NULL (N>0)
    if (N > 0) {
        info_result = slicot_mb05md(BALANC, N, DELTA, nullptr, LDA, V_data.data(), LDV, Y_data.data(), LDY, VALR_data.data(), VALI_data.data(), 0);
        EXPECT_EQ(info_result, -4);
    }
    // LDA invalid (col-major)
    if (N > 0) {
        info_result = slicot_mb05md(BALANC, N, DELTA, A_data.data(), 0, V_data.data(), LDV, Y_data.data(), LDY, VALR_data.data(), VALI_data.data(), 0);
        EXPECT_EQ(info_result, -5);
    }
    // LDA invalid (row-major, lda < N)
    if (N > 1) { // Ensure N-1 is valid for LDA
        info_result = slicot_mb05md(BALANC, N, DELTA, A_data.data(), N - 1, V_data.data(), LDV, Y_data.data(), LDY, VALR_data.data(), VALI_data.data(), 1);
        EXPECT_EQ(info_result, -5);
    }


    // V is NULL (N>0)
    if (N > 0) {
        info_result = slicot_mb05md(BALANC, N, DELTA, A_data.data(), LDA, nullptr, LDV, Y_data.data(), LDY, VALR_data.data(), VALI_data.data(), 0);
        EXPECT_EQ(info_result, -6);
    }
    // Y is NULL (N>0)
    if (N > 0) {
        info_result = slicot_mb05md(BALANC, N, DELTA, A_data.data(), LDA, V_data.data(), LDV, nullptr, LDY, VALR_data.data(), VALI_data.data(), 0);
        EXPECT_EQ(info_result, -8);
    }
    // VALR is NULL (N>0)
    if (N > 0) {
        info_result = slicot_mb05md(BALANC, N, DELTA, A_data.data(), LDA, V_data.data(), LDV, Y_data.data(), LDY, nullptr, VALI_data.data(), 0);
        EXPECT_EQ(info_result, -10);
    }
    // VALI is NULL (N>0)
    if (N > 0) {
        info_result = slicot_mb05md(BALANC, N, DELTA, A_data.data(), LDA, V_data.data(), LDV, Y_data.data(), LDY, VALR_data.data(), nullptr, 0);
        EXPECT_EQ(info_result, -11);
    }
}

TEST_F(MB05MDTest, ZeroDimensionN) {
    InitializeData(0); 
    info_result = slicot_mb05md(BALANC, N, DELTA, 
                                nullptr, LDA, nullptr, LDV, nullptr, LDY, 
                                nullptr, nullptr, 0); 
    EXPECT_EQ(info_result, 0);
}

TEST_F(MB05MDTest, DocExample_ColMajor) {
    InitializeData(4, 'N', 1.0, true); // N=4, BALANC='N', DELTA=1.0, use example data

    info_result = slicot_mb05md(BALANC, N, DELTA, A_data.data(), LDA, 
                                V_data.data(), LDV, Y_data.data(), LDY, 
                                VALR_data.data(), VALI_data.data(), 0); // Column-major
    EXPECT_EQ(info_result, 0);

    if (info_result == 0 && N > 0) { // Proceed with checks only if call was successful
        for(int j=0; j<N; ++j) { // col
            for(int i=0; i<N; ++i) { // row
                EXPECT_NEAR(A_data[get_2d_col_major_idx(i,j,LDA)], A_expected_data[get_2d_col_major_idx(i,j,N)], check_tol) << "A mismatch at (" << i << "," << j << ")";
                EXPECT_NEAR(V_data[get_2d_col_major_idx(i,j,LDV)], V_expected_data[get_2d_col_major_idx(i,j,N)], check_tol) << "V mismatch at (" << i << "," << j << ")";
                EXPECT_NEAR(Y_data[get_2d_col_major_idx(i,j,LDY)], Y_expected_data[get_2d_col_major_idx(i,j,N)], check_tol) << "Y mismatch at (" << i << "," << j << ")";
            }
        }
        for(int i=0; i<N; ++i) {
            EXPECT_NEAR(VALR_data[i], VALR_expected_data[i], check_tol) << "VALR mismatch at " << i;
            EXPECT_NEAR(VALI_data[i], VALI_expected_data[i], check_tol) << "VALI mismatch at " << i;
        }
    }
}

TEST_F(MB05MDTest, DocExample_RowMajor) {
    InitializeData(4, 'N', 1.0, true); // N=4, BALANC='N', DELTA=1.0, use example data

    // Transpose input A_data (which is col-major from InitializeData) to row-major for the test input
    std::vector<double> A_input_rm(A_data.size());
    if (N > 0) { // only transpose if N > 0
      slicot_transpose_to_c_with_ld(A_data.data(), A_input_rm.data(), N, N, N, N, sizeof(double));
    }
    
    // Output arrays for row-major call
    std::vector<double> V_output_rm(V_data.size());
    std::vector<double> Y_output_rm(Y_data.size());
    std::vector<double> VALR_output_rm(VALR_data.size());
    std::vector<double> VALI_output_rm(VALI_data.size());


    info_result = slicot_mb05md(BALANC, N, DELTA, A_input_rm.data(), N, // For row-major, LDA is num_cols (which is N)
                                V_output_rm.data(), N, Y_output_rm.data(), N, // LDV, LDY are also num_cols (N)
                                VALR_output_rm.data(), VALI_output_rm.data(), 1); // Row-major
    EXPECT_EQ(info_result, 0);

    if (info_result == 0 && N > 0) { // Proceed with checks only if call was successful
        for(int i=0; i<N; ++i) { // row
            for(int j=0; j<N; ++j) { // col
                EXPECT_NEAR(A_input_rm[get_2d_row_major_idx(i,j,N)], A_expected_data[get_2d_col_major_idx(i,j,N)], check_tol) << "A_rm mismatch at (" << i << "," << j << ")";
                EXPECT_NEAR(V_output_rm[get_2d_row_major_idx(i,j,N)], V_expected_data[get_2d_col_major_idx(i,j,N)], check_tol) << "V_rm mismatch at (" << i << "," << j << ")";
                EXPECT_NEAR(Y_output_rm[get_2d_row_major_idx(i,j,N)], Y_expected_data[get_2d_col_major_idx(i,j,N)], check_tol) << "Y_rm mismatch at (" << i << "," << j << ")";
            }
        }
        for(int i=0; i<N; ++i) {
            EXPECT_NEAR(VALR_output_rm[i], VALR_expected_data[i], check_tol) << "VALR_rm mismatch at " << i;
            EXPECT_NEAR(VALI_output_rm[i], VALI_expected_data[i], check_tol) << "VALI_rm mismatch at " << i;
        }
    }
}


TEST_F(MB05MDTest, Functionality_N1_BALANC_N) {
    InitializeData(1, 'N', 2.0); 
    if (N > 0 && !A_data.empty()) A_data[0] = 3.0; // A = [[3.0]]
    
    double expected_expA = 0.0;
    if (N > 0) expected_expA = std::exp(6.0); // exp(A*delta) = exp(3.0 * 2.0)

    info_result = slicot_mb05md(BALANC, N, DELTA, A_data.data(), LDA, 
                                V_data.data(), LDV, Y_data.data(), LDY, 
                                VALR_data.data(), VALI_data.data(), 0); // Col-major
    EXPECT_EQ(info_result, 0);

    if (info_result == 0 && N==1 && !A_data.empty() && !VALR_data.empty() && !VALI_data.empty() && !V_data.empty() && !Y_data.empty()) {
        EXPECT_NEAR(A_data[0], expected_expA, check_tol);
        EXPECT_NEAR(VALR_data[0], 3.0, check_tol); 
        EXPECT_NEAR(VALI_data[0], 0.0, check_tol);
        // Eigenvector for N=1 can be 1.0 or -1.0.
        // If V[0] is -1.0, then Y[0] would be -exp(6.0) for V*Y = exp(6.0)
        // The SLICOT example shows positive V.
        EXPECT_NEAR(std::abs(V_data[0]), 1.0, check_tol); 
        if (std::abs(V_data[0] - 1.0) < check_tol) {
            EXPECT_NEAR(Y_data[0], expected_expA, check_tol); 
        } else if (std::abs(V_data[0] + 1.0) < check_tol) {
            EXPECT_NEAR(Y_data[0], -expected_expA, check_tol);
        }
    }
}

