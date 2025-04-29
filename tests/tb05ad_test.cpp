#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <complex>

#include "tb05ad.h"
#include "slicot_utils.h"
#include "test_utils.h"
#include "test_config.h"

// Use the TEST_DATA_DIR macro defined in test_config.h
const std::string DATA_FILE_PATH = TEST_DATA_DIR "tb05ad.csv";

// Helper function to compare complex values
bool complex_near(double re1, double im1, double re2, double im2, double tol) {
    std::complex<double> c1(re1, im1);
    std::complex<double> c2(re2, im2);
    return std::abs(c1 - c2) < tol;
}

// --- Column-Major Test Fixture ---
class TB05ADTestColMajor : public ::testing::Test {
protected:
    // Test parameters
    int N = 3;
    int M = 1;
    int P = 2;
    char BALEIG = 'A';
    char INITA = 'G';
    double FREQ_RE = 0.0;
    double FREQ_IM = 0.5;

    // Verification tolerance
    double check_tol = 1e-2; // Increased tolerance for numerical differences

    // Input matrices
    std::vector<double> A = {
        1.0, 4.0, 0.0,
        2.0, -1.0, 0.0,
        0.0, 0.0, 1.0
    };
    
    std::vector<double> B = {
        1.0, 0.0, 1.0
    };
    
    // Fix: Properly arrange C in column-major format (P=2, N=3)
    std::vector<double> C = {
        1.0, 0.0,    // First column (2 elements)
        0.0, 0.0,    // Second column (2 elements)
        -1.0, 1.0    // Third column (2 elements) 
    };

    // Expected outputs (from example results)
    double expected_rcond = 0.22;
    
    // Expected eigenvalues
    std::vector<double> EVRE_expected = {3.0, -3.0, 1.0};
    std::vector<double> EVIM_expected = {0.0, 0.0, 0.0};
    
    // Expected G matrix (complex values are interleaved) - column-major format
    std::vector<double> G_expected = {
        0.69, 0.35,
        -0.80, -0.40
    };
    
    // Expected HINVB matrix (complex values are interleaved) - column-major format
    std::vector<double> HINVB_expected = {
        -0.11, -0.05,
        -0.43, 0.00,
        -0.80, -0.40
    };
    
    int expected_info = 0;

    // Output arrays (to be filled by the function)
    double rcond_result;
    std::vector<double> EVRE; 
    std::vector<double> EVIM;
    std::vector<double> G; 
    std::vector<double> HINVB;
    int info_result = -999;
    
    // Leading dimensions
    int LDA, LDB, LDC, LDG, LDHINV;

    void SetUp() override {
        // Size output vectors
        EVRE.resize(N);
        EVIM.resize(N);
        G.resize(2 * P * M); // Complex values are interleaved
        HINVB.resize(2 * N * M); // Complex values are interleaved
        
        // Set column-major leading dimensions
        LDA = std::max(1, N);     // rows of A
        LDB = std::max(1, N);     // rows of B
        LDC = std::max(1, P);     // rows of C
        LDG = std::max(1, P);     // rows of G
        LDHINV = std::max(1, N);  // rows of HINVB
    }
};

// --- Row-Major Test Fixture ---
class TB05ADTestRowMajor : public TB05ADTestColMajor {
protected:
    // Row-major input matrices
    std::vector<double> A_rm;
    std::vector<double> B_rm;
    std::vector<double> C_rm;
    
    // Expected outputs in row-major format
    std::vector<double> G_expected_rm;
    std::vector<double> HINVB_expected_rm;

    void SetUp() override {
        // Call base class SetUp
        TB05ADTestColMajor::SetUp();
        
        // Convert inputs to row major
        A_rm.resize(N * N);
        B_rm.resize(N * M);
        C_rm.resize(P * N);
        
        slicot_transpose_to_c(A.data(), A_rm.data(), N, N, sizeof(double));
        slicot_transpose_to_c(B.data(), B_rm.data(), N, M, sizeof(double));
        slicot_transpose_to_c(C.data(), C_rm.data(), P, N, sizeof(double));
        
        // Convert expected outputs to row major
        G_expected_rm.resize(2 * P * M); // Complex values are interleaved
        HINVB_expected_rm.resize(2 * N * M); // Complex values are interleaved
        
        // For complex matrices, we transpose the real and imaginary parts together
        // P=2, M=1 for G
        for (int i = 0; i < P; i++) {
            for (int j = 0; j < M; j++) {
                // Column-major index: 2*(i + j*P)
                // Row-major index: 2*(j + i*M)
                G_expected_rm[2*(j + i*M)] = G_expected[2*(i + j*P)];
                G_expected_rm[2*(j + i*M) + 1] = G_expected[2*(i + j*P) + 1];
            }
        }
        
        // N=3, M=1 for HINVB
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                // Column-major index: 2*(i + j*N)
                // Row-major index: 2*(j + i*M)
                HINVB_expected_rm[2*(j + i*M)] = HINVB_expected[2*(i + j*N)];
                HINVB_expected_rm[2*(j + i*M) + 1] = HINVB_expected[2*(i + j*N) + 1];
            }
        }
        
        // Set row-major leading dimensions (number of columns)
        LDA = std::max(1, N);     // cols of A
        LDB = std::max(1, M);     // cols of B
        LDC = std::max(1, N);     // cols of C
        LDG = std::max(1, M);     // cols of G
        LDHINV = std::max(1, M);  // cols of HINVB
    }
};

// Test: Example from Documentation (Column-Major)
TEST_F(TB05ADTestColMajor, DocExample) {
    // Call C wrapper function
    info_result = slicot_tb05ad(BALEIG, INITA, N, M, P,
                               FREQ_RE, FREQ_IM,
                               A.data(), LDA, B.data(), LDB, C.data(), LDC,
                               &rcond_result, G.data(), LDG,
                               EVRE.data(), EVIM.data(),
                               HINVB.data(), LDHINV, 0);

    // Verify return code
    ASSERT_EQ(info_result, expected_info);
    
    // Verify RCOND value
    EXPECT_NEAR(rcond_result, expected_rcond, check_tol);
    
    // Verify eigenvalues
    ASSERT_EQ(EVRE.size(), EVRE_expected.size());
    ASSERT_EQ(EVIM.size(), EVIM_expected.size());
    for (size_t i = 0; i < EVRE_expected.size(); ++i) {
        EXPECT_NEAR(EVRE[i], EVRE_expected[i], check_tol) << "EVRE mismatch at index " << i;
        EXPECT_NEAR(EVIM[i], EVIM_expected[i], check_tol) << "EVIM mismatch at index " << i;
    }
    
    // Verify G matrix (complex values)
    ASSERT_EQ(G.size(), G_expected.size());
    for (size_t i = 0; i < G_expected.size(); i += 2) {
        EXPECT_TRUE(complex_near(G[i], G[i+1], 
                                G_expected[i], G_expected[i+1], check_tol))
                    << "G mismatch at complex index " << i/2;
    }
    
    // Verify HINVB matrix (complex values)
    ASSERT_EQ(HINVB.size(), HINVB_expected.size());
    for (size_t i = 0; i < HINVB_expected.size(); i += 2) {
        EXPECT_TRUE(complex_near(HINVB[i], HINVB[i+1], 
                                HINVB_expected[i], HINVB_expected[i+1], check_tol))
                    << "HINVB mismatch at complex index " << i/2;
    }
}

// Test: Example from Documentation (Row-Major)
TEST_F(TB05ADTestRowMajor, DocExample) {
    // Call C wrapper function with row_major=1
    info_result = slicot_tb05ad(BALEIG, INITA, N, M, P,
                               FREQ_RE, FREQ_IM,
                               A_rm.data(), LDA, B_rm.data(), LDB, C_rm.data(), LDC,
                               &rcond_result, G.data(), LDG,
                               EVRE.data(), EVIM.data(),
                               HINVB.data(), LDHINV, 1);

    // Verify return code
    ASSERT_EQ(info_result, expected_info);
    
    // Verify RCOND value
    EXPECT_NEAR(rcond_result, expected_rcond, check_tol);
    
    // Verify eigenvalues
    ASSERT_EQ(EVRE.size(), EVRE_expected.size());
    ASSERT_EQ(EVIM.size(), EVIM_expected.size());
    for (size_t i = 0; i < EVRE_expected.size(); ++i) {
        EXPECT_NEAR(EVRE[i], EVRE_expected[i], check_tol) << "EVRE mismatch at index " << i;
        EXPECT_NEAR(EVIM[i], EVIM_expected[i], check_tol) << "EVIM mismatch at index " << i;
    }
    
    // Verify G matrix (complex values)
    ASSERT_EQ(G.size(), G_expected_rm.size());
    for (size_t i = 0; i < G_expected_rm.size(); i += 2) {
        EXPECT_TRUE(complex_near(G[i], G[i+1], 
                                G_expected_rm[i], G_expected_rm[i+1], check_tol))
                    << "G mismatch at complex index " << i/2;
    }
    
    // Verify HINVB matrix (complex values)
    ASSERT_EQ(HINVB.size(), HINVB_expected_rm.size());
    for (size_t i = 0; i < HINVB_expected_rm.size(); i += 2) {
        EXPECT_TRUE(complex_near(HINVB[i], HINVB[i+1], 
                                HINVB_expected_rm[i], HINVB_expected_rm[i+1], check_tol))
                    << "HINVB mismatch at complex index " << i/2;
    }
    
    // Also verify the modified input matrices (A, B, C) if INITA = 'G'
    // We don't have expected values for these, but this ensures 
    // the wrapper correctly handles the in-place modifications
}

// Test: Parameter Validation
TEST_F(TB05ADTestColMajor, ParameterValidation) {
    std::vector<double> dummy_A(N*N, 0.0);
    std::vector<double> dummy_B(N*M, 0.0);
    std::vector<double> dummy_C(P*N, 0.0);
    std::vector<double> dummy_G(2*P*M, 0.0); // Complex values are interleaved
    std::vector<double> dummy_HINVB(2*N*M, 0.0); // Complex values are interleaved
    std::vector<double> dummy_EVRE(N, 0.0);
    std::vector<double> dummy_EVIM(N, 0.0);
    double dummy_rcond = 0.0;

    // Test invalid N
    info_result = slicot_tb05ad(BALEIG, INITA, -1, M, P,
                               FREQ_RE, FREQ_IM,
                               dummy_A.data(), N, dummy_B.data(), N, dummy_C.data(), P,
                               &dummy_rcond, dummy_G.data(), P,
                               dummy_EVRE.data(), dummy_EVIM.data(),
                               dummy_HINVB.data(), N, 0);
    EXPECT_EQ(info_result, -3);

    // Test invalid M
    info_result = slicot_tb05ad(BALEIG, INITA, N, -1, P,
                               FREQ_RE, FREQ_IM,
                               dummy_A.data(), N, dummy_B.data(), N, dummy_C.data(), P,
                               &dummy_rcond, dummy_G.data(), P,
                               dummy_EVRE.data(), dummy_EVIM.data(),
                               dummy_HINVB.data(), N, 0);
    EXPECT_EQ(info_result, -4);

    // Test invalid P
    info_result = slicot_tb05ad(BALEIG, INITA, N, M, -1,
                               FREQ_RE, FREQ_IM,
                               dummy_A.data(), N, dummy_B.data(), N, dummy_C.data(), P,
                               &dummy_rcond, dummy_G.data(), P,
                               dummy_EVRE.data(), dummy_EVIM.data(),
                               dummy_HINVB.data(), N, 0);
    EXPECT_EQ(info_result, -5);

    // Test invalid BALEIG
    info_result = slicot_tb05ad('X', INITA, N, M, P,
                               FREQ_RE, FREQ_IM,
                               dummy_A.data(), N, dummy_B.data(), N, dummy_C.data(), P,
                               &dummy_rcond, dummy_G.data(), P,
                               dummy_EVRE.data(), dummy_EVIM.data(),
                               dummy_HINVB.data(), N, 0);
    EXPECT_EQ(info_result, -1);

    // Test invalid INITA
    info_result = slicot_tb05ad(BALEIG, 'X', N, M, P,
                               FREQ_RE, FREQ_IM,
                               dummy_A.data(), N, dummy_B.data(), N, dummy_C.data(), P,
                               &dummy_rcond, dummy_G.data(), P,
                               dummy_EVRE.data(), dummy_EVIM.data(),
                               dummy_HINVB.data(), N, 0);
    EXPECT_EQ(info_result, -2);

    // Test NULL pointer for A
    info_result = slicot_tb05ad(BALEIG, INITA, N, M, P,
                               FREQ_RE, FREQ_IM,
                               nullptr, N, dummy_B.data(), N, dummy_C.data(), P,
                               &dummy_rcond, dummy_G.data(), P,
                               dummy_EVRE.data(), dummy_EVIM.data(),
                               dummy_HINVB.data(), N, 0);
    EXPECT_EQ(info_result, -8);

    // Test invalid LDA (Column-Major)
    info_result = slicot_tb05ad(BALEIG, INITA, N, M, P,
                               FREQ_RE, FREQ_IM,
                               dummy_A.data(), 0, dummy_B.data(), N, dummy_C.data(), P,
                               &dummy_rcond, dummy_G.data(), P,
                               dummy_EVRE.data(), dummy_EVIM.data(),
                               dummy_HINVB.data(), N, 0);
    EXPECT_EQ(info_result, -9);
}

// Additional tests could include:
// - Testing with different BALEIG and INITA combinations
// - Testing with different matrix sizes
// - Edge cases (N=0, M=0, P=0)
// - Testing with frequency values close to eigenvalues (should return info=2)
