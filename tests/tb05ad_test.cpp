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

// Test: 4x4 Known System from Python test
TEST_F(TB05ADTestColMajor, KnownSystem4x4) {
    // Define a 4x4 system from the Python test
    int N_4x4 = 4;
    int M_4x4 = 2;
    int P_4x4 = 3;
    
    // Define matrices (column-major format)
    std::vector<double> A_4x4 = {
        -0.5, 0.0, 1.0, 0.0,
        0.0, -1.0, 0.0, 1.0,
        0.0, 0.0, -0.5, 0.0,
        0.0, 0.0, 0.0, -1.0
    };
    
    std::vector<double> B_4x4 = {
        1.0, 0.0, 0.0, 0.0, 
        0.0, 1.0, 0.0, 0.0
    };
    
    std::vector<double> C_4x4 = {
        0.0, 0.0, 0.0, 
        1.0, 1.0, 1.0,
        1.0, 0.0, 1.0,
        0.0, 1.0, 1.0
    };
    
    // Output arrays
    double rcond_result;
    std::vector<double> EVRE(N_4x4);
    std::vector<double> EVIM(N_4x4);
    std::vector<double> G(2 * P_4x4 * M_4x4);
    std::vector<double> HINVB(2 * N_4x4 * M_4x4);
    
    // Set frequency to 10j (similar to Python test)
    double FREQ_RE = 0.0;
    double FREQ_IM = 10.0;
    
    // Set parameters
    char BALEIG = 'A';  // corresponds to job='AG' in Python
    char INITA = 'G';
    
    // Leading dimensions
    int LDA = N_4x4;
    int LDB = N_4x4;
    int LDC = P_4x4;
    int LDG = P_4x4;
    int LDHINV = N_4x4;
    
    // Call function
    int info_result = slicot_tb05ad(BALEIG, INITA, N_4x4, M_4x4, P_4x4,
                                  FREQ_RE, FREQ_IM,
                                  A_4x4.data(), LDA, B_4x4.data(), LDB, C_4x4.data(), LDC,
                                  &rcond_result, G.data(), LDG,
                                  EVRE.data(), EVIM.data(),
                                  HINVB.data(), LDHINV, 0);
    
    // Check return code
    EXPECT_EQ(info_result, 0);
    
    // Expected eigenvalues should be around -0.5, -0.5, -1, -1
    std::vector<double> expected_eigenvalues = {-1.0, -1.0, -0.5, -0.5};
    
    // Sort the computed eigenvalues for comparison
    std::vector<std::pair<double, double>> eig_pairs;
    for (int i = 0; i < N_4x4; i++) {
        eig_pairs.push_back({EVRE[i], EVIM[i]});
    }
    std::sort(eig_pairs.begin(), eig_pairs.end());
    std::vector<double> sorted_eigs;
    for (const auto& pair : eig_pairs) {
        sorted_eigs.push_back(pair.first);
    }
    
    // Verify eigenvalues
    for (int i = 0; i < N_4x4; i++) {
        EXPECT_NEAR(sorted_eigs[i], expected_eigenvalues[i], check_tol);
    }
}

// Test: Different Job Parameters
TEST_F(TB05ADTestColMajor, DifferentJobParameters) {
    // Test with BALEIG = 'N', INITA = 'G' (corresponds to job='NG' in Python)
    char baleig_ng = 'N';
    char inita_g = 'G';
    
    // Output arrays for 'NG' job
    double rcond_ng;
    std::vector<double> EVRE_ng(N); 
    std::vector<double> EVIM_ng(N);
    std::vector<double> G_ng(2 * P * M);
    std::vector<double> HINVB_ng(2 * N * M);
    
    // Create copies of the original matrices for the NG job
    std::vector<double> A_ng = A;
    std::vector<double> B_ng = B;
    std::vector<double> C_ng = C;
    
    int info_ng = slicot_tb05ad(baleig_ng, inita_g, N, M, P,
                              FREQ_RE, FREQ_IM,
                              A_ng.data(), LDA, B_ng.data(), LDB, C_ng.data(), LDC,
                              &rcond_ng, G_ng.data(), LDG,
                              EVRE_ng.data(), EVIM_ng.data(),
                              HINVB_ng.data(), LDHINV, 0);
    
    EXPECT_EQ(info_ng, 0);
    
    // Test with BALEIG = 'A', INITA = 'G' (corresponds to job='AG' in Python)
    char baleig_ag = 'A';
    char inita_ag = 'G';
    
    // Output arrays for 'AG' job
    double rcond_ag;
    std::vector<double> EVRE_ag(N); 
    std::vector<double> EVIM_ag(N);
    std::vector<double> G_ag(2 * P * M);
    std::vector<double> HINVB_ag(2 * N * M);
    
    // Create copies of the original matrices for the AG job
    std::vector<double> A_ag = A;
    std::vector<double> B_ag = B;
    std::vector<double> C_ag = C;
    
    int info_ag = slicot_tb05ad(baleig_ag, inita_ag, N, M, P,
                              FREQ_RE, FREQ_IM,
                              A_ag.data(), LDA, B_ag.data(), LDB, C_ag.data(), LDC,
                              &rcond_ag, G_ag.data(), LDG,
                              EVRE_ag.data(), EVIM_ag.data(),
                              HINVB_ag.data(), LDHINV, 0);
    
    EXPECT_EQ(info_ag, 0);
    
    // Both methods should produce the same G matrix since they're computing
    // the same transfer function
    for (size_t i = 0; i < G_ng.size(); i += 2) {
        EXPECT_TRUE(complex_near(G_ng[i], G_ng[i+1], 
                                G_ag[i], G_ag[i+1], check_tol))
                    << "G matrix mismatch between NG and AG jobs at index " << i/2;
    }
    
    // NOTE: We don't compare eigenvalues because BALEIG='N' doesn't compute them
    // or computes them differently than BALEIG='A'
}

// Test: Balancing
TEST_F(TB05ADTestColMajor, Balancing) {
    // Skip actual balancing comparison and just verify the base class test works
    // This avoids complicating with a special test case
    
    std::vector<double> A_copy = A;
    std::vector<double> B_copy = B;
    std::vector<double> C_copy = C;
    
    // Use the original test fixture data and parameters
    double rcond_result;
    std::vector<double> EVRE_bal(N);
    std::vector<double> EVIM_bal(N);
    std::vector<double> G_balanced(2 * P * M);
    std::vector<double> HINVB_bal(2 * N * M);
    
    // Call function with balancing (BALEIG='A')
    int info_balanced = slicot_tb05ad('A', INITA, N, M, P,
                                     FREQ_RE, FREQ_IM,
                                     A_copy.data(), LDA, B_copy.data(), LDB, C_copy.data(), LDC,
                                     &rcond_result, G_balanced.data(), LDG,
                                     EVRE_bal.data(), EVIM_bal.data(),
                                     HINVB_bal.data(), LDHINV, 0);
    
    // Verify the function executes without error
    EXPECT_EQ(info_balanced, 0) << "Balancing with BALEIG='A' should succeed";
    
    // Make a fresh copy of data for the non-balanced test
    A_copy = A;
    B_copy = B;
    C_copy = C;
    
    // Run without balancing (BALEIG='N')
    double rcond_nobal;
    std::vector<double> EVRE_nobal(N);
    std::vector<double> EVIM_nobal(N);
    std::vector<double> G_nobalanced(2 * P * M);
    std::vector<double> HINVB_nobal(2 * N * M);
    
    // Call function without balancing (BALEIG='N')
    int info_nobalanced = slicot_tb05ad('N', INITA, N, M, P,
                                       FREQ_RE, FREQ_IM,
                                       A_copy.data(), LDA, B_copy.data(), LDB, C_copy.data(), LDC,
                                       &rcond_nobal, G_nobalanced.data(), LDG,
                                       EVRE_nobal.data(), EVIM_nobal.data(),
                                       HINVB_nobal.data(), LDHINV, 0);
    
    // Verify the function executes without error
    EXPECT_EQ(info_nobalanced, 0) << "Execution without balancing (BALEIG='N') should succeed";
    
    // If both executions succeeded, verify the G matrices are similar
    if (info_balanced == 0 && info_nobalanced == 0) {
        for (size_t i = 0; i < G_balanced.size(); i += 2) {
            EXPECT_TRUE(complex_near(G_balanced[i], G_balanced[i+1],
                                    G_nobalanced[i], G_nobalanced[i+1], check_tol))
                << "G matrix mismatch between balanced and non-balanced at index " << i/2;
        }
    }
}

// Test: Complex Frequency
TEST_F(TB05ADTestColMajor, ComplexFrequency) {
    // Set a complex frequency (both real and imaginary parts non-zero)
    double FREQ_RE_complex = 5.0;
    double FREQ_IM_complex = 8.0;
    
    // Output arrays
    double rcond_result;
    std::vector<double> EVRE(N);
    std::vector<double> EVIM(N);
    std::vector<double> G(2 * P * M);
    std::vector<double> HINVB(2 * N * M);
    
    // Make a copy of original matrices to avoid modifying test fixture data
    std::vector<double> A_copy = A;
    std::vector<double> B_copy = B;
    std::vector<double> C_copy = C;
    
    int info_result = slicot_tb05ad(BALEIG, INITA, N, M, P,
                                  FREQ_RE_complex, FREQ_IM_complex,
                                  A_copy.data(), LDA, B_copy.data(), LDB, C_copy.data(), LDC,
                                  &rcond_result, G.data(), LDG,
                                  EVRE.data(), EVIM.data(),
                                  HINVB.data(), LDHINV, 0);
    
    EXPECT_EQ(info_result, 0);
}

// Test: Resonance Condition
TEST_F(TB05ADTestColMajor, ResonanceCondition) {
    // Create a system with known eigenvalues
    std::vector<double> A_resonance = {
        0.0, -1.0,
        1.0, 0.0
    };
    
    std::vector<double> B_resonance = {
        1.0, 0.0
    };
    
    std::vector<double> C_resonance = {
        0.0, 1.0
    };
    
    // Set frequency close to an eigenvalue
    double FREQ_RE_res = 0.0;
    double FREQ_IM_res = 1.0; // The eigenvalues are ±i
    
    // Output arrays
    double rcond_res;
    std::vector<double> EVRE_res(2);
    std::vector<double> EVIM_res(2);
    std::vector<double> G_res(2); // 1x1 complex number
    std::vector<double> HINVB_res(4); // 2x1 complex vector
    
    int info_res = slicot_tb05ad('N', 'G', 2, 1, 1,
                               FREQ_RE_res, FREQ_IM_res,
                               A_resonance.data(), 2, B_resonance.data(), 2, C_resonance.data(), 1,
                               &rcond_res, G_res.data(), 1,
                               EVRE_res.data(), EVIM_res.data(),
                               HINVB_res.data(), 2, 0);
    
    // We expect info = 2 for resonance condition (frequency too close to eigenvalue)
    EXPECT_EQ(info_res, 2);
}
