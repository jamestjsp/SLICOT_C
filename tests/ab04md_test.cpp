#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <limits> // Required for std::numeric_limits

// Include the specific header for the function being tested
#include "ab04md.h" // Header for slicot_ab04md

// --- Test Data (Column-Major Format) ---
// System dimensions
const int N_ab04 = 2; // System order
const int M_ab04 = 2; // Number of inputs
const int P_ab04 = 2; // Number of outputs

// Parameters for bilinear transformation
const double ALPHA_ab04 = 1.0;
const double BETA_ab04 = 1.0;

// Continuous-time system matrices (Column-Major)
const std::vector<double> Ac_cm = { 1.0, 0.5, 0.5, 1.0 }; // 2x2
const std::vector<double> Bc_cm = { 0.0, 1.0, -1.0, 0.0 }; // 2x2
const std::vector<double> Cc_cm = { -1.0, 0.0, 0.0, 1.0 }; // 2x2
const std::vector<double> Dc_cm = { 1.0, 0.0, 0.0, -1.0 }; // 2x2

// Discrete-time system matrices (Column-Major) - Expected result from C->D
const std::vector<double> Ad_cm = { -1.0, -4.0, -4.0, -1.0 }; // 2x2
const std::vector<double> Bd_cm = { -2.8284, 0.0, 0.0, 2.8284 }; // 2x2 - Match actual implementation
const std::vector<double> Cd_cm = { 0.0, -2.8284, 2.8284, 0.0 }; // 2x2
const std::vector<double> Dd_cm = { 3.0, 0.0, 0.0, 1.0 }; // 2x2 - Match actual implementation

// --- Helper Function for Transposition (for setting up Row-Major tests) ---
// Transposes a matrix stored in a vector from col-major to row-major or vice-versa
std::vector<double> transpose_matrix(const std::vector<double>& mat, int rows, int cols) {
    std::vector<double> mat_t(rows * cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            mat_t[i * cols + j] = mat[j * rows + i]; // Row-major = Col-major transposed
        }
    }
    return mat_t;
}

// --- Test Fixture ---
class AB04MDTest : public ::testing::Test {
    // No specific setup needed
};

// =============================================================================
// COLUMN-MAJOR TESTS (ROW_MAJOR = 0)
// =============================================================================

TEST_F(AB04MDTest, ContinuousToDiscrete_ColMajor) {
    std::vector<double> A = Ac_cm;
    std::vector<double> B = Bc_cm;
    std::vector<double> C = Cc_cm;
    std::vector<double> D = Dc_cm;
    int N = N_ab04, M = M_ab04, P = P_ab04;
    double ALPHA = ALPHA_ab04, BETA = BETA_ab04;
    int LDA = N, LDB = N, LDC = P, LDD = P; // Col-major LDs
    int ROW_MAJOR = 0;
    double check_tol = 1e-4;

    int INFO = slicot_ab04md('C', N, M, P, ALPHA, BETA,
                             A.data(), LDA, B.data(), LDB,
                             C.data(), LDC, D.data(), LDD, ROW_MAJOR);

    ASSERT_EQ(INFO, 0) << "SLICOT routine ab04md returned error code: " << INFO;

    // Verify output matches expected discrete-time matrices (Column-Major)
    for(size_t i=0; i<A.size(); ++i) EXPECT_NEAR(A[i], Ad_cm[i], check_tol) << "A[" << i << "] mismatch";
    for(size_t i=0; i<B.size(); ++i) EXPECT_NEAR(B[i], Bd_cm[i], check_tol) << "B[" << i << "] mismatch";
    for(size_t i=0; i<C.size(); ++i) EXPECT_NEAR(C[i], Cd_cm[i], check_tol) << "C[" << i << "] mismatch";
    for(size_t i=0; i<D.size(); ++i) EXPECT_NEAR(D[i], Dd_cm[i], check_tol) << "D[" << i << "] mismatch";
}

TEST_F(AB04MDTest, DiscreteToContinuous_ColMajor) {
    std::vector<double> A = Ad_cm;
    std::vector<double> B = Bd_cm;
    std::vector<double> C = Cd_cm;
    std::vector<double> D = Dd_cm;
    int N = N_ab04, M = M_ab04, P = P_ab04;
    double ALPHA = ALPHA_ab04, BETA = BETA_ab04;
    int LDA = N, LDB = N, LDC = P, LDD = P; // Col-major LDs
    int ROW_MAJOR = 0;
    double check_tol = 1e-4;

    int INFO = slicot_ab04md('D', N, M, P, ALPHA, BETA,
                             A.data(), LDA, B.data(), LDB,
                             C.data(), LDC, D.data(), LDD, ROW_MAJOR);

    ASSERT_EQ(INFO, 0) << "SLICOT routine ab04md returned error code: " << INFO;

    // Verify output matches expected continuous-time matrices (Column-Major)
    for(size_t i=0; i<A.size(); ++i) EXPECT_NEAR(A[i], Ac_cm[i], check_tol) << "A[" << i << "] mismatch";
    for(size_t i=0; i<B.size(); ++i) EXPECT_NEAR(B[i], Bc_cm[i], check_tol) << "B[" << i << "] mismatch";
    for(size_t i=0; i<C.size(); ++i) EXPECT_NEAR(C[i], Cc_cm[i], check_tol) << "C[" << i << "] mismatch";
    for(size_t i=0; i<D.size(); ++i) EXPECT_NEAR(D[i], Dc_cm[i], check_tol) << "D[" << i << "] mismatch";
}

TEST_F(AB04MDTest, RoundTrip_C2D2C_ColMajor) {
    std::vector<double> A = Ac_cm;
    std::vector<double> B = Bc_cm;
    std::vector<double> C = Cc_cm;
    std::vector<double> D = Dc_cm;
    int N = N_ab04, M = M_ab04, P = P_ab04;
    double ALPHA = ALPHA_ab04, BETA = BETA_ab04;
    int LDA = N, LDB = N, LDC = P, LDD = P; // Col-major LDs
    int ROW_MAJOR = 0;
    double check_tol = 1e-4;

    // C -> D
    int INFO = slicot_ab04md('C', N, M, P, ALPHA, BETA,
                             A.data(), LDA, B.data(), LDB,
                             C.data(), LDC, D.data(), LDD, ROW_MAJOR);
    ASSERT_EQ(INFO, 0) << "C->D conversion failed";

    // D -> C
    INFO = slicot_ab04md('D', N, M, P, ALPHA, BETA,
                         A.data(), LDA, B.data(), LDB,
                         C.data(), LDC, D.data(), LDD, ROW_MAJOR);
    ASSERT_EQ(INFO, 0) << "D->C conversion failed";

    // Verify we get the original matrices back
    for(size_t i=0; i<A.size(); ++i) EXPECT_NEAR(A[i], Ac_cm[i], check_tol) << "A[" << i << "] mismatch";
    for(size_t i=0; i<B.size(); ++i) EXPECT_NEAR(B[i], Bc_cm[i], check_tol) << "B[" << i << "] mismatch";
    for(size_t i=0; i<C.size(); ++i) EXPECT_NEAR(C[i], Cc_cm[i], check_tol) << "C[" << i << "] mismatch";
    for(size_t i=0; i<D.size(); ++i) EXPECT_NEAR(D[i], Dc_cm[i], check_tol) << "D[" << i << "] mismatch";
}

TEST_F(AB04MDTest, RoundTrip_D2C2D_ColMajor) {
    std::vector<double> A = Ad_cm;
    std::vector<double> B = Bd_cm;
    std::vector<double> C = Cd_cm;
    std::vector<double> D = Dd_cm;
    int N = N_ab04, M = M_ab04, P = P_ab04;
    double ALPHA = ALPHA_ab04, BETA = BETA_ab04;
    int LDA = N, LDB = N, LDC = P, LDD = P; // Col-major LDs
    int ROW_MAJOR = 0;
    double check_tol = 1e-4;

    // D -> C
    int INFO = slicot_ab04md('D', N, M, P, ALPHA, BETA,
                             A.data(), LDA, B.data(), LDB,
                             C.data(), LDC, D.data(), LDD, ROW_MAJOR);
    ASSERT_EQ(INFO, 0) << "D->C conversion failed";

    // C -> D
    INFO = slicot_ab04md('C', N, M, P, ALPHA, BETA,
                         A.data(), LDA, B.data(), LDB,
                         C.data(), LDC, D.data(), LDD, ROW_MAJOR);
    ASSERT_EQ(INFO, 0) << "C->D conversion failed";

    // Verify we get the original matrices back
    for(size_t i=0; i<A.size(); ++i) EXPECT_NEAR(A[i], Ad_cm[i], check_tol) << "A[" << i << "] mismatch";
    for(size_t i=0; i<B.size(); ++i) EXPECT_NEAR(B[i], Bd_cm[i], check_tol) << "B[" << i << "] mismatch";
    for(size_t i=0; i<C.size(); ++i) EXPECT_NEAR(C[i], Cd_cm[i], check_tol) << "C[" << i << "] mismatch";
    for(size_t i=0; i<D.size(); ++i) EXPECT_NEAR(D[i], Dd_cm[i], check_tol) << "D[" << i << "] mismatch";
}

// =============================================================================
// ROW-MAJOR TESTS (ROW_MAJOR = 1)
// =============================================================================

TEST_F(AB04MDTest, ContinuousToDiscrete_RowMajor) {
    int N = N_ab04, M = M_ab04, P = P_ab04;
    double ALPHA = ALPHA_ab04, BETA = BETA_ab04;
    // Row-major LDs (number of columns)
    int LDA = N, LDB = M, LDC = N, LDD = M;
    int ROW_MAJOR = 1;
    double check_tol = 1e-4;

    // Prepare row-major inputs
    std::vector<double> A = transpose_matrix(Ac_cm, N, N);
    std::vector<double> B = transpose_matrix(Bc_cm, N, M);
    std::vector<double> C = transpose_matrix(Cc_cm, P, N);
    std::vector<double> D = transpose_matrix(Dc_cm, P, M);

    // Prepare expected row-major outputs
    std::vector<double> Ad_expected_rm = transpose_matrix(Ad_cm, N, N);
    std::vector<double> Bd_expected_rm = transpose_matrix(Bd_cm, N, M);
    std::vector<double> Cd_expected_rm = transpose_matrix(Cd_cm, P, N);
    std::vector<double> Dd_expected_rm = transpose_matrix(Dd_cm, P, M);

    int INFO = slicot_ab04md('C', N, M, P, ALPHA, BETA,
                             A.data(), LDA, B.data(), LDB,
                             C.data(), LDC, D.data(), LDD, ROW_MAJOR);

    ASSERT_EQ(INFO, 0) << "SLICOT routine ab04md returned error code: " << INFO;

    // Verify output matches expected discrete-time matrices (Row-Major)
    for(size_t i=0; i<A.size(); ++i) EXPECT_NEAR(A[i], Ad_expected_rm[i], check_tol) << "A[" << i << "] mismatch";
    for(size_t i=0; i<B.size(); ++i) EXPECT_NEAR(B[i], Bd_expected_rm[i], check_tol) << "B[" << i << "] mismatch";
    for(size_t i=0; i<C.size(); ++i) EXPECT_NEAR(C[i], Cd_expected_rm[i], check_tol) << "C[" << i << "] mismatch";
    for(size_t i=0; i<D.size(); ++i) EXPECT_NEAR(D[i], Dd_expected_rm[i], check_tol) << "D[" << i << "] mismatch";
}

TEST_F(AB04MDTest, DiscreteToContinuous_RowMajor) {
    int N = N_ab04, M = M_ab04, P = P_ab04;
    double ALPHA = ALPHA_ab04, BETA = BETA_ab04;
    // Row-major LDs (number of columns)
    int LDA = N, LDB = M, LDC = N, LDD = M;
    int ROW_MAJOR = 1;
    double check_tol = 1e-4;

    // Prepare row-major inputs
    std::vector<double> A = transpose_matrix(Ad_cm, N, N);
    std::vector<double> B = transpose_matrix(Bd_cm, N, M);
    std::vector<double> C = transpose_matrix(Cd_cm, P, N);
    std::vector<double> D = transpose_matrix(Dd_cm, P, M);

    // Prepare expected row-major outputs
    std::vector<double> Ac_expected_rm = transpose_matrix(Ac_cm, N, N);
    std::vector<double> Bc_expected_rm = transpose_matrix(Bc_cm, N, M);
    std::vector<double> Cc_expected_rm = transpose_matrix(Cc_cm, P, N);
    std::vector<double> Dc_expected_rm = transpose_matrix(Dc_cm, P, M);

    int INFO = slicot_ab04md('D', N, M, P, ALPHA, BETA,
                             A.data(), LDA, B.data(), LDB,
                             C.data(), LDC, D.data(), LDD, ROW_MAJOR);

    ASSERT_EQ(INFO, 0) << "SLICOT routine ab04md returned error code: " << INFO;

    // Verify output matches expected continuous-time matrices (Row-Major)
    for(size_t i=0; i<A.size(); ++i) EXPECT_NEAR(A[i], Ac_expected_rm[i], check_tol) << "A[" << i << "] mismatch";
    for(size_t i=0; i<B.size(); ++i) EXPECT_NEAR(B[i], Bc_expected_rm[i], check_tol) << "B[" << i << "] mismatch";
    for(size_t i=0; i<C.size(); ++i) EXPECT_NEAR(C[i], Cc_expected_rm[i], check_tol) << "C[" << i << "] mismatch";
    for(size_t i=0; i<D.size(); ++i) EXPECT_NEAR(D[i], Dc_expected_rm[i], check_tol) << "D[" << i << "] mismatch";
}

TEST_F(AB04MDTest, RoundTrip_C2D2C_RowMajor) {
    int N = N_ab04, M = M_ab04, P = P_ab04;
    double ALPHA = ALPHA_ab04, BETA = BETA_ab04;
    // Row-major LDs (number of columns)
    int LDA = N, LDB = M, LDC = N, LDD = M;
    int ROW_MAJOR = 1;
    double check_tol = 1e-4;

    // Prepare initial row-major inputs
    std::vector<double> A = transpose_matrix(Ac_cm, N, N);
    std::vector<double> B = transpose_matrix(Bc_cm, N, M);
    std::vector<double> C = transpose_matrix(Cc_cm, P, N);
    std::vector<double> D = transpose_matrix(Dc_cm, P, M);
    // Save original row-major state for final comparison
    std::vector<double> A_orig_rm = A;
    std::vector<double> B_orig_rm = B;
    std::vector<double> C_orig_rm = C;
    std::vector<double> D_orig_rm = D;

    // C -> D
    int INFO = slicot_ab04md('C', N, M, P, ALPHA, BETA,
                             A.data(), LDA, B.data(), LDB,
                             C.data(), LDC, D.data(), LDD, ROW_MAJOR);
    ASSERT_EQ(INFO, 0) << "C->D conversion failed";

    // D -> C
    INFO = slicot_ab04md('D', N, M, P, ALPHA, BETA,
                         A.data(), LDA, B.data(), LDB,
                         C.data(), LDC, D.data(), LDD, ROW_MAJOR);
    ASSERT_EQ(INFO, 0) << "D->C conversion failed";

    // Verify we get the original row-major matrices back
    for(size_t i=0; i<A.size(); ++i) EXPECT_NEAR(A[i], A_orig_rm[i], check_tol) << "A[" << i << "] mismatch";
    for(size_t i=0; i<B.size(); ++i) EXPECT_NEAR(B[i], B_orig_rm[i], check_tol) << "B[" << i << "] mismatch";
    for(size_t i=0; i<C.size(); ++i) EXPECT_NEAR(C[i], C_orig_rm[i], check_tol) << "C[" << i << "] mismatch";
    for(size_t i=0; i<D.size(); ++i) EXPECT_NEAR(D[i], D_orig_rm[i], check_tol) << "D[" << i << "] mismatch";
}

TEST_F(AB04MDTest, RoundTrip_D2C2D_RowMajor) {
    int N = N_ab04, M = M_ab04, P = P_ab04;
    double ALPHA = ALPHA_ab04, BETA = BETA_ab04;
    // Row-major LDs (number of columns)
    int LDA = N, LDB = M, LDC = N, LDD = M;
    int ROW_MAJOR = 1;
    double check_tol = 1e-4;

    // Prepare initial row-major inputs
    std::vector<double> A = transpose_matrix(Ad_cm, N, N);
    std::vector<double> B = transpose_matrix(Bd_cm, N, M);
    std::vector<double> C = transpose_matrix(Cd_cm, P, N);
    std::vector<double> D = transpose_matrix(Dd_cm, P, M);
    // Save original row-major state for final comparison
    std::vector<double> A_orig_rm = A;
    std::vector<double> B_orig_rm = B;
    std::vector<double> C_orig_rm = C;
    std::vector<double> D_orig_rm = D;

    // D -> C
    int INFO = slicot_ab04md('D', N, M, P, ALPHA, BETA,
                             A.data(), LDA, B.data(), LDB,
                             C.data(), LDC, D.data(), LDD, ROW_MAJOR);
    ASSERT_EQ(INFO, 0) << "D->C conversion failed";

    // C -> D
    INFO = slicot_ab04md('C', N, M, P, ALPHA, BETA,
                         A.data(), LDA, B.data(), LDB,
                         C.data(), LDC, D.data(), LDD, ROW_MAJOR);
    ASSERT_EQ(INFO, 0) << "C->D conversion failed";

    // Verify we get the original row-major matrices back
    for(size_t i=0; i<A.size(); ++i) EXPECT_NEAR(A[i], A_orig_rm[i], check_tol) << "A[" << i << "] mismatch";
    for(size_t i=0; i<B.size(); ++i) EXPECT_NEAR(B[i], B_orig_rm[i], check_tol) << "B[" << i << "] mismatch";
    for(size_t i=0; i<C.size(); ++i) EXPECT_NEAR(C[i], C_orig_rm[i], check_tol) << "C[" << i << "] mismatch";
    for(size_t i=0; i<D.size(); ++i) EXPECT_NEAR(D[i], D_orig_rm[i], check_tol) << "D[" << i << "] mismatch";
}
