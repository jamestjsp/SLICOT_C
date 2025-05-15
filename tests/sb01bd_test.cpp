#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm> // For std::max
#include <stdexcept> // For std::runtime_error

#include "sb01bd.h"
#include "slicot_utils.h"

// --- Column-Major Test Fixture ---
class SB01BDTestColMajor : public ::testing::Test {
protected:
    char DICO_in = 'C';    // Continuous-time system
    int N_in = 4;          // Order of the state matrix A
    int M_in = 2;          // Number of inputs
    int NP_in = 2;         // Number of poles to assign
    double ALPHA_in = -2.5; // Stability margin for eigenvalues (NFP should be 2 with A_in)
    double TOL_in = 1e-6;  // Tolerance for numerical computations
    int IWARN_out = -1;    // Warning indicator
    int INFO_out = -999;   // Result from slicot_sb01bd call

    // Input matrices (column-major format)
    std::vector<double> A_in;
    std::vector<double> B_in;
    std::vector<double> WR_in;
    std::vector<double> WI_in;

    // Outputs
    std::vector<double> F_out;
    std::vector<double> Z_out;
    int NFP_out = -1;
    int NAP_out = -1;
    int NUP_out = -1;

    // Expected results
    std::vector<double> F_expected;
    std::vector<double> Z_expected;
    int EXPECTED_NFP = 2;  // With A_in=diag(-1,-2,-3,-4) and ALPHA_in=-2.5, eigenvalues -3,-4 are "fixed"
    int EXPECTED_NAP = 2;  // N-NFP = 4-2=2. MIN(NP_in, N-NFP) = MIN(2,2)=2.
    int EXPECTED_NUP = 0;
    int EXPECTED_INFO = 0;
    int EXPECTED_IWARN = 0;

    // Leading dimensions
    int LDA_in = 0;
    int LDB_in = 0;
    int LDF_out = 0;
    int LDZ_out = 0;

    // Verification tolerance
    double check_tol = 1e-6;

    void SetUp() override {
        // Set leading dimensions for Column Major (Fortran-style: LD is number of rows)
        LDA_in = std::max(1, N_in);
        LDB_in = std::max(1, N_in); // B is N x M, Fortran LDB >= N
        LDF_out = std::max(1, M_in); // F is M x N, Fortran LDF >= M
        LDZ_out = std::max(1, N_in); // Z is N x N, Fortran LDZ >= N

        // Initialize input matrices (column-major format)
        A_in = { // diag(-1, -2, -3, -4)
            -1.0, 0.0, 0.0, 0.0,
             0.0, -2.0, 0.0, 0.0,
             0.0,  0.0, -3.0, 0.0,
             0.0,  0.0,  0.0, -4.0
        };
        B_in = { // B = [1 0; 0 1; 0 0; 0 0] (col-major)
            1.0, 0.0, 0.0, 0.0, // col 1
            0.0, 1.0, 0.0, 0.0  // col 2
        };
        // Desired eigenvalues (first NP_in are used)
        WR_in = {-1.0, -2.0, -0.5, -0.7}; // Assign -1, -2
        WI_in = {0.0, 0.0, 0.0, 0.0};

        // Resize output matrices
        F_out.resize((size_t)LDF_out * N_in); // M x N
        Z_out.resize((size_t)LDZ_out * N_in); // N x N

        // Expected results for A=diag(-1,-2,-3,-4), B=[I(2);0], ALPHA=-2.5, NP=2, WR/WI={-1,0;-2,0}
        // A11=diag(-3,-4) (fixed), A22=diag(-1,-2) (to modify)
        // Desired poles for A22 are -1, -2, which are already its poles.
        // So, F should be close to zero. Z should be identity.
        F_expected.assign((size_t)LDF_out * N_in, 0.0); // Expect F to be zero
        
        // Z_expected should be the permutation matrix that reorders eigenvalues of A_in
        // from diag(-1, -2, -3, -4) to diag(-3, -4, -2, -1).
        // Fixed poles (-3, -4) come first.
        // Assigned poles (-1, -2) are placed by the routine, potentially reordered.
        // The actual Z produced by SB01BD for this case results in assigned poles as (-2, -1).
        // Z_out corresponds to [e3, e4, e2, -e1]
        // Z_out matrix form:
        // [ 0  0  0 -1 ]
        // [ 0  0  1  0 ]
        // [ 1  0  0  0 ]
        // [ 0  1  0  0 ]
        // In column-major:
        Z_expected = { 
            0.0,  0.0,  1.0,  0.0, // col 0 (e3)
            0.0,  0.0,  0.0,  1.0, // col 1 (e4)
            0.0,  1.0,  0.0,  0.0, // col 2 (e2)
           -1.0,  0.0,  0.0,  0.0  // col 3 (-e1)
        };
    }
};

// --- Row-Major Test Fixture ---
class SB01BDTestRowMajor : public SB01BDTestColMajor {
protected:
    std::vector<double> A_rm_in;
    std::vector<double> B_rm_in;
    std::vector<double> F_rm_out;
    std::vector<double> Z_rm_out;

    std::vector<double> F_expected_rm;
    std::vector<double> Z_expected_rm;

    void SetUp() override {
        SB01BDTestColMajor::SetUp(); // Base class SetUp called first

        // For row-major C arrays, the 'leading dimension' passed to the wrapper
        // is the number of columns.
        // A is N x N. C LDA (cols) = N.
        LDA_in = std::max(1, N_in);
        // B is N x M. C LDB (cols) = M.
        LDB_in = std::max(1, M_in);
        // F is M x N. C LDF (cols) = N.
        LDF_out = std::max(1, N_in); 
        // Z is N x N. C LDZ (cols) = N.
        LDZ_out = std::max(1, N_in);

        // Resize row-major matrices based on their C dimensions (rows x C_LD_cols)
        A_rm_in.resize((size_t)N_in * LDA_in); 
        B_rm_in.resize((size_t)N_in * LDB_in); 
        F_rm_out.resize((size_t)M_in * LDF_out); 
        Z_rm_out.resize((size_t)N_in * LDZ_out); 

        // Transpose column-major inputs (from base SetUp) to row-major
        // A_in (N_in x N_in), col-major from base. A_rm_in (N_in x LDA_in), row-major.
        if (!A_in.empty()) slicot_transpose_to_c_with_ld(A_in.data(), A_rm_in.data(), N_in, N_in, std::max(1,N_in), LDA_in, sizeof(double));
        // B_in (N_in x M_in), col-major from base. B_rm_in (N_in x LDB_in), row-major.
        if (!B_in.empty()) slicot_transpose_to_c_with_ld(B_in.data(), B_rm_in.data(), N_in, M_in, std::max(1,N_in), LDB_in, sizeof(double));

        // Transpose expected results (from base SetUp) to row-major
        // F_expected (M_in x N_in), col-major. F_expected_rm (M_in x LDF_out), row-major.
        F_expected_rm.resize((size_t)M_in * LDF_out); 
        if (!F_expected.empty()) slicot_transpose_to_c_with_ld(F_expected.data(), F_expected_rm.data(), M_in, N_in, std::max(1,M_in), LDF_out, sizeof(double));
        
        // Z_expected (N_in x N_in), col-major. Z_expected_rm (N_in x LDZ_out), row-major.
        Z_expected_rm.resize((size_t)N_in * LDZ_out); 
        if (!Z_expected.empty()) slicot_transpose_to_c_with_ld(Z_expected.data(), Z_expected_rm.data(), N_in, N_in, std::max(1,N_in), LDZ_out, sizeof(double));
    }
};

// --- Test Cases ---

// Test: Documentation Example (Column-Major)
TEST_F(SB01BDTestColMajor, DocExample) {
    // LDB_in is already correctly set in SetUp for ColMajor
    // LDF_out is already correctly set in SetUp for ColMajor

    INFO_out = slicot_sb01bd(
        DICO_in, N_in, M_in, NP_in, ALPHA_in,
        A_in.data(), LDA_in, B_in.data(), LDB_in, // Use LDB_in from SetUp
        WR_in.data(), WI_in.data(),
        &NFP_out, &NAP_out, &NUP_out,
        F_out.data(), LDF_out, Z_out.data(), LDZ_out, // Use LDF_out from SetUp
        TOL_in, &IWARN_out, 0 /* column-major */
    );

    ASSERT_EQ(INFO_out, EXPECTED_INFO) << "slicot_sb01bd call failed with INFO = " << INFO_out;
    EXPECT_EQ(NFP_out, EXPECTED_NFP);
    EXPECT_EQ(NAP_out, EXPECTED_NAP);
    EXPECT_EQ(NUP_out, EXPECTED_NUP);
    EXPECT_EQ(IWARN_out, EXPECTED_IWARN);

    for (size_t i = 0; i < F_expected.size(); ++i) {
        EXPECT_NEAR(F_out[i], F_expected[i], check_tol) << "F mismatch at index " << i;
    }
    for (size_t i = 0; i < Z_expected.size(); ++i) {
        EXPECT_NEAR(Z_out[i], Z_expected[i], check_tol) << "Z mismatch at index " << i;
    }
}

// Test: Documentation Example (Row-Major)
TEST_F(SB01BDTestRowMajor, DocExample) {
    // ALPHA_in is inherited and set in base SetUp. No need to change it here.
    // LDA_in, LDB_in, LDF_out, LDZ_out are correctly set in RowMajor::SetUp for C row-major calls.

    INFO_out = slicot_sb01bd(
        DICO_in, N_in, M_in, NP_in, ALPHA_in,
        A_rm_in.data(), LDA_in, B_rm_in.data(), LDB_in,
        WR_in.data(), WI_in.data(),
        &NFP_out, &NAP_out, &NUP_out,
        F_rm_out.data(), LDF_out, Z_rm_out.data(), LDZ_out,
        TOL_in, &IWARN_out, 1 /* row-major */
    );

    ASSERT_EQ(INFO_out, EXPECTED_INFO) << "slicot_sb01bd call failed with INFO = " << INFO_out;
    EXPECT_EQ(NFP_out, EXPECTED_NFP);
    EXPECT_EQ(NAP_out, EXPECTED_NAP);
    EXPECT_EQ(NUP_out, EXPECTED_NUP);
    EXPECT_EQ(IWARN_out, EXPECTED_IWARN);

    for (size_t i = 0; i < F_expected_rm.size(); ++i) {
        EXPECT_NEAR(F_rm_out[i], F_expected_rm[i], check_tol) << "F_rm mismatch at index " << i;
    }
    for (size_t i = 0; i < Z_expected_rm.size(); ++i) {
        EXPECT_NEAR(Z_rm_out[i], Z_expected_rm[i], check_tol) << "Z_rm mismatch at index " << i;
    }
}

// Test: Zero Dimensions
TEST_F(SB01BDTestColMajor, ZeroDimensions) {
    int n_zero = 0, m_zero = 0, np_zero = 0;
    int lda_z = 1, ldb_z = 1, ldf_z = 1, ldz_z = 1; // Min leading dims
    int nfp_z = -1, nap_z = -1, nup_z = -1;
    int iwarn_z = -1;
    double alpha_z = 0.0; // Valid alpha for DICO='C' or 'D'
    char dico_z = 'C';

    INFO_out = slicot_sb01bd(
        dico_z, n_zero, m_zero, np_zero, alpha_z,
        nullptr, lda_z, nullptr, ldb_z,
        nullptr, nullptr,
        &nfp_z, &nap_z, &nup_z,
        nullptr, ldf_z, nullptr, ldz_z,
        TOL_in, &iwarn_z, 0 /* column-major */
    );

    EXPECT_EQ(INFO_out, 0); // Expect success for zero dimensions
    EXPECT_EQ(nfp_z, 0);
    EXPECT_EQ(nap_z, 0);
    EXPECT_EQ(nup_z, 0);
    EXPECT_EQ(iwarn_z, 0);
}

// Test: Parameter Validation
TEST_F(SB01BDTestColMajor, ParameterValidation) {
    // Test invalid N
    INFO_out = slicot_sb01bd(
        DICO_in, -1, M_in, NP_in, ALPHA_in,
        nullptr, 1, nullptr, 1, // Dummy A, B, LDA, LDB
        nullptr, nullptr,           // Dummy WR, WI
        &NFP_out, &NAP_out, &NUP_out,
        nullptr, 1, nullptr, 1, // Dummy F, Z, LDF, LDZ
        TOL_in, &IWARN_out, 0 
    );
    EXPECT_EQ(INFO_out, -2);

    // Test invalid M
    INFO_out = slicot_sb01bd(DICO_in, N_in, -1, NP_in, ALPHA_in, A_in.data(), LDA_in, B_in.data(), LDB_in, WR_in.data(), WI_in.data(), &NFP_out, &NAP_out, &NUP_out, F_out.data(), LDF_out, Z_out.data(), LDZ_out, TOL_in, &IWARN_out, 0);
    EXPECT_EQ(INFO_out, -3);
    
    // Test invalid NP (<0)
    INFO_out = slicot_sb01bd(DICO_in, N_in, M_in, -1, ALPHA_in, A_in.data(), LDA_in, B_in.data(), LDB_in, WR_in.data(), WI_in.data(), &NFP_out, &NAP_out, &NUP_out, F_out.data(), LDF_out, Z_out.data(), LDZ_out, TOL_in, &IWARN_out, 0);
    EXPECT_EQ(INFO_out, -4);

    // Test invalid NP (>N)
    INFO_out = slicot_sb01bd(DICO_in, N_in, M_in, N_in + 1, ALPHA_in, A_in.data(), LDA_in, B_in.data(), LDB_in, WR_in.data(), WI_in.data(), &NFP_out, &NAP_out, &NUP_out, F_out.data(), LDF_out, Z_out.data(), LDZ_out, TOL_in, &IWARN_out, 0);
    EXPECT_EQ(INFO_out, -4);
    
    // Test invalid ALPHA for DICO='D'
    INFO_out = slicot_sb01bd('D', N_in, M_in, NP_in, -0.1, A_in.data(), LDA_in, B_in.data(), LDB_in, WR_in.data(), WI_in.data(), &NFP_out, &NAP_out, &NUP_out, F_out.data(), LDF_out, Z_out.data(), LDZ_out, TOL_in, &IWARN_out, 0);
    EXPECT_EQ(INFO_out, -5);

    // Test invalid DICO
    INFO_out = slicot_sb01bd(
        'X', N_in, M_in, NP_in, ALPHA_in,
        A_in.data(), LDA_in, B_in.data(), LDB_in,
        WR_in.data(), WI_in.data(),
        &NFP_out, &NAP_out, &NUP_out,
        F_out.data(), LDF_out, Z_out.data(), LDZ_out,
        TOL_in, &IWARN_out, 0 /* column-major */
    );
    EXPECT_EQ(INFO_out, -1);

    // Test invalid LDA (col-major)
    if (N_in > 0) {
      INFO_out = slicot_sb01bd( DICO_in, N_in, M_in, NP_in, ALPHA_in, A_in.data(), 0, B_in.data(), LDB_in, WR_in.data(), WI_in.data(), &NFP_out, &NAP_out, &NUP_out, F_out.data(), LDF_out, Z_out.data(), LDZ_out, TOL_in, &IWARN_out, 0);
      EXPECT_EQ(INFO_out, -7);
    }

    // Test invalid LDB (col-major)
    if (N_in > 0) { // LDB depends on N for col-major
      INFO_out = slicot_sb01bd( DICO_in, N_in, M_in, NP_in, ALPHA_in, A_in.data(), LDA_in, B_in.data(), 0, WR_in.data(), WI_in.data(), &NFP_out, &NAP_out, &NUP_out, F_out.data(), LDF_out, Z_out.data(), LDZ_out, TOL_in, &IWARN_out, 0);
      EXPECT_EQ(INFO_out, -9);
    }

    // Test invalid LDF (col-major)
    if (M_in > 0) { // LDF depends on M for col-major
      INFO_out = slicot_sb01bd( DICO_in, N_in, M_in, NP_in, ALPHA_in, A_in.data(), LDA_in, B_in.data(), LDB_in, WR_in.data(), WI_in.data(), &NFP_out, &NAP_out, &NUP_out, F_out.data(), 0, Z_out.data(), LDZ_out, TOL_in, &IWARN_out, 0);
      EXPECT_EQ(INFO_out, -16);
    }
    
    // Test invalid LDZ (col-major)
    if (N_in > 0) {
      INFO_out = slicot_sb01bd( DICO_in, N_in, M_in, NP_in, ALPHA_in, A_in.data(), LDA_in, B_in.data(), LDB_in, WR_in.data(), WI_in.data(), &NFP_out, &NAP_out, &NUP_out, F_out.data(), LDF_out, Z_out.data(), 0, TOL_in, &IWARN_out, 0);
      EXPECT_EQ(INFO_out, -18);
    }

    // Test NULL A when N > 0
    if (N_in > 0) {
        INFO_out = slicot_sb01bd(DICO_in, N_in, M_in, NP_in, ALPHA_in, nullptr, LDA_in, B_in.data(), LDB_in, WR_in.data(), WI_in.data(), &NFP_out, &NAP_out, &NUP_out, F_out.data(), LDF_out, Z_out.data(), LDZ_out, TOL_in, &IWARN_out, 0);
        EXPECT_EQ(INFO_out, -6);
    }
    // Test NULL B when N > 0 and M > 0
    if (N_in > 0 && M_in > 0) {
        INFO_out = slicot_sb01bd(DICO_in, N_in, M_in, NP_in, ALPHA_in, A_in.data(), LDA_in, nullptr, LDB_in, WR_in.data(), WI_in.data(), &NFP_out, &NAP_out, &NUP_out, F_out.data(), LDF_out, Z_out.data(), LDZ_out, TOL_in, &IWARN_out, 0);
        EXPECT_EQ(INFO_out, -8);
    }
    // Test NULL WR when NP > 0
    if (NP_in > 0) {
        INFO_out = slicot_sb01bd(DICO_in, N_in, M_in, NP_in, ALPHA_in, A_in.data(), LDA_in, B_in.data(), LDB_in, nullptr, WI_in.data(), &NFP_out, &NAP_out, &NUP_out, F_out.data(), LDF_out, Z_out.data(), LDZ_out, TOL_in, &IWARN_out, 0);
        EXPECT_EQ(INFO_out, -10); // Assuming WR is 10th, WI is 11th
    }
     // Test NULL F when N > 0 and M > 0
    if (N_in > 0 && M_in > 0) {
        INFO_out = slicot_sb01bd(DICO_in, N_in, M_in, NP_in, ALPHA_in, A_in.data(), LDA_in, B_in.data(), LDB_in, WR_in.data(), WI_in.data(), &NFP_out, &NAP_out, &NUP_out, nullptr, LDF_out, Z_out.data(), LDZ_out, TOL_in, &IWARN_out, 0);
        EXPECT_EQ(INFO_out, -15);
    }
    // Test NULL Z when N > 0
    if (N_in > 0) {
        INFO_out = slicot_sb01bd(DICO_in, N_in, M_in, NP_in, ALPHA_in, A_in.data(), LDA_in, B_in.data(), LDB_in, WR_in.data(), WI_in.data(), &NFP_out, &NAP_out, &NUP_out, F_out.data(), LDF_out, nullptr, LDZ_out, TOL_in, &IWARN_out, 0);
        EXPECT_EQ(INFO_out, -17);
    }
    // Test NULL IWARN
    INFO_out = slicot_sb01bd(DICO_in, N_in, M_in, NP_in, ALPHA_in, A_in.data(), LDA_in, B_in.data(), LDB_in, WR_in.data(), WI_in.data(), &NFP_out, &NAP_out, &NUP_out, F_out.data(), LDF_out, Z_out.data(), LDZ_out, TOL_in, nullptr, 0);
    EXPECT_EQ(INFO_out, -21);

    // Test NULL NFP
    INFO_out = slicot_sb01bd(DICO_in, N_in, M_in, NP_in, ALPHA_in, A_in.data(), LDA_in, B_in.data(), LDB_in, WR_in.data(), WI_in.data(), nullptr, &NAP_out, &NUP_out, F_out.data(), LDF_out, Z_out.data(), LDZ_out, TOL_in, &IWARN_out, 0);
    EXPECT_EQ(INFO_out, -12); // NFP is 12th argument
}
