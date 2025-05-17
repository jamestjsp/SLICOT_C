#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm> // For std::max
#include <stdexcept> // For std::runtime_error
#include <iomanip>   // For std::fixed and std::setprecision

#include "sb02mt.h"
#include "slicot_utils.h"

// Helper to print matrix for debugging
void printMatrix(const std::string& name, const std::vector<double>& matrix, int rows, int cols, int ld, bool is_row_major) {
    std::cout << name << " (" << rows << "x" << cols << ", LD=" << ld << (is_row_major ? ", RowMajor" : ", ColMajor") << "):\n";
    if (matrix.empty() && rows > 0 && cols > 0) {
        std::cout << "  <empty or uninitialized>\n";
        return;
    }
     if (rows == 0 || cols == 0) {
        std::cout << "  <zero dimension>\n";
        return;
    }
    std::cout << std::fixed << std::setprecision(7);
    for (int i = 0; i < rows; ++i) {
        std::cout << "  ";
        for (int j = 0; j < cols; ++j) {
            if (is_row_major) {
                 if ( (size_t)(i * ld + j) < matrix.size()) std::cout << matrix[i * ld + j] << "\t"; else std::cout << "OOB\t";
            } else {
                 if ( (size_t)(j * ld + i) < matrix.size()) std::cout << matrix[j * ld + i] << "\t"; else std::cout << "OOB\t";
            }
        }
        std::cout << "\n";
    }
    std::cout << std::defaultfloat << std::setprecision(6);
}


// --- Column-Major Test Fixture ---
class SB02MTTestColMajor : public ::testing::Test {
protected:
    char JOBG_in = 'G';
    char JOBL_in = 'Z';
    char FACT_in = 'N';
    char UPLO_in = 'U';
    int N_in = 3;
    int M_in = 2;

    double check_tol = 1e-6;

    int INFO_out = -999;
    int OUFACT_out = 0;
    int EXPECTED_OUFACT = 1; // From Slycot

    std::vector<double> A_io_orig; // Store original A for checking non-modification when JOBL='Z'
    std::vector<double> Q_io_orig; // Store original Q for checking non-modification when JOBL='Z'
    std::vector<double> B_io_orig; // Store original B for checking B_b
    std::vector<double> R_io_orig; // Store original R for checking R_b


    std::vector<double> A_io;
    std::vector<double> B_io;
    std::vector<double> Q_io;
    std::vector<double> R_io;
    std::vector<double> L_io;
    std::vector<int>    IPIV_io;
    std::vector<double> G_out;

    std::vector<double> A_expected;
    std::vector<double> Q_expected;
    std::vector<double> G_expected;
    std::vector<double> B_expected_modified; // B*chol(R)^-1
    std::vector<double> R_expected_modified; // chol(R)


    int LDA_in = 0;
    int LDB_in = 0;
    int LDQ_in = 0;
    int LDR_in = 0;
    int LDL_in = 0;
    int LDG_out = 0;

    void SetUp() override {
        LDA_in = std::max(1, N_in); LDB_in = std::max(1, N_in); LDQ_in = std::max(1, N_in);
        LDR_in = std::max(1, M_in); LDL_in = std::max(1, N_in); LDG_out = std::max(1, N_in);

        A_io_orig = { // This is column-major for A = [[0.5, 0.2, 0.0], [0.0, 0.5, 0.0], [0.1, 0.1, 0.8]]
            0.5, 0.0, 0.1,
            0.2, 0.5, 0.1,
            0.0, 0.0, 0.8
        };
        A_io = A_io_orig;
        if (N_in > 0) A_io.resize((size_t)LDA_in * N_in); else A_io.assign((size_t)LDA_in * N_in, 0.0);


        B_io_orig = { // Column-major for B = [[1.0, 0.0], [0.0, 1.0], [0.5, 0.0]]
            1.0, 0.0, 0.5,
            0.0, 1.0, 0.0
        };
        B_io = B_io_orig;
        if (N_in > 0 && M_in > 0) B_io.resize((size_t)LDB_in * M_in); else B_io.assign((size_t)LDB_in*M_in, 0.0);

        Q_io_orig = { // Column-major for Q = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0
        };
        Q_io = Q_io_orig;
        if (N_in > 0) Q_io.resize((size_t)LDQ_in * N_in); else Q_io.assign((size_t)LDQ_in*N_in, 0.0);

        R_io_orig = { // Column-major for R = [[2.0, 0.0], [0.0, 2.0]]
            2.0, 0.0,
            0.0, 2.0
        };
        R_io = R_io_orig;
        if (M_in > 0) R_io.resize((size_t)LDR_in * M_in); else R_io.assign((size_t)LDR_in*M_in, 0.0);

        if (N_in > 0 && M_in > 0) L_io.resize((size_t)LDL_in * M_in, 0.0); else L_io.assign((size_t)LDL_in*M_in, 0.0);
        if (JOBL_in == 'Z' && LDL_in == 0 && N_in==0 && M_in==0) L_io.clear();


        IPIV_io.resize(std::max(1, M_in), 0);

        if (JOBG_in == 'G' && N_in > 0) G_out.resize((size_t)LDG_out * N_in);
        else if (JOBG_in == 'G' && N_in == 0) G_out.assign((size_t)LDG_out * N_in, 0.0);
        else G_out.clear();

        // Expected results from Slycot
        // A should be unchanged as JOBL='Z'
        A_expected = A_io_orig;
        if (N_in > 0) A_expected.resize((size_t)LDA_in * N_in); else A_expected.assign((size_t)LDA_in*N_in,0.0);


        // Q should be unchanged as JOBL='Z'
        Q_expected = Q_io_orig;
        if (N_in > 0) Q_expected.resize((size_t)LDQ_in * N_in); else Q_expected.assign((size_t)LDQ_in*N_in,0.0);


        if (JOBG_in == 'G') {
            G_expected = { 0.5, 0.0, 0.25, 0.0, 0.5, 0.0, 0.25, 0.0, 0.125 }; // Column Major
            if (N_in > 0) G_expected.resize((size_t)LDG_out * N_in); else G_expected.assign((size_t)LDG_out*N_in,0.0);
        }

        // Expected B_modified = B*chol(R)^-1 (if oufact=1)
        // B_b (Fortran order from Slycot): [[0.70710678 0. ], [0. 0.70710678], [0.35355339 0. ]]
        B_expected_modified = {0.70710678, 0.0, 0.35355339, 0.0, 0.70710678, 0.0};
        if (N_in > 0 && M_in > 0) B_expected_modified.resize((size_t)LDB_in*M_in); else B_expected_modified.assign((size_t)LDB_in*M_in,0.0);


        // Expected R_modified = chol(R) (if oufact=1)
        // R_b (Fortran order from Slycot): [[1.41421356 0. ], [0. 1.41421356]]
        R_expected_modified = {1.41421356, 0.0, 0.0, 1.41421356};
        if (M_in>0) R_expected_modified.resize((size_t)LDR_in*M_in); else R_expected_modified.assign((size_t)LDR_in*M_in,0.0);

        EXPECTED_OUFACT = 1;
    }
};

// --- Row-Major Test Fixture ---
class SB02MTTestRowMajor : public SB02MTTestColMajor {
protected:
    std::vector<double> A_rm_io; std::vector<double> B_rm_io; std::vector<double> Q_rm_io;
    std::vector<double> R_rm_io; std::vector<double> L_rm_io; std::vector<double> G_rm_out;

    std::vector<double> A_expected_rm; std::vector<double> Q_expected_rm; std::vector<double> G_expected_rm;
    std::vector<double> B_expected_modified_rm; std::vector<double> R_expected_modified_rm;


    void SetUp() override {
        SB02MTTestColMajor::SetUp();

        int lda_c = (N_in > 0) ? N_in : 1; int ldb_c = (M_in > 0) ? M_in : 1; int ldq_c = (N_in > 0) ? N_in : 1;
        int ldr_c = (M_in > 0) ? M_in : 1; int ldl_c = (M_in > 0) ? M_in : 1; int ldg_c = (N_in > 0) ? N_in : 1;

        auto resize_and_transpose = [&](
            const std::vector<double>& col_major_src, std::vector<double>& row_major_dest,
            int rows, int cols, int ld_col, int ld_row) {
            if (rows > 0 && cols > 0) {
                row_major_dest.resize((size_t)rows * ld_row);
                if (!col_major_src.empty()) { // Ensure source is not empty before transposing
                     slicot_transpose_to_c_with_ld(col_major_src.data(), row_major_dest.data(), rows, cols, ld_col, ld_row, sizeof(double));
                } else { // if source is empty (e.g. G_expected if JOBG='N') fill with placeholder
                    std::fill(row_major_dest.begin(), row_major_dest.end(), -999.0); // Or some other indicator
                }
            } else {
                row_major_dest.clear();
            }
        };

        resize_and_transpose(A_io, A_rm_io, N_in, N_in, LDA_in, lda_c);
        resize_and_transpose(B_io, B_rm_io, N_in, M_in, LDB_in, ldb_c);
        resize_and_transpose(Q_io, Q_rm_io, N_in, N_in, LDQ_in, ldq_c);
        resize_and_transpose(R_io, R_rm_io, M_in, M_in, LDR_in, ldr_c);
        resize_and_transpose(L_io, L_rm_io, N_in, M_in, LDL_in, ldl_c);

        if (JOBG_in == 'G') { G_rm_out.resize((size_t)N_in * ldg_c); } else { G_rm_out.clear(); }


        resize_and_transpose(A_expected, A_expected_rm, N_in, N_in, LDA_in, lda_c);
        resize_and_transpose(Q_expected, Q_expected_rm, N_in, N_in, LDQ_in, ldq_c);
        if (JOBG_in == 'G') { resize_and_transpose(G_expected, G_expected_rm, N_in, N_in, LDG_out, ldg_c); }
        else { G_expected_rm.clear(); }

        resize_and_transpose(B_expected_modified, B_expected_modified_rm, N_in, M_in, LDB_in, ldb_c);
        resize_and_transpose(R_expected_modified, R_expected_modified_rm, M_in, M_in, LDR_in, ldr_c);


        LDA_in = lda_c; LDB_in = ldb_c; LDQ_in = ldq_c;
        LDR_in = ldr_c; LDL_in = ldl_c; LDG_out = ldg_c;
    }
};

// --- Test Cases ---
TEST_F(SB02MTTestColMajor, BasicFunctionality) {
    // Make copies of inputs that might be modified if we want to check their original vs modified state
    std::vector<double> A_io_copy = A_io;
    std::vector<double> Q_io_copy = Q_io;
    std::vector<double> B_io_copy = B_io;
    std::vector<double> R_io_copy = R_io;

    INFO_out = slicot_sb02mt(
        JOBG_in, JOBL_in, FACT_in, UPLO_in, N_in, M_in,
        A_io.data(), LDA_in, B_io.data(), LDB_in, Q_io.data(), LDQ_in,
        R_io.data(), LDR_in, L_io.data(), LDL_in, IPIV_io.data(),
        G_out.data(), LDG_out, &OUFACT_out, 0 );

    ASSERT_EQ(INFO_out, 0) << "slicot_sb02mt failed with INFO = " << INFO_out;
    ASSERT_EQ(OUFACT_out, EXPECTED_OUFACT);

    // Verify A_io (should be unchanged if JOBL_in == 'Z')
    if (JOBL_in == 'Z' && N_in > 0) {
        for (size_t i = 0; i < A_io_orig.size(); ++i) { // Compare A_io (potentially modified) with A_expected (original A)
            EXPECT_NEAR(A_io[i], A_expected[i], check_tol) << "A_io mismatch at index " << i << " (expected original A)";
        }
    }
    // Verify Q_io (should be unchanged if JOBL_in == 'Z')
    if (JOBL_in == 'Z' && N_in > 0) {
         for (size_t i = 0; i < Q_io_orig.size(); ++i) {
            EXPECT_NEAR(Q_io[i], Q_expected[i], check_tol) << "Q_io mismatch at index " << i << " (expected original Q)";
        }
    }
    // Verify G_out if computed
    if (JOBG_in == 'G' && N_in > 0 && !G_out.empty() && !G_expected.empty()) {
        for (int j = 0; j < N_in; ++j) {
            for (int i = 0; i < N_in; ++i) {
                if ((UPLO_in == 'U' && i <= j) || (UPLO_in == 'L' && i >= j)) {
                    EXPECT_NEAR(G_out[j * LDG_out + i], G_expected[j * LDG_out + i], check_tol)
                        << "G_out mismatch at (" << i << "," << j << ") (stored part)";
                }
            }
        }
    }
    // Verify modified B_io and R_io if OUFACT indicates modification
    if (OUFACT_out == 1) { // Cholesky factor stored in R_io, B_io updated
        if (M_in > 0 && !R_io.empty() && !R_expected_modified.empty()){
            for (int j=0; j < M_in; ++j) for (int i=0; i < M_in; ++i) {
                 if ((UPLO_in == 'U' && i <= j) || (UPLO_in == 'L' && i >= j)) {
                    EXPECT_NEAR(R_io[j*LDR_in+i], R_expected_modified[j*LDR_in+i], check_tol) << "R_io (chol(R)) mismatch at ("<<i<<","<<j<<")";
                 }
            }
        }
        if (N_in > 0 && M_in > 0 && !B_io.empty() && !B_expected_modified.empty()){
             for (size_t i=0; i < B_expected_modified.size(); ++i) {
                EXPECT_NEAR(B_io[i], B_expected_modified[i], check_tol) << "B_io (B*inv(chol(R))) mismatch at index "<<i;
             }
        }
    }
}

TEST_F(SB02MTTestRowMajor, BasicFunctionality) {
    std::vector<double> A_rm_io_copy = A_rm_io; // For checking original values if needed
    std::vector<double> Q_rm_io_copy = Q_rm_io;

    INFO_out = slicot_sb02mt(
        JOBG_in, JOBL_in, FACT_in, UPLO_in, N_in, M_in,
        A_rm_io.data(), LDA_in, B_rm_io.data(), LDB_in, Q_rm_io.data(), LDQ_in,
        R_rm_io.data(), LDR_in, L_rm_io.data(), LDL_in, IPIV_io.data(),
        G_rm_out.data(), LDG_out, &OUFACT_out, 1 );

    ASSERT_EQ(INFO_out, 0) << "slicot_sb02mt failed with INFO = " << INFO_out;
    ASSERT_EQ(OUFACT_out, EXPECTED_OUFACT);

    if (JOBL_in == 'Z' && N_in > 0 && !A_rm_io.empty() && !A_expected_rm.empty()) {
        for (size_t i = 0; i < A_expected_rm.size(); ++i) {
            EXPECT_NEAR(A_rm_io[i], A_expected_rm[i], check_tol) << "A_rm_io mismatch at index " << i;
        }
    }
    if (JOBL_in == 'Z' && N_in > 0 && !Q_rm_io.empty() && !Q_expected_rm.empty()) {
        for (size_t i = 0; i < Q_expected_rm.size(); ++i) {
            EXPECT_NEAR(Q_rm_io[i], Q_expected_rm[i], check_tol) << "Q_rm_io mismatch at index " << i;
        }
    }

    if (JOBG_in == 'G' && N_in > 0 && !G_rm_out.empty() && !G_expected_rm.empty()) {
        for (int i = 0; i < N_in; ++i) { // Row index
            for (int j = 0; j < N_in; ++j) { // Col index
                if ((UPLO_in == 'U' && i <= j) || (UPLO_in == 'L' && i >= j)) {
                     EXPECT_NEAR(G_rm_out[i * LDG_out + j], G_expected_rm[i * LDG_out + j], check_tol)
                        << "G_rm_out mismatch at (" << i << "," << j << ") (stored part)";
                }
            }
        }
    }
     if (OUFACT_out == 1) {
        if (M_in > 0 && !R_rm_io.empty() && !R_expected_modified_rm.empty()){
            for (int i=0; i < M_in; ++i) for (int j=0; j < M_in; ++j) {
                 if ((UPLO_in == 'U' && i <= j) || (UPLO_in == 'L' && i >= j)) {
                    EXPECT_NEAR(R_rm_io[i*LDR_in+j], R_expected_modified_rm[i*LDR_in+j], check_tol) << "R_rm_io (chol(R)) mismatch at ("<<i<<","<<j<<")";
                 }
            }
        }
        if (N_in > 0 && M_in > 0 && !B_rm_io.empty() && !B_expected_modified_rm.empty()){
             for (size_t i=0; i < B_expected_modified_rm.size(); ++i) {
                EXPECT_NEAR(B_rm_io[i], B_expected_modified_rm[i], check_tol) << "B_rm_io (B*inv(chol(R))) mismatch at index "<<i;
             }
        }
    }
}

TEST_F(SB02MTTestColMajor, ZeroDimensions) {
    int n_zero = 0, m_zero = 0;
    int lda_z = 1, ldb_z = 1, ldq_z = 1, ldr_z = 1, ldl_z = 1, ldg_z = 1;
    double dummy_val = 0.0;
    std::vector<int> ipiv_z(1,0);
    int oufact_z = 0;
    int expected_oufact_z = 0; // If M=0, oufact=0

    INFO_out = slicot_sb02mt(
        'N', 'Z', 'N', 'U',
        n_zero, m_zero,
        &dummy_val, lda_z, &dummy_val, ldb_z, &dummy_val, ldq_z,
        &dummy_val, ldr_z, &dummy_val, ldl_z, ipiv_z.data(),
        &dummy_val, ldg_z, &oufact_z, 0 );
    EXPECT_EQ(INFO_out, 0);
    EXPECT_EQ(oufact_z, expected_oufact_z);
}

TEST_F(SB02MTTestColMajor, ParameterValidation) {
    std::vector<int> ipiv_val(std::max(1,M_in),0);
    int oufact_val = 0;

    INFO_out = slicot_sb02mt('X', JOBL_in, FACT_in, UPLO_in, N_in, M_in,
                           A_io.data(), LDA_in, B_io.data(), LDB_in, Q_io.data(), LDQ_in,
                           R_io.data(), LDR_in, L_io.data(), LDL_in, ipiv_val.data(),
                           G_out.data(), LDG_out, &oufact_val, 0);
    EXPECT_EQ(INFO_out, -1); // Invalid JOBG

    INFO_out = slicot_sb02mt(JOBG_in, JOBL_in, FACT_in, UPLO_in, -1, M_in, // Invalid N
                           A_io.data(), LDA_in, B_io.data(), LDB_in, Q_io.data(), LDQ_in,
                           R_io.data(), LDR_in, L_io.data(), LDL_in, ipiv_val.data(),
                           G_out.data(), LDG_out, &oufact_val, 0);
    EXPECT_EQ(INFO_out, -5); // N mapped to arg 5

    char orig_jobl = JOBL_in; JOBL_in = 'N'; // A_io is used
    if (N_in > 0) {
        INFO_out = slicot_sb02mt(JOBG_in, JOBL_in, FACT_in, UPLO_in, N_in, M_in,
                               A_io.data(), 0, B_io.data(), LDB_in, Q_io.data(), LDQ_in, // Invalid LDA
                               R_io.data(), LDR_in, L_io.data(), LDL_in, ipiv_val.data(),
                               G_out.data(), LDG_out, &oufact_val, 0);
        EXPECT_EQ(INFO_out, -8); // LDA mapped to arg 8
    }
    JOBL_in = orig_jobl;

    JOBL_in = 'N';
    if (N_in > 0) {
        INFO_out = slicot_sb02mt(JOBG_in, JOBL_in, FACT_in, UPLO_in, N_in, M_in,
                               nullptr, LDA_in, B_io.data(), LDB_in, Q_io.data(), LDQ_in, // NULL A_io
                               R_io.data(), LDR_in, L_io.data(), LDL_in, ipiv_val.data(),
                               G_out.data(), LDG_out, &oufact_val, 0);
        EXPECT_EQ(INFO_out, -7); // A_io mapped to arg 7
    }
    JOBL_in = orig_jobl;
}