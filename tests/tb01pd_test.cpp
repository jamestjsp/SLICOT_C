#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm> // For std::max, std::min
#include <iostream>
#include <iomanip>   // For std::fixed, std::setprecision

#include "tb01pd.h"
#include "slicot_utils.h" // For transpose functions

// Helper to print matrices for debugging (optional, can be commented out)
void print_matrix_tb01pd_test(const std::string& name, const std::vector<double>& mat, int r, int c, int ld_c, bool is_row_major) {
    std::cout << name << " (" << r << "x" << c << ") " << (is_row_major ? "RowMajor" : "ColMajor") << ", LD=" << ld_c << ":\n";
    if (mat.empty() || r == 0 || c == 0) {
        std::cout << "  (empty or zero dim)\n";
        return;
    }
    std::cout << std::fixed << std::setprecision(4);
    for (int i = 0; i < r; ++i) {
        std::cout << "  [";
        for (int j = 0; j < c; ++j) {
            double val = is_row_major ? mat[i * ld_c + j] : mat[i + j * ld_c];
            std::cout << val << (j == c - 1 ? "" : ", ");
        }
        std::cout << "]\n";
    }
    std::cout << std::resetiosflags(std::ios_base::fixed) << std::setprecision(6); // Reset to default
}

// --- Test Fixture Base ---
class TB01PDTest : public ::testing::Test {
protected:
    // Parameters for the SLICOT routine
    char JOB_param = 'M';
    char EQUIL_param = 'N';
    int N_param = 0;
    int M_param = 0;
    int P_param = 0;
    double TOL_param = 0.0;

    // Input data (stored in row-major for easy initialization in tests)
    std::vector<double> A_data_rm, B_data_rm, C_data_rm;

    // Buffers for C function call (these will be copies/transpositions of _data_rm)
    std::vector<double> A_test_buf, B_test_buf, C_test_buf;
    int NR_out_val; // Output: order of the reduced system

    // Expected results (stored in row-major for easy definition)
    std::vector<double> Ar_expected_rm, Br_expected_rm, Cr_expected_rm;
    int NR_expected;
    int info_expected = 0; // Expected return code from the wrapper

    // Leading dimensions for C arrays passed to the wrapper
    int LDA, LDB, LDC;

    // Tolerance for floating point comparisons
    double check_tol_matrix = 1e-4; // Based on example output precision

    // Method to set up data based on the SLICOT TB01PD.html example
    void setup_slicot_example_data() {
        JOB_param = 'M';
        EQUIL_param = 'N';
        N_param = 3;
        M_param = 1;
        P_param = 2;
        TOL_param = 0.0; // Use default tolerance in the routine

        // Original matrices from SLICOT documentation example
        A_data_rm = {1.0, 2.0, 0.0,
                     4.0, -1.0, 0.0,
                     0.0, 0.0, 1.0}; // N x N = 3 x 3

        B_data_rm = {1.0,
                     0.0,
                     1.0}; // N x M = 3 x 1

        C_data_rm = {0.0, 1.0, -1.0,
                     0.0, 0.0, 1.0}; // P x N = 2 x 3

        // Expected results from SLICOT documentation example
        NR_expected = 3;
        Ar_expected_rm = {1.0000, -1.4142,  1.4142,
                         -2.8284,  0.0000,  1.0000,
                          2.8284,  1.0000,  0.0000}; // NR x NR = 3 x 3

        Br_expected_rm = {-1.0000,
                           0.7071,
                           0.7071}; // NR x M = 3 x 1

        Cr_expected_rm = {0.0000,  0.0000, -1.4142,
                          0.0000,  0.7071,  0.7071}; // P x NR = 2 x 3
        info_expected = 0;
    }

    // Helper to prepare C call buffers for row-major tests
    void prepare_rowmajor_buffers_for_slicot_example() {
        // For row-major C arrays, LDA is number of columns of A
        LDA = N_param;

        // LDB for C row-major is number of columns of B_test_buf.
        // The Fortran routine uses B for workspace if JOB='M' or 'O',
        // requiring its column dimension to be MAX(M,P).
        // The C wrapper expects the input B array (b_io) to have N_param rows
        // and LDB columns, where LDB must accommodate this workspace requirement.
        LDB = (JOB_param == 'C') ? M_param : std::max(M_param, P_param);
        LDB = std::max(1, LDB); // Ensure LDB is at least 1

        // LDC for C row-major is number of columns of C_test_buf (which is N_param).
        LDC = N_param;
        LDC = std::max(1,LDC);


        // Copy A_data_rm to A_test_buf (no transpose needed for A's content initially)
        A_test_buf = A_data_rm;

        // Size B_test_buf for N_param rows and LDB columns (C definition for row-major)
        B_test_buf.assign((size_t)N_param * LDB, 0.0); // Initialize with zeros
        if (M_param > 0 && N_param > 0) {
            // Copy original B_data_rm (N_param x M_param) into B_test_buf
            for (int i = 0; i < N_param; ++i) {
                for (int j = 0; j < M_param; ++j) {
                    // B_data_rm has M_param columns. B_test_buf has LDB columns.
                    if ((size_t)(i * M_param + j) < B_data_rm.size()) {
                         if ((size_t)(i * LDB + j) < B_test_buf.size()) {
                            B_test_buf[i * LDB + j] = B_data_rm[i * M_param + j];
                         }
                    }
                }
            }
        }

        // C_test_buf is P_param x N_param. LDC is N_param.
        C_test_buf = C_data_rm;
    }

    // Helper to prepare C call buffers for column-major tests
    void prepare_colmajor_buffers_for_slicot_example() {
        // For column-major C arrays (Fortran-style), LDA is number of rows of A
        LDA = std::max(1, N_param);
        LDB = std::max(1, N_param); // Fortran LDB must be >= N

        // Fortran LDC for C matrix:
        // From TB01PD.html: LDC >= MAX(1,M,P) if N > 0. LDC >= 1 if N = 0.
        // This is the leading dimension of the C array in Fortran.
        LDC = (N_param == 0) ? 1 : std::max(1, std::max(M_param, P_param));


        // Transpose A_data_rm (N_param x N_param, row-major) to A_test_buf (col-major)
        A_test_buf.resize((size_t)LDA * N_param);
        if (N_param > 0) {
            slicot_transpose_to_fortran_with_ld(A_data_rm.data(), A_test_buf.data(), N_param, N_param, N_param, LDA, sizeof(double));
        } else {
            A_test_buf.clear();
        }

        // Determine Fortran's required column dimension for B
        int b_fortran_cols = (JOB_param == 'C') ? M_param : std::max(M_param, P_param);
        b_fortran_cols = std::max(0, b_fortran_cols); // Ensure non-negative

        // Size B_test_buf for LDB rows and b_fortran_cols columns (Fortran definition)
        if (N_param > 0 && b_fortran_cols > 0) {
            B_test_buf.resize((size_t)LDB * b_fortran_cols);
            std::fill(B_test_buf.begin(), B_test_buf.end(), 0.0); // Initialize with zeros
            if (M_param > 0) {
                 // Transpose B_data_rm (N_param x M_param, row-major) to a temporary col-major B
                 std::vector<double> b_temp_cm((size_t)N_param * M_param);
                 slicot_transpose_to_fortran_with_ld(B_data_rm.data(), b_temp_cm.data(), N_param, M_param, M_param, N_param, sizeof(double));
                 // Copy the M_param columns from b_temp_cm into B_test_buf
                 for (int j_f = 0; j_f < M_param; ++j_f) {
                    for (int i_f = 0; i_f < N_param; ++i_f) {
                        if (((size_t)i_f + (size_t)j_f * LDB) < B_test_buf.size() && ((size_t)i_f + (size_t)j_f * N_param) < b_temp_cm.size())
                           B_test_buf[i_f + j_f * LDB] = b_temp_cm[i_f + j_f * N_param];
                    }
                 }
            }
        } else {
            B_test_buf.clear();
        }

        // Transpose C_data_rm (P_param x N_param, row-major) to C_test_buf (col-major)
        // C_test_buf will be LDC rows x N_param columns (Fortran definition)
        if (P_param > 0 && N_param > 0) {
            C_test_buf.resize((size_t)LDC * N_param);
            std::fill(C_test_buf.begin(), C_test_buf.end(), 0.0); // Initialize with zeros
            slicot_transpose_to_fortran_with_ld(C_data_rm.data(), C_test_buf.data(), P_param, N_param, N_param, LDC, sizeof(double));
        } else if (N_param == 0 && P_param > 0) { // Case for C when N=0, P>0
             C_test_buf.resize((size_t)LDC * std::max(1,N_param)); // LDC x N (N could be 0, so max(1,N))
             std::fill(C_test_buf.begin(), C_test_buf.end(), 0.0);
        }
        else {
            C_test_buf.clear();
        }
    }
};

TEST_F(TB01PDTest, SlicotDocExampleRowMajor) {
    setup_slicot_example_data();
    prepare_rowmajor_buffers_for_slicot_example();

    // std::cout << "Row Major Test - Before call:\n";
    // print_matrix_tb01pd_test("A_test_buf (in)", A_test_buf, N_param, N_param, LDA, true);
    // print_matrix_tb01pd_test("B_test_buf (in)", B_test_buf, N_param, LDB, LDB, true); // LDB is cols of B_test_buf
    // print_matrix_tb01pd_test("C_test_buf (in)", C_test_buf, P_param, N_param, LDC, true);


    int info = slicot_tb01pd(JOB_param, EQUIL_param, N_param, M_param, P_param,
                             A_test_buf.data(), LDA,
                             B_test_buf.empty() ? nullptr : B_test_buf.data(), LDB,
                             C_test_buf.empty() ? nullptr : C_test_buf.data(), LDC,
                             &NR_out_val, TOL_param,
                             1 /* row_major = true */);

    // std::cout << "Row Major Test - After call (info=" << info << ", NR_out_val=" << NR_out_val << "):\n";
    // print_matrix_tb01pd_test("A_test_buf (Ar_rm)", A_test_buf, NR_out_val, NR_out_val, LDA, true);
    // print_matrix_tb01pd_test("B_test_buf (Br_rm)", B_test_buf, NR_out_val, M_param, LDB, true);
    // print_matrix_tb01pd_test("C_test_buf (Cr_rm)", C_test_buf, P_param, NR_out_val, LDC, true);


    ASSERT_EQ(info, info_expected);
    ASSERT_EQ(NR_out_val, NR_expected);

    // Verify Ar (output A_test_buf is NR_expected x NR_expected, row-major, LDA is original N_param cols)
    // The wrapper transposes back into the leading NR_expected x NR_expected part.
    // The C array A_test_buf has LDA columns.
    for (int i_row = 0; i_row < NR_expected; ++i_row) {
        for (int j_col = 0; j_col < NR_expected; ++j_col) {
            double actual_val = A_test_buf[i_row * LDA + j_col];
            double expected_val = Ar_expected_rm[i_row * NR_expected + j_col];
            EXPECT_NEAR(actual_val, expected_val, check_tol_matrix)
                << "Mismatch at Ar_rm(row " << i_row << ", col " << j_col << ")";
        }
    }

    // Verify Br (output B_test_buf is NR_expected x M_param, row-major, LDB is original LDB cols)
    // The C array B_test_buf has LDB columns.
    if (M_param > 0 && NR_expected > 0) {
        for (int i_row = 0; i_row < NR_expected; ++i_row) {
            for (int j_col = 0; j_col < M_param; ++j_col) {
                double actual_val = B_test_buf[i_row * LDB + j_col];
                double expected_val = Br_expected_rm[i_row * M_param + j_col];
                EXPECT_NEAR(actual_val, expected_val, check_tol_matrix)
                    << "Mismatch at Br_rm(row " << i_row << ", col " << j_col << ")";
            }
        }
    }

    // Verify Cr (output C_test_buf is P_param x NR_expected, row-major, LDC is original N_param cols)
    // The C array C_test_buf has LDC columns.
    if (P_param > 0 && NR_expected > 0) {
        for (int i_row = 0; i_row < P_param; ++i_row) {
            for (int j_col = 0; j_col < NR_expected; ++j_col) {
                double actual_val = C_test_buf[i_row * LDC + j_col]; // LDC is N_param
                double expected_val = Cr_expected_rm[i_row * NR_expected + j_col];
                EXPECT_NEAR(actual_val, expected_val, check_tol_matrix)
                    << "Mismatch at Cr_rm(row " << i_row << ", col " << j_col << ")\n"
                    << "  Actual C_test_buf[" << i_row << " * LDC(" << LDC << ") + " << j_col << "] = " << actual_val << "\n"
                    << "  Expected Cr_expected_rm[" << i_row << " * NR_expected(" << NR_expected << ") + " << j_col << "] = " << expected_val;
            }
        }
    }
}

TEST_F(TB01PDTest, SlicotDocExampleColMajor) {
    setup_slicot_example_data();
    prepare_colmajor_buffers_for_slicot_example();

    // std::cout << "Col Major Test - Before call:\n";
    // print_matrix_tb01pd_test("A_test_buf (in_cm)", A_test_buf, N_param, N_param, LDA, false);
    // print_matrix_tb01pd_test("B_test_buf (in_cm)", B_test_buf, N_param, (JOB_param == 'C') ? M_param : std::max(M_param, P_param), LDB, false);
    // print_matrix_tb01pd_test("C_test_buf (in_cm)", C_test_buf, P_param, N_param, LDC, false);

    int info = slicot_tb01pd(JOB_param, EQUIL_param, N_param, M_param, P_param,
                             A_test_buf.data(), LDA,
                             B_test_buf.empty() ? nullptr : B_test_buf.data(), LDB,
                             C_test_buf.empty() ? nullptr : C_test_buf.data(), LDC,
                             &NR_out_val, TOL_param,
                             0 /* row_major = false */);

    // std::cout << "Col Major Test - After call (info=" << info << ", NR_out_val=" << NR_out_val << "):\n";
    // print_matrix_tb01pd_test("A_test_buf (Ar_cm)", A_test_buf, NR_out_val, NR_out_val, LDA, false);
    // print_matrix_tb01pd_test("B_test_buf (Br_cm)", B_test_buf, NR_out_val, M_param, LDB, false);
    // print_matrix_tb01pd_test("C_test_buf (Cr_cm)", C_test_buf, P_param, NR_out_val, LDC, false);


    ASSERT_EQ(info, info_expected);
    ASSERT_EQ(NR_out_val, NR_expected);

    // Expected Ar (transpose Ar_expected_rm to column-major for comparison)
    // Ar_expected_rm is NR_expected x NR_expected, row-major.
    // Ar_exp_cm will be NR_expected x NR_expected, col-major, LD_cm = NR_expected.
    std::vector<double> Ar_exp_cm((size_t)NR_expected * NR_expected);
    if (NR_expected > 0) {
        slicot_transpose_to_fortran_with_ld(Ar_expected_rm.data(), Ar_exp_cm.data(), NR_expected, NR_expected, NR_expected, NR_expected, sizeof(double));
    }
    // A_test_buf contains Ar, is col-major, NR_expected x NR_expected content, actual LD in memory is LDA.
    for (int j_col = 0; j_col < NR_expected; ++j_col) {
        for (int i_row = 0; i_row < NR_expected; ++i_row) {
            EXPECT_NEAR(A_test_buf[i_row + j_col * LDA], Ar_exp_cm[i_row + j_col * NR_expected], check_tol_matrix)
                << "Mismatch at Ar_cm(row " << i_row << ", col " << j_col << ")";
        }
    }

    // Expected Br
    // Br_expected_rm is NR_expected x M_param, row-major.
    // Br_exp_cm will be NR_expected x M_param, col-major, LD_cm = NR_expected.
    if (M_param > 0 && NR_expected > 0) {
        std::vector<double> Br_exp_cm((size_t)NR_expected * M_param);
        slicot_transpose_to_fortran_with_ld(Br_expected_rm.data(), Br_exp_cm.data(), NR_expected, M_param, M_param, NR_expected, sizeof(double));
        // B_test_buf contains Br, is col-major, NR_expected x M_param content, actual LD in memory is LDB.
        for (int j_col = 0; j_col < M_param; ++j_col) {
            for (int i_row = 0; i_row < NR_expected; ++i_row) {
                EXPECT_NEAR(B_test_buf[i_row + j_col * LDB], Br_exp_cm[i_row + j_col * NR_expected], check_tol_matrix)
                    << "Mismatch at Br_cm(row " << i_row << ", col " << j_col << ")";
            }
        }
    }

    // Expected Cr
    // Cr_expected_rm is P_param x NR_expected, row-major.
    // Cr_exp_cm will be P_param x NR_expected, col-major, LD_cm = P_param.
    if (P_param > 0 && NR_expected > 0) {
        std::vector<double> Cr_exp_cm((size_t)P_param * NR_expected);
        slicot_transpose_to_fortran_with_ld(Cr_expected_rm.data(), Cr_exp_cm.data(), P_param, NR_expected, NR_expected, P_param, sizeof(double));

        // C_test_buf contains Cr, is col-major, P_param x NR_expected content, actual LD in memory is LDC.
        for (int j_col = 0; j_col < NR_expected; ++j_col) {
            for (int i_row = 0; i_row < P_param; ++i_row) {
                double actual_val = C_test_buf[i_row + j_col * LDC];
                double expected_val = Cr_exp_cm[i_row + j_col * P_param]; // Use P_param as LD for Cr_exp_cm
                EXPECT_NEAR(actual_val, expected_val, check_tol_matrix)
                    << "Mismatch at Cr_cm(row " << i_row << ", col " << j_col << ")\n"
                    << "  Actual C_test_buf[" << i_row << " + " << j_col << " * LDC(" << LDC << ")] = " << actual_val << "\n"
                    << "  Expected Cr_exp_cm[" << i_row << " + " << j_col << " * P_param(" << P_param << ")] = " << expected_val;
            }
        }
    }
}


TEST_F(TB01PDTest, ParameterValidation) {
    setup_slicot_example_data(); // Initialize some default valid parameters
    // Override N, M, P for specific validation checks if needed, or use example's
    int temp_n = 1, temp_m = 1, temp_p = 1;
    std::vector<double> dummy_A(1,0.0), dummy_B(1,0.0), dummy_C(1,0.0);
    int dummy_lda=1, dummy_ldb=1, dummy_ldc=1;
    int temp_nr_out;

    // Test invalid JOB
    int info = slicot_tb01pd('X', EQUIL_param, temp_n, temp_m, temp_p, dummy_A.data(), dummy_lda, dummy_B.data(), dummy_ldb, dummy_C.data(), dummy_ldc, &temp_nr_out, TOL_param, 1);
    EXPECT_EQ(info, -1); // JOB is 1st arg

    // Test invalid EQUIL
    info = slicot_tb01pd(JOB_param, 'Y', temp_n, temp_m, temp_p, dummy_A.data(), dummy_lda, dummy_B.data(), dummy_ldb, dummy_C.data(), dummy_ldc, &temp_nr_out, TOL_param, 1);
    EXPECT_EQ(info, -2); // EQUIL is 2nd arg

    // Test invalid N
    info = slicot_tb01pd(JOB_param, EQUIL_param, -1, temp_m, temp_p, nullptr, dummy_lda, nullptr, dummy_ldb, nullptr, dummy_ldc, &temp_nr_out, TOL_param, 1);
    EXPECT_EQ(info, -3); // N is 3rd arg

    // Test invalid M
    info = slicot_tb01pd(JOB_param, EQUIL_param, temp_n, -1, temp_p, dummy_A.data(), dummy_lda, nullptr, dummy_ldb, dummy_C.data(), dummy_ldc, &temp_nr_out, TOL_param, 1);
    EXPECT_EQ(info, -4); // M is 4th arg

    // Test invalid P
    info = slicot_tb01pd(JOB_param, EQUIL_param, temp_n, temp_m, -1, dummy_A.data(), dummy_lda, dummy_B.data(), dummy_ldb, nullptr, dummy_ldc, &temp_nr_out, TOL_param, 1);
    EXPECT_EQ(info, -5); // P is 5th arg

    // Test NULL A with N > 0
    if (temp_n > 0) {
        info = slicot_tb01pd(JOB_param, EQUIL_param, temp_n, temp_m, temp_p, nullptr, dummy_lda, dummy_B.data(), dummy_ldb, dummy_C.data(), dummy_ldc, &temp_nr_out, TOL_param, 1);
        EXPECT_EQ(info, -6); // A is 6th arg
    }

    // Test invalid LDA (row-major: LDA < N_param (cols))
    if (temp_n > 0) {
        info = slicot_tb01pd(JOB_param, EQUIL_param, temp_n, temp_m, temp_p, dummy_A.data(), 0, dummy_B.data(), dummy_ldb, dummy_C.data(), dummy_ldc, &temp_nr_out, TOL_param, 1);
        EXPECT_EQ(info, -7); // LDA is 7th arg
    }
     // Test invalid LDA (col-major: LDA < N_param (rows))
    if (temp_n > 0) {
        info = slicot_tb01pd(JOB_param, EQUIL_param, temp_n, temp_m, temp_p, dummy_A.data(), 0, dummy_B.data(), dummy_ldb, dummy_C.data(), dummy_ldc, &temp_nr_out, TOL_param, 0);
        EXPECT_EQ(info, -7); // LDA is 7th arg
    }

    // Test NULL NR_out
    info = slicot_tb01pd(JOB_param, EQUIL_param, temp_n, temp_m, temp_p, dummy_A.data(), dummy_lda, dummy_B.data(), dummy_ldb, dummy_C.data(), dummy_ldc, nullptr, TOL_param, 1);
    EXPECT_EQ(info, -12); // NR is 12th arg
}

TEST_F(TB01PDTest, ZeroDimensionN) {
    setup_slicot_example_data(); // Initialize JOB_param, EQUIL_param, TOL_param
    N_param = 0; M_param = 1; P_param = 1; // Override N
    NR_expected = 0; // Expected NR for N=0

    // For N=0, A, B, C can be NULL. Leading dimensions can be 1.
    LDA=1; LDB=1; LDC=1;

    int info_rm = slicot_tb01pd(JOB_param, EQUIL_param, N_param, M_param, P_param,
                             nullptr, LDA,
                             nullptr, LDB,
                             nullptr, LDC, // C can be NULL if N=0, even if P > 0
                             &NR_out_val, TOL_param,
                             1 /* row_major = true */);
    ASSERT_EQ(info_rm, 0);
    EXPECT_EQ(NR_out_val, NR_expected);

    // Test with column-major as well
    // For N=0, Fortran LDC is MAX(1,M,P) if N>0, or 1 if N=0.
    // So if N=0, LDC for Fortran is 1.
    // The C wrapper will calculate Fortran LDC for C as MAX(1, MAX(M,P)) if N_param=0.
    // This is fine, the C_test_buf can be NULL.
    int col_major_ldc = (N_param == 0) ? 1 : std::max(1, std::max(M_param, P_param));

    int info_cm = slicot_tb01pd(JOB_param, EQUIL_param, N_param, M_param, P_param,
                             nullptr, 1, /* LDA for N=0 */
                             nullptr, 1, /* LDB for N=0 */
                             nullptr, col_major_ldc, /* LDC for N=0, P=1, M=1 -> LDC=1 */
                             &NR_out_val, TOL_param,
                             0 /* row_major = false */);
    ASSERT_EQ(info_cm, 0);
    EXPECT_EQ(NR_out_val, NR_expected);
}

TEST_F(TB01PDTest, ZeroDimensionNMP) {
    setup_slicot_example_data();
    N_param = 0; M_param = 0; P_param = 0;
    NR_expected = 0;
    LDA=1; LDB=1; LDC=1;

    int info_rm = slicot_tb01pd(JOB_param, EQUIL_param, N_param, M_param, P_param,
                             nullptr, LDA, nullptr, LDB, nullptr, LDC,
                             &NR_out_val, TOL_param, 1 /* row_major = true */);
    ASSERT_EQ(info_rm, 0);
    EXPECT_EQ(NR_out_val, NR_expected);
    
    int info_cm = slicot_tb01pd(JOB_param, EQUIL_param, N_param, M_param, P_param,
                             nullptr, LDA, nullptr, LDB, nullptr, LDC,
                             &NR_out_val, TOL_param, 0 /* row_major = false */);
    ASSERT_EQ(info_cm, 0);
    EXPECT_EQ(NR_out_val, NR_expected);
}

