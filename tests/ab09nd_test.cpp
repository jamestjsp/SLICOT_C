#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <limits> // Required for std::numeric_limits
#include <algorithm> // Required for std::min, std::max
#include <string> // Required for std::to_string
#include <iostream> // For printing matrix values during debugging (optional)

// Include the specific header for the function being tested
#include "ab09nd.h" // Header for slicot_ab09nd

// Helper function to print matrix (optional, for debugging)
void printMatrixD(const std::string& name, const double* data, int rows, int cols, int ld, bool rowMajor) {
    std::cout << name << " (" << rows << "x" << cols << ", ld=" << ld << ", " << (rowMajor ? "RowMajor" : "ColMajor") << "):\n";
    for (int i = 0; i < rows; ++i) {
        std::cout << "  [";
        for (int j = 0; j < cols; ++j) {
            double val = rowMajor ? data[i * ld + j] : data[i + j * ld];
            // Format output for better readability
            std::cout << std::fixed << std::setprecision(4) << std::setw(9) << val << (j == cols - 1 ? "" : ", ");
        }
        std::cout << "]\n";
    }
     std::cout << std::defaultfloat << std::setprecision(6); // Reset default formatting
}


// =============================================================================
// Test Fixture and Test Case for Column-Major Input (ROW_MAJOR = 0)
// =============================================================================

class AB09NDTestColMajor : public ::testing::Test {
protected:
    // Define common variables for tests if needed
    int N = 7, M = 2, P = 3;
    char DICO = 'C', JOB = 'N', EQUIL = 'N', ORDSEL = 'A';
    double ALPHA = -0.6, TOL1 = 0.1, TOL2 = 1.0e-14;
    int ROW_MAJOR = 0;
    double check_tol = 1e-4; // Tolerance for checking results

    // Expected results from documentation
    int expected_NR = 5;
    int expected_NS = 5;
    int expected_IWARN = 0;
    std::vector<double> expected_HSV = {1.9178, 0.8621, 0.7666, 0.0336, 0.0246};
    std::vector<double> expected_Ar = {
        -0.5181,  8.8157,  0.0000,  0.0000,  0.0000, // Col 1
        -1.1084, -0.5181,  0.0000,  0.0000,  0.0000, // Col 2
         0.0000,  0.0000,  0.5847,  0.0000, -4.3823, // Col 3
         0.0000,  0.0000,  0.0000, -1.6606,  0.0000, // Col 4
         0.0000,  0.0000,  1.9230,  0.0000, -3.2922  // Col 5
    }; // 5x5, Column Major
    std::vector<double> expected_Br = {
        -1.2837, -0.7522, -0.6379,  2.0656, -3.9315, // Col 1
         1.2837,  0.7522, -0.6379, -2.0656, -3.9315  // Col 2
    }; // 5x2, Column Major
    std::vector<double> expected_Cr = {
        -0.1380,  0.6246,  0.1380, // Col 1
        -0.6445,  0.0196,  0.6445, // Col 2
        -0.6416,  0.0000, -0.6416, // Col 3
        -0.6293,  0.4107,  0.6293, // Col 4
         0.2526,  0.0000,  0.2526  // Col 5
    }; // 3x5, Column Major
    std::vector<double> expected_Dr = {
         0.0582,  0.0015, -0.0090, // Col 1
        -0.0090, -0.0015,  0.0582  // Col 2
    }; // 3x2, Column Major
};

TEST_F(AB09NDTestColMajor, DocExampleColMajor) {

    // --- Leading Dimensions for Column-Major (Fortran style) ---
    int LDA = N; // >= N
    int LDB = N; // >= N
    int LDC = P; // >= P
    int LDD = P; // >= P

    // --- Input Matrices (Column-major order, N=7, M=2, P=3) ---
    // A (7x7, LDA=7)
    std::vector<double> A = {
        -0.04165, -5.2100,  0.0000,  0.5450,  0.0000,  0.0000,  0.0000, // Col 1
         0.0000, -12.500,  3.3300,  0.0000,  0.0000,  0.0000,  0.0000, // Col 2
         4.9200,  0.0000, -3.3300,  0.0000,  0.0000,  0.0000,  0.0000, // Col 3
        -4.9200,  0.0000,  0.0000,  0.0000,  4.9200,  0.0000,  0.0000, // Col 4
         0.0000,  0.0000,  0.0000, -0.5450, -0.04165, -5.2100,  0.0000, // Col 5
         0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -12.500,  3.3300, // Col 6
         0.0000,  0.0000,  0.0000,  0.0000,  4.9200,  0.0000, -3.3300  // Col 7
    };
    ASSERT_EQ(A.size(), LDA * N);

    // B (7x2, LDB=7)
    std::vector<double> B = {
         0.0000, 12.500,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, // Col 1
         0.0000,  0.0000,  0.0000,  0.0000,  0.0000, 12.500,  0.0000  // Col 2
    };
    ASSERT_EQ(B.size(), LDB * M);

    // C (3x7, LDC=3)
    std::vector<double> C = {
        1.0000, 0.0000, 0.0000, // Col 1
        0.0000, 0.0000, 0.0000, // Col 2
        0.0000, 0.0000, 0.0000, // Col 3
        0.0000, 1.0000, 0.0000, // Col 4
        0.0000, 0.0000, 1.0000, // Col 5
        0.0000, 0.0000, 0.0000, // Col 6
        0.0000, 0.0000, 0.0000  // Col 7
    };
    ASSERT_EQ(C.size(), LDC * N);

    // D (3x2, LDD=3)
    std::vector<double> D(LDD * M, 0.0); // Zero matrix
    ASSERT_EQ(D.size(), LDD * M);

    // Input/Output NR for ORDSEL='A' (Input value doesn't matter much)
    int nr_io = 0;

    // Output variables
    int NS_out;
    int IWARN_out;
    std::vector<double> hsv_out(N);

    // --- Call the SLICOT C wrapper ---
    int INFO = slicot_ab09nd(DICO, JOB, EQUIL, ORDSEL, N, M, P, &nr_io, ALPHA,
                            A.data(), LDA, B.data(), LDB, C.data(), LDC, D.data(), LDD,
                            &NS_out, hsv_out.data(), TOL1, TOL2,
                            &IWARN_out, ROW_MAJOR);

    // --- Check results ---
    ASSERT_EQ(INFO, 0) << "SLICOT routine ab09nd returned error code: " << INFO;

    // Check scalar outputs
    EXPECT_EQ(nr_io, expected_NR) << "Reduced order NR mismatch";
    EXPECT_EQ(NS_out, expected_NS) << "Stable dimension NS mismatch";
    EXPECT_EQ(IWARN_out, expected_IWARN) << "Warning flag IWARN mismatch";

    // Check Hankel Singular Values (first NS_out values)
    ASSERT_GE(hsv_out.size(), NS_out);
    ASSERT_GE(expected_HSV.size(), NS_out);
    for (int i = 0; i < NS_out; ++i) {
        EXPECT_NEAR(hsv_out[i], expected_HSV[i], check_tol) << "HSV[" << i << "] mismatch";
    }

    // Check reduced matrices (overwritten in A, B, C, D)
    // Ar (NR x NR)
    ASSERT_GE(A.size(), LDA * nr_io); // Check size before access
    ASSERT_GE(expected_Ar.size(), expected_NR * expected_NR);
    for (int j = 0; j < expected_NR; ++j) { // Col
        for (int i = 0; i < expected_NR; ++i) { // Row
            EXPECT_NEAR(A[i + j * LDA], expected_Ar[i + j * expected_NR], check_tol)
                << "Ar[" << i << "," << j << "] mismatch (ColMajor)";
        }
    }

    // Br (NR x M)
    ASSERT_GE(B.size(), LDB * M);
    ASSERT_GE(expected_Br.size(), expected_NR * M);
     for (int j = 0; j < M; ++j) { // Col
        for (int i = 0; i < expected_NR; ++i) { // Row
            EXPECT_NEAR(B[i + j * LDB], expected_Br[i + j * expected_NR], check_tol)
                << "Br[" << i << "," << j << "] mismatch (ColMajor)";
        }
    }

    // Cr (P x NR)
    ASSERT_GE(C.size(), LDC * nr_io);
    ASSERT_GE(expected_Cr.size(), P * expected_NR);
     for (int j = 0; j < expected_NR; ++j) { // Col
        for (int i = 0; i < P; ++i) { // Row
            EXPECT_NEAR(C[i + j * LDC], expected_Cr[i + j * P], check_tol)
                << "Cr[" << i << "," << j << "] mismatch (ColMajor)";
        }
    }

    // Dr (P x M)
    ASSERT_GE(D.size(), LDD * M);
    ASSERT_GE(expected_Dr.size(), P * M);
     for (int j = 0; j < M; ++j) { // Col
        for (int i = 0; i < P; ++i) { // Row
            EXPECT_NEAR(D[i + j * LDD], expected_Dr[i + j * P], check_tol)
                << "Dr[" << i << "," << j << "] mismatch (ColMajor)";
        }
    }
}

// =============================================================================
// Test Fixture and Test Case for Row-Major Input (ROW_MAJOR = 1)
// =============================================================================

// Inherit expected values and parameters from ColMajor fixture
class AB09NDTestRowMajor : public AB09NDTestColMajor {
public:
    AB09NDTestRowMajor() {
        ROW_MAJOR = 1; // Override ROW_MAJOR for this fixture
    }
};


TEST_F(AB09NDTestRowMajor, DocExampleRowMajor) {

    // --- Leading Dimensions for Row-Major (C style) ---
    int LDA = N; // Cols for A (NxN)
    int LDB = M; // Cols for B (NxM)
    int LDC = N; // Cols for C (PxN)
    int LDD = M; // Cols for D (PxM)

    // --- Input Matrices (Row-major order, N=7, M=2, P=3) ---
    // A (7x7, LDA=7)
     std::vector<double> A = {
        -0.04165,  0.0000,  4.9200, -4.9200,  0.0000,  0.0000,  0.0000, // Row 0
        -5.2100, -12.500,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, // Row 1
         0.0000,  3.3300, -3.3300,  0.0000,  0.0000,  0.0000,  0.0000, // Row 2
         0.5450,  0.0000,  0.0000,  0.0000, -0.5450,  0.0000,  0.0000, // Row 3
         0.0000,  0.0000,  0.0000,  4.9200, -0.04165,  0.0000,  4.9200, // Row 4
         0.0000,  0.0000,  0.0000,  0.0000, -5.2100, -12.500,  0.0000, // Row 5
         0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  3.3300, -3.3300  // Row 6
    };
    ASSERT_EQ(A.size(), N * LDA);

    // B (7x2, LDB=2)
    std::vector<double> B = {
         0.0000,  0.0000, // Row 0
        12.500,  0.0000, // Row 1
         0.0000,  0.0000, // Row 2
         0.0000,  0.0000, // Row 3
         0.0000,  0.0000, // Row 4
         0.0000, 12.500, // Row 5
         0.0000,  0.0000  // Row 6
    };
    ASSERT_EQ(B.size(), N * LDB);

    // C (3x7, LDC=7)
    std::vector<double> C = {
        1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, // Row 0
        0.0000, 0.0000, 0.0000, 1.0000, 0.0000, 0.0000, 0.0000, // Row 1
        0.0000, 0.0000, 0.0000, 0.0000, 1.0000, 0.0000, 0.0000  // Row 2
    };
    ASSERT_EQ(C.size(), P * LDC);

    // D (3x2, LDD=2)
    std::vector<double> D(P * LDD, 0.0); // Zero matrix
    ASSERT_EQ(D.size(), P * LDD);

    // Input/Output NR for ORDSEL='A'
    int nr_io = 0;

    // Output variables
    int NS_out;
    int IWARN_out;
    std::vector<double> hsv_out(N);

    // --- Call the SLICOT C wrapper ---
    int INFO = slicot_ab09nd(DICO, JOB, EQUIL, ORDSEL, N, M, P, &nr_io, ALPHA,
                            A.data(), LDA, B.data(), LDB, C.data(), LDC, D.data(), LDD,
                            &NS_out, hsv_out.data(), TOL1, TOL2,
                            &IWARN_out, ROW_MAJOR);

    // --- Check results ---
    ASSERT_EQ(INFO, 0) << "SLICOT routine ab09nd returned error code: " << INFO;

    // Check scalar outputs
    EXPECT_EQ(nr_io, expected_NR) << "Reduced order NR mismatch";
    EXPECT_EQ(NS_out, expected_NS) << "Stable dimension NS mismatch";
    EXPECT_EQ(IWARN_out, expected_IWARN) << "Warning flag IWARN mismatch";

    // Check Hankel Singular Values (first NS_out values)
    ASSERT_GE(hsv_out.size(), NS_out);
    ASSERT_GE(expected_HSV.size(), NS_out);
    for (int i = 0; i < NS_out; ++i) {
        EXPECT_NEAR(hsv_out[i], expected_HSV[i], check_tol) << "HSV[" << i << "] mismatch";
    }

    // --- Convert expected results to Row Major for comparison ---
    std::vector<double> expected_Ar_rm(expected_NR * expected_NR);
    std::vector<double> expected_Br_rm(expected_NR * M);
    std::vector<double> expected_Cr_rm(P * expected_NR);
    std::vector<double> expected_Dr_rm(P * M);

    slicot_transpose_to_c(expected_Ar.data(), expected_Ar_rm.data(), expected_NR, expected_NR, sizeof(double));
    slicot_transpose_to_c(expected_Br.data(), expected_Br_rm.data(), expected_NR, M, sizeof(double));
    slicot_transpose_to_c(expected_Cr.data(), expected_Cr_rm.data(), P, expected_NR, sizeof(double));
    slicot_transpose_to_c(expected_Dr.data(), expected_Dr_rm.data(), P, M, sizeof(double));


    // Check reduced matrices (overwritten in A, B, C, D)
    // Ar (NR x NR)
    ASSERT_GE(A.size(), nr_io * LDA); // Check size before access (N rows, LDA cols)
    ASSERT_GE(expected_Ar_rm.size(), expected_NR * expected_NR);
    // printMatrixD("Actual Ar (RowMajor)", A.data(), expected_NR, expected_NR, LDA, true); // Optional debug print
    for (int i = 0; i < expected_NR; ++i) { // Row
        for (int j = 0; j < expected_NR; ++j) { // Col
            EXPECT_NEAR(A[i * LDA + j], expected_Ar_rm[i * expected_NR + j], check_tol)
                << "Ar[" << i << "," << j << "] mismatch (RowMajor)";
        }
    }

    // Br (NR x M)
    ASSERT_GE(B.size(), nr_io * LDB);
    ASSERT_GE(expected_Br_rm.size(), expected_NR * M);
    // printMatrixD("Actual Br (RowMajor)", B.data(), expected_NR, M, LDB, true); // Optional debug print
     for (int i = 0; i < expected_NR; ++i) { // Row
        for (int j = 0; j < M; ++j) { // Col
            EXPECT_NEAR(B[i * LDB + j], expected_Br_rm[i * M + j], check_tol)
                << "Br[" << i << "," << j << "] mismatch (RowMajor)";
        }
    }

    // Cr (P x NR)
    ASSERT_GE(C.size(), P * LDC); // C is P rows, LDC cols
    ASSERT_GE(expected_Cr_rm.size(), P * expected_NR);
    // printMatrixD("Actual Cr (RowMajor)", C.data(), P, expected_NR, LDC, true); // Optional debug print
     for (int i = 0; i < P; ++i) { // Row
        for (int j = 0; j < expected_NR; ++j) { // Col
             EXPECT_NEAR(C[i * LDC + j], expected_Cr_rm[i * expected_NR + j], check_tol)
                << "Cr[" << i << "," << j << "] mismatch (RowMajor)";
        }
    }

    // Dr (P x M)
    ASSERT_GE(D.size(), P * LDD);
    ASSERT_GE(expected_Dr_rm.size(), P * M);
    // printMatrixD("Actual Dr (RowMajor)", D.data(), P, M, LDD, true); // Optional debug print
     for (int i = 0; i < P; ++i) { // Row
        for (int j = 0; j < M; ++j) { // Col
            EXPECT_NEAR(D[i * LDD + j], expected_Dr_rm[i * M + j], check_tol)
                << "Dr[" << i << "," << j << "] mismatch (RowMajor)";
        }
    }
}
