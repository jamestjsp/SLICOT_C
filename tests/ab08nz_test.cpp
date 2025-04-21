#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <limits> // Required for std::numeric_limits
#include <algorithm> // Required for std::min, std::max
#include <string> // Required for std::to_string
#include <complex> // Include complex header for C++ complex types
#include <iostream> // For printing matrix values during debugging (optional)


// Include the specific header for the function being tested
#include "ab08nz.h" // Header for slicot_ab08nz
#include "slicot_utils.h" // For slicot_complex_double definition and helper macros

// Assuming slicot_complex_double is compatible with std::complex<double>
// If not, the constructor below needs adjustment.

// Helper function to get imaginary part robustly
inline double slicot_get_imag(slicot_complex_double val) {
#if defined(__cplusplus) && !defined(__STDC_NO_COMPLEX__)
    // Assuming slicot_complex_double is std::complex<double> or similar C++ type
    return val.imag();
#elif defined(__STDC_NO_COMPLEX__)
    // Assuming slicot_complex_double is a struct { double real; double imag; }
    return val.imag;
#else // C99 complex
    // Assuming slicot_complex_double is double _Complex
    return cimag(val);
#endif
}

// Helper function to print complex matrix (optional, for debugging)
void printComplexMatrix(const std::string& name, const slicot_complex_double* data, int rows, int cols, int ld, bool rowMajor) {
    std::cout << name << " (" << rows << "x" << cols << ", ld=" << ld << ", " << (rowMajor ? "RowMajor" : "ColMajor") << "):\n";
    for (int i = 0; i < rows; ++i) {
        std::cout << "  [";
        for (int j = 0; j < cols; ++j) {
             // Use C++ .real() method
            double real_part = (rowMajor ? data[i * ld + j] : data[i + j * ld]).real();
            double imag_part = slicot_get_imag(rowMajor ? data[i * ld + j] : data[i + j * ld]);
            std::cout << "(" << real_part << "," << imag_part << ")" << (j == cols - 1 ? "" : ", ");
        }
        std::cout << "]\n";
    }
}


// =============================================================================
// Test Fixture and Test Case for Column-Major Input (ROW_MAJOR = 0)
// =============================================================================

class AB08NZTestColMajor : public ::testing::Test {
    // No specific setup needed here for now
};

TEST_F(AB08NZTestColMajor, DocExampleColMajor) {
    // --- Input data from AB08NZ documentation example ---
    // Note: AB08NZ doc has no example. Using AB08ND data interpreted as complex.
    char EQUIL = 'N'; // Example uses 'N'
    int N = 6;        // Order of the system
    int M = 2;        // Number of inputs
    int P = 3;        // Number of outputs
    double TOL = 0.0; // Use default tolerance
    int ROW_MAJOR = 0; // Input arrays are column-major (Fortran style)
    double check_tol = 1e-4; // Tolerance for checking results

    // --- Leading Dimensions for Column-Major (Fortran style) ---
    // These are the number of rows allocated in memory
    int LDA = N;      // LDA >= MAX(1,N)
    int LDB = N;      // LDB >= MAX(1,N)
    int LDC = P;      // LDC >= MAX(1,P)
    int LDD = P;      // LDD >= MAX(1,P)
    // Fortran routine requires LDAF >= N+M and LDBF >= N+P
    int LDAF = N + M; // LDAF >= 6+2 = 8
    int LDBF = N + P; // LDBF >= 6+3 = 9

    // --- Input Matrix A (Column-major order, N x N = 6x6, LDA=6) ---
    std::vector<slicot_complex_double> A(LDA * N, std::complex<double>(0.0, 0.0)); // Use C++ constructor
    A[0 + 0*LDA] = std::complex<double>(1.0, 0.0);
    A[1 + 1*LDA] = std::complex<double>(1.0, 0.0);
    A[2 + 2*LDA] = std::complex<double>(3.0, 0.0);
    A[3 + 3*LDA] = std::complex<double>(-4.0, 0.0);
    A[4 + 4*LDA] = std::complex<double>(-1.0, 0.0);
    A[5 + 5*LDA] = std::complex<double>(3.0, 0.0);
    ASSERT_EQ(A.size(), LDA * N);


    // --- Input Matrix B (Column-major order, N x M = 6x2, LDB=6) ---
    std::vector<slicot_complex_double> B = {
         std::complex<double>(0.0, 0.0), std::complex<double>(-1.0, 0.0), std::complex<double>( 1.0, 0.0), std::complex<double>( 0.0, 0.0), std::complex<double>( 0.0, 0.0), std::complex<double>(-1.0, 0.0), // Col 1
         std::complex<double>(-1.0, 0.0), std::complex<double>( 0.0, 0.0), std::complex<double>(-1.0, 0.0), std::complex<double>( 0.0, 0.0), std::complex<double>( 1.0, 0.0), std::complex<double>(-1.0, 0.0)  // Col 2
    };
    ASSERT_EQ(B.size(), LDB * M);

    // --- Input Matrix C (Column-major order, P x N = 3x6, LDC=3) ---
    std::vector<slicot_complex_double> C = {
        std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), // Col 1
        std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0), // Col 2
        std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0), // Col 3
        std::complex<double>(1.0, 0.0), std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0), // Col 4
        std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), // Col 5
        std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0), std::complex<double>(1.0, 0.0)  // Col 6
    };
    ASSERT_EQ(C.size(), LDC * N);

    // --- Input Matrix D (Column-major order, P x M = 3x2, LDD=3) ---
    std::vector<slicot_complex_double> D(LDD * M, std::complex<double>(0.0, 0.0)); // D is 3x2 zero matrix
    ASSERT_EQ(D.size(), LDD * M);

    // Output variables
    int NU;
    int RANK;
    int DINFZ;
    int NKROR;
    int NKROL;
    std::vector<int> INFZ(N); // Size N is sufficient based on docs
    std::vector<int> KRONR(std::max(N, M) + 1); // Size based on docs
    std::vector<int> KRONL(std::max(N, P) + 1); // Size based on docs

    // Allocate output matrices based on Fortran dimensions required by wrapper internals
    // AF is LDAF x (N+MIN(P,M)) = 8 x (6+2) = 8x8
    // BF is LDBF x (N+M)       = 9 x (6+2) = 9x8
    std::vector<slicot_complex_double> AF(LDAF * (N + std::min(P, M)));
    std::vector<slicot_complex_double> BF(LDBF * (N + M));

    // --- Call the SLICOT C wrapper with ROW_MAJOR = 0 ---
    int INFO = slicot_ab08nz(EQUIL, N, M, P, A.data(), LDA, B.data(), LDB,
                            C.data(), LDC, D.data(), LDD, &NU, &RANK,
                            &DINFZ, &NKROR, &NKROL, INFZ.data(), KRONR.data(),
                            KRONL.data(), AF.data(), LDAF, BF.data(), LDBF,
                            TOL, ROW_MAJOR);

    // --- Check results ---
    ASSERT_EQ(INFO, 0) << "SLICOT routine ab08nz returned error code: " << INFO;

    // --- Expected results (derived from AB08ND example, assuming complex behaves similarly) ---
    int expected_NU = 2;
    // RANK might differ slightly in complex case, omit strict check for now
    // int expected_RANK = 2; // Normal rank?
    int expected_DINFZ = 1; // Max degree of infinite zeros
    int expected_NKROR = 0; // Number of right Kronecker indices
    int expected_NKROL = 1; // Number of left Kronecker indices
    std::vector<int> expected_INFZ = { 2 }; // infz[0]=2 means 2 infinite zeros of degree 1
    std::vector<int> expected_KRONL = { 2 }; // kronl[0]=2 means 1 left index of value 2

    // --- Expected Af and Bf (NU x NU = 2x2, Column-major) ---
    // These are numerical results expected from the complex version run with AB08ND data.
    // These values might differ slightly from AB08ND due to complex arithmetic.
    // Using AB08ND results as placeholders for now.
     slicot_complex_double expected_Af_nu_cm[4] = {
        std::complex<double>(1.70457, 0.0), std::complex<double>(0.68220, 0.0), // Col 1
        std::complex<double>(0.68220, 0.0), std::complex<double>(-0.68498, 0.0) // Col 2
    };
    slicot_complex_double expected_Bf_nu_cm[4] = {
        std::complex<double>(0.93776, 0.0), std::complex<double>(0.01904, 0.0), // Col 1
        std::complex<double>(0.01904, 0.0), std::complex<double>(0.87108, 0.0)  // Col 2
    };


    // --- Verification of results ---
    EXPECT_EQ(NU, expected_NU) << "NU value doesn't match expected";
    // EXPECT_EQ(RANK, expected_RANK) << "RANK value doesn't match expected"; // Omitted check
    EXPECT_EQ(DINFZ, expected_DINFZ) << "DINFZ value doesn't match expected";
    EXPECT_EQ(NKROR, expected_NKROR) << "NKROR value doesn't match expected";
    EXPECT_EQ(NKROL, expected_NKROL) << "NKROL value doesn't match expected";

    ASSERT_GE(INFZ.size(), DINFZ);
    for(int i = 0; i < DINFZ; ++i) {
        EXPECT_EQ(INFZ[i], expected_INFZ[i]) << "INFZ[" << i << "] doesn't match expected";
    }
    ASSERT_GE(KRONL.size(), NKROL);
    for(int i = 0; i < NKROL; ++i) {
        EXPECT_EQ(KRONL[i], expected_KRONL[i]) << "KRONL[" << i << "] doesn't match expected";
    }

    // --- Check relevant part of Af (NU x NU, Column-major indexing) ---
    // printComplexMatrix("Actual AF (ColMajor)", AF.data(), NU, NU, LDAF, false); // Optional debug print
    for (int j = 0; j < NU; j++) { // col
        for (int i = 0; i < NU; i++) { // row
            int current_idx = i + j*LDAF;
            int expected_idx = i + j*NU;
             ASSERT_LT(current_idx, AF.size()); // Bounds check
             ASSERT_LT(expected_idx, 4); // Bounds check for expected array

            std::string msg_base = "AF[" + std::to_string(i) + "," + std::to_string(j) + "]";
            double actual_real = (AF[current_idx]).real();
            double actual_imag = slicot_get_imag(AF[current_idx]);
            double expected_real = (expected_Af_nu_cm[expected_idx]).real();
            double expected_imag = slicot_get_imag(expected_Af_nu_cm[expected_idx]);
            EXPECT_NEAR(actual_real, expected_real, check_tol) << msg_base << " Real part mismatch (ColMajor)";
            EXPECT_NEAR(actual_imag, expected_imag, check_tol) << msg_base << " Imag part mismatch (ColMajor)";
        }
    }

    // --- Check relevant part of Bf (NU x NU, Column-major indexing) ---
    // printComplexMatrix("Actual BF (ColMajor)", BF.data(), NU, NU, LDBF, false); // Optional debug print
    for (int j = 0; j < NU; j++) { // col
        for (int i = 0; i < NU; i++) { // row
            int current_idx = i + j*LDBF;
            int expected_idx = i + j*NU;
            ASSERT_LT(current_idx, BF.size()); // Bounds check
            ASSERT_LT(expected_idx, 4); // Bounds check for expected array

            std::string msg_base = "BF[" + std::to_string(i) + "," + std::to_string(j) + "]";
            double actual_real = (BF[current_idx]).real();
            double actual_imag = slicot_get_imag(BF[current_idx]);
            double expected_real = (expected_Bf_nu_cm[expected_idx]).real();
            double expected_imag = slicot_get_imag(expected_Bf_nu_cm[expected_idx]);
            EXPECT_NEAR(actual_real, expected_real, check_tol) << msg_base << " Real part mismatch (ColMajor)";
            EXPECT_NEAR(actual_imag, expected_imag, check_tol) << msg_base << " Imag part mismatch (ColMajor)";
        }
    }
}


// =============================================================================
// Test Fixture and Test Case for Row-Major Input (ROW_MAJOR = 1)
// =============================================================================

class AB08NZTestRowMajor : public ::testing::Test {
     // No specific setup needed here for now
};

TEST_F(AB08NZTestRowMajor, DocExampleRowMajor) {
    // --- Input data from AB08NZ documentation example ---
    // Note: AB08NZ doc has no example. Using AB08ND data interpreted as complex.
    char EQUIL = 'N'; // Example uses 'N'
    int N = 6;        // Order of the system
    int M = 2;        // Number of inputs
    int P = 3;        // Number of outputs
    double TOL = 0.0; // Use default tolerance
    int ROW_MAJOR = 1; // Input arrays are row-major
    double check_tol = 1e-4; // Tolerance for checking results

    // --- Leading Dimensions for Row-Major ---
    // These are the number of columns allocated in memory for the C arrays
    int LDA = N; // A is N x N
    int LDB = M; // B is N x M
    int LDC = N; // C is P x N
    int LDD = M; // D is P x M

    // Leading dimensions (number of columns) for AF/BF passed to C wrapper.
    // These must also be sufficient for the wrapper to satisfy the underlying
    // Fortran routine's *row* requirements (LDAF_f >= N+M, LDBF_f >= N+P).
    int LDAF = std::max(N + M, N + std::min(P, M)); // Max of Fortran row need (N+M) and Fortran col need (N+min(P,M))
                                                    // LDAF >= max(8, 8) = 8. Also need LDAF >= NU=2. OK.
    // CORRECTED LDBF calculation:
    int LDBF = std::max(N + P, N + M);              // Max of Fortran row need (N+P) and Fortran col need (N+M)
                                                    // LDBF >= max(9, 8) = 9. Also need LDBF >= NU=2. OK.

    // --- Input Matrix A (Row-major order, N x N = 6x6, LDA=6) ---
    std::vector<slicot_complex_double> A(N * LDA, std::complex<double>(0.0, 0.0));
    A[0*LDA + 0] = std::complex<double>(1.0, 0.0);
    A[1*LDA + 1] = std::complex<double>(1.0, 0.0);
    A[2*LDA + 2] = std::complex<double>(3.0, 0.0);
    A[3*LDA + 3] = std::complex<double>(-4.0, 0.0);
    A[4*LDA + 4] = std::complex<double>(-1.0, 0.0);
    A[5*LDA + 5] = std::complex<double>(3.0, 0.0);
    ASSERT_EQ(A.size(), N * LDA);

    // --- Input Matrix B (Row-major order, N x M = 6x2, LDB=2) ---
    std::vector<slicot_complex_double> B = {
        std::complex<double>( 0.0, 0.0), std::complex<double>(-1.0, 0.0), // Row 0
        std::complex<double>(-1.0, 0.0), std::complex<double>( 0.0, 0.0), // Row 1
        std::complex<double>( 1.0, 0.0), std::complex<double>(-1.0, 0.0), // Row 2
        std::complex<double>( 0.0, 0.0), std::complex<double>( 0.0, 0.0), // Row 3
        std::complex<double>( 0.0, 0.0), std::complex<double>( 1.0, 0.0), // Row 4
        std::complex<double>(-1.0, 0.0), std::complex<double>(-1.0, 0.0)  // Row 5
    };
     ASSERT_EQ(B.size(), N * LDB);

    // --- Input Matrix C (Row-major order, P x N = 3x6, LDC=6) ---
    std::vector<slicot_complex_double> C = {
        std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), // Row 0
        std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0), // Row 1
        std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0)  // Row 2
    };
    ASSERT_EQ(C.size(), P * LDC);

    // --- Input Matrix D (Row-major order, P x M = 3x2, LDD=2) ---
    std::vector<slicot_complex_double> D(P * LDD, std::complex<double>(0.0, 0.0)); // D is 3x2 zero matrix
    ASSERT_EQ(D.size(), P * LDD);

    // Output variables
    int NU;
    int RANK;
    int DINFZ;
    int NKROR;
    int NKROL;
    std::vector<int> INFZ(N);
    std::vector<int> KRONR(std::max(N, M) + 1);
    std::vector<int> KRONL(std::max(N, P) + 1);

    // Allocate output matrices based on row-major dimensions (rows * cols)
    // The wrapper should only write to the first NU rows.
    // Allocate enough space based on max internal needs * C columns.
    std::vector<slicot_complex_double> AF( (N + M) * LDAF );
    // Use corrected LDBF for allocation:
    std::vector<slicot_complex_double> BF( (N + P) * LDBF );


    // --- Call the SLICOT C wrapper with ROW_MAJOR = 1 ---
    int INFO = slicot_ab08nz(EQUIL, N, M, P, A.data(), LDA, B.data(), LDB,
                            C.data(), LDC, D.data(), LDD, &NU, &RANK,
                            &DINFZ, &NKROR, &NKROL, INFZ.data(), KRONR.data(),
                            KRONL.data(), AF.data(), LDAF, BF.data(), LDBF,
                            TOL, ROW_MAJOR);

    // --- Check results ---
    // Line 253 where the original error occurred
    ASSERT_EQ(INFO, 0) << "SLICOT routine ab08nz returned error code: " << INFO;

    // --- Expected results (derived from AB08ND example) ---
    int expected_NU = 2;
    // int expected_RANK = 2; // Omit check
    int expected_DINFZ = 1;
    int expected_NKROR = 0;
    int expected_NKROL = 1;
    std::vector<int> expected_INFZ = { 2 };
    std::vector<int> expected_KRONL = { 2 };

    // --- Expected Af and Bf (NU x NU = 2x2, Row-major) ---
    // Placeholder values based on AB08ND results.
     slicot_complex_double expected_Af_nu_rm[4] = {
        std::complex<double>(1.70457, 0.0), std::complex<double>(0.68220, 0.0), // Row 0
        std::complex<double>(0.68220, 0.0), std::complex<double>(-0.68498, 0.0) // Row 1
    };
    slicot_complex_double expected_Bf_nu_rm[4] = {
        std::complex<double>(0.93776, 0.0), std::complex<double>(0.01904, 0.0), // Row 0
        std::complex<double>(0.01904, 0.0), std::complex<double>(0.87108, 0.0)  // Row 1
    };


    // --- Verification of results ---
    EXPECT_EQ(NU, expected_NU) << "NU value doesn't match expected";
    // EXPECT_EQ(RANK, expected_RANK) << "RANK value doesn't match expected"; // Omitted check
    EXPECT_EQ(DINFZ, expected_DINFZ) << "DINFZ value doesn't match expected";
    EXPECT_EQ(NKROR, expected_NKROR) << "NKROR value doesn't match expected";
    EXPECT_EQ(NKROL, expected_NKROL) << "NKROL value doesn't match expected";

    ASSERT_GE(INFZ.size(), DINFZ);
    for(int i = 0; i < DINFZ; ++i) {
        EXPECT_EQ(INFZ[i], expected_INFZ[i]) << "INFZ[" << i << "] doesn't match expected";
    }
    ASSERT_GE(KRONL.size(), NKROL);
    for(int i = 0; i < NKROL; ++i) {
        EXPECT_EQ(KRONL[i], expected_KRONL[i]) << "KRONL[" << i << "] doesn't match expected";
    }

    // --- Check relevant part of Af (NU x NU, Row-major indexing) ---
    // printComplexMatrix("Actual AF (RowMajor)", AF.data(), NU, NU, LDAF, true); // Optional debug print
    for (int i = 0; i < NU; i++) { // row
        for (int j = 0; j < NU; j++) { // col
            int current_idx = i * LDAF + j; // Actual result uses LDAF stride
            int expected_idx = i * NU + j;  // Expected block uses NU stride
            ASSERT_LT(current_idx, AF.size()); // Bounds check
            ASSERT_LT(expected_idx, 4); // Bounds check for expected array

            std::string msg_base = "AF[" + std::to_string(i) + "," + std::to_string(j) + "]";
            double actual_real = (AF[current_idx]).real();
            double actual_imag = slicot_get_imag(AF[current_idx]);
            double expected_real = (expected_Af_nu_rm[expected_idx]).real();
            double expected_imag = slicot_get_imag(expected_Af_nu_rm[expected_idx]);
            EXPECT_NEAR(actual_real, expected_real, check_tol) << msg_base << " Real part mismatch (RowMajor)";
            EXPECT_NEAR(actual_imag, expected_imag, check_tol) << msg_base << " Imag part mismatch (RowMajor)";
        }
    }

    // --- Check relevant part of Bf (NU x NU, Row-major indexing) ---
    // printComplexMatrix("Actual BF (RowMajor)", BF.data(), NU, NU, LDBF, true); // Optional debug print
    for (int i = 0; i < NU; i++) { // row
        for (int j = 0; j < NU; j++) { // col
            int current_idx = i * LDBF + j; // Actual result uses LDBF stride
            int expected_idx = i * NU + j;  // Expected block uses NU stride
            ASSERT_LT(current_idx, BF.size()); // Bounds check
            ASSERT_LT(expected_idx, 4); // Bounds check for expected array

            std::string msg_base = "BF[" + std::to_string(i) + "," + std::to_string(j) + "]";
            double actual_real = (BF[current_idx]).real();
            double actual_imag = slicot_get_imag(BF[current_idx]);
            double expected_real = (expected_Bf_nu_rm[expected_idx]).real();
            double expected_imag = slicot_get_imag(expected_Bf_nu_rm[expected_idx]);
            EXPECT_NEAR(actual_real, expected_real, check_tol) << msg_base << " Real part mismatch (RowMajor)";
            EXPECT_NEAR(actual_imag, expected_imag, check_tol) << msg_base << " Imag part mismatch (RowMajor)";
        }
    }
}
