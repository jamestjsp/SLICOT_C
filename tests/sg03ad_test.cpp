#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm> // For std::max, std::min
#include <iostream> 
#include <iomanip>   

#include "sg03ad.h"       
#include "slicot_utils.h" 

// --- Test Fixture Base ---
class SG03ADTest : public ::testing::Test {
protected:
    // Common parameters, can be overridden in specific tests
    char DICO = 'C';
    char JOB = 'B';
    char FACT = 'N';
    char TRANS = 'N';
    char UPLO = 'L'; // Python test_sg03ad_b1 uses 'L' for Y input. HTML example uses 'U'.
    int N_param = 0;

    // Input matrices (row-major)
    std::vector<double> A_data, E_data, Q_data_in, Z_data_in, Y_data;
    // Output matrices/scalars
    std::vector<double> A_out, E_out, Q_out, Z_out, X_out;
    double SCALE_out, SEP_out, FERR_out;
    std::vector<double> ALPHAR_out, ALFAI_out, BETA_out;

    // Expected results
    std::vector<double> X_expected_rm;
    int expected_info = 0;

    // Leading dimensions (will be set based on N_param and row_major in tests)
    int LDA, LDE, LDQ, LDZ, LDX;
    
    double check_tol_X = 1e-5; // Default tolerance for X comparison

    void print_rm_matrix(const std::string& name, const std::vector<double>& mat, int r, int c, int ld_rm_cols) {
        std::cout << name << " (" << r << "x" << c << ", ld=" << ld_rm_cols << ") row-major:\n";
        if (mat.empty() && r > 0 && c > 0) { std::cout << "  (empty)\n"; return; }
        if (r == 0 || c == 0) { std::cout << "  (zero dim)\n"; return; }
        for (int i = 0; i < r; ++i) {
            std::cout << "  [";
            for (int j = 0; j < c; ++j) {
                std::cout << std::fixed << std::setprecision(8) << mat[i * ld_rm_cols + j] << (j == c - 1 ? "" : ", ");
            }
            std::cout << "]\n";
        }
    }


    void setup_common_buffers(int n) {
        N_param = n;
        if (n > 0) {
            A_out.resize((size_t)n*n); E_out.resize((size_t)n*n);
            Q_out.resize((size_t)n*n); Z_out.resize((size_t)n*n);
            X_out.resize((size_t)n*n);
            ALPHAR_out.resize(n); ALFAI_out.resize(n); BETA_out.resize(n);
        } else { // N=0
            A_out.assign(1,0.0); E_out.assign(1,0.0); Q_out.assign(1,0.0); 
            Z_out.assign(1,0.0); X_out.assign(1,0.0);
            ALPHAR_out.assign(1,0.0); ALFAI_out.assign(1,0.0); BETA_out.assign(1,0.0);
        }
    }
};


// Test based on SLICOT Documentation Example for SG03AD
TEST_F(SG03ADTest, SlicotDocExample) {
    DICO = 'C'; JOB = 'B'; FACT = 'N'; TRANS = 'N'; UPLO = 'U'; // HTML example uses UPLO='U' for Y
    N_param = 3;
    setup_common_buffers(N_param);

    A_data = {3.0, 1.0, 1.0, 1.0, 3.0, 0.0, 1.0, 0.0, 2.0};
    E_data = {1.0, 3.0, 0.0, 3.0, 2.0, 1.0, 1.0, 0.0, 1.0};
    // Y (input to X parameter) from example data. It's upper triangular.
    Y_data = {-64.0, -73.0, -28.0, 0.0, -70.0, -25.0, 0.0, 0.0, -18.0}; 
    
    // Q_data_in and Z_data_in are not strictly input if FACT='N', but buffers must be passed
    Q_data_in.assign((size_t)N_param*N_param, 0.0); 
    Z_data_in.assign((size_t)N_param*N_param, 0.0);

    X_expected_rm = {-2.0000, -1.0000,  0.0000,
                     -1.0000, -3.0000, -1.0000,
                      0.0000, -1.0000, -3.0000};
    check_tol_X = 1e-4; // Doc example output is rounded

    // For row_major = 1 (C style)
    LDA = N_param; LDE = N_param; LDQ = N_param; LDZ = N_param; LDX = N_param;
    
    // Copy input data to output buffers as they are modified
    A_out = A_data; E_out = E_data; Q_out = Q_data_in; Z_out = Z_data_in; X_out = Y_data;

    int info = slicot_sg03ad(DICO, JOB, FACT, TRANS, UPLO, N_param,
                             A_out.data(), LDA, E_out.data(), LDE,
                             Q_out.data(), LDQ, Z_out.data(), LDZ,
                             X_out.data(), LDX,
                             &SCALE_out, &SEP_out, &FERR_out,
                             ALPHAR_out.data(), ALFAI_out.data(), BETA_out.data(),
                             1 /* row_major = true */);
    ASSERT_EQ(info, 0);
    // Compare X_out with X_expected_rm
    for(size_t i=0; i < X_out.size(); ++i) {
        EXPECT_NEAR(X_out[i], X_expected_rm[i], check_tol_X) << "X_out[" << i << "]";
    }
    // Optionally check SCALE, SEP, FERR if their expected values are known
    // From example output: SEP ~ 0.29, FERR ~ 0.40D-13, SCALE ~ 0.10D+01 (means 1.0)
    EXPECT_NEAR(SCALE_out, 1.0, 1e-1); // SCALE is 0.10D+01 in example
    EXPECT_NEAR(SEP_out, 0.29, 1e-2);  // SEP is 0.29D+00
    // FERR is very small, check if it's close to zero
    EXPECT_LT(std::abs(FERR_out), 1e-12); 
}

// Test based on Python test_sg03ad_b1 (which is also SLICOT doc example B.1)
class SG03ADTestPythonB1 : public SG03ADTest {
protected:
    void SetUp() override {
        DICO = 'C'; JOB = 'B'; FACT = 'N'; TRANS = 'N'; UPLO = 'L'; // Python test uses 'L'
        N_param = 3;
        setup_common_buffers(N_param);

        A_data = {3.0, 1.0, 1.0, 1.0, 3.0, 0.0, 1.0, 0.0, 2.0};
        E_data = {1.0, 3.0, 0.0, 3.0, 2.0, 1.0, 1.0, 0.0, 1.0};
        // Y from python test (it's symmetric, -Y is passed to sg03ad)
        std::vector<double> Y_python = {64.0, 73.0, 28.0, 73.0, 70.0, 25.0, 28.0, 25.0, 18.0};
        Y_data.resize(Y_python.size());
        for(size_t i=0; i<Y_python.size(); ++i) Y_data[i] = -Y_python[i]; // Pass -Y

        Q_data_in.assign((size_t)N_param*N_param, 0.0); 
        Z_data_in.assign((size_t)N_param*N_param, 0.0);

        X_expected_rm = {-2.0000, -1.0000,  0.0000,
                         -1.0000, -3.0000, -1.0000,
                          0.0000, -1.0000, -3.0000};
        check_tol_X = 1e-5; // Python uses assert_almost_equal (default decimal=7)
    }
};

TEST_F(SG03ADTestPythonB1, PythonB1RowMajor) {
    LDA = N_param; LDE = N_param; LDQ = N_param; LDZ = N_param; LDX = N_param;
    A_out = A_data; E_out = E_data; Q_out = Q_data_in; Z_out = Z_data_in; X_out = Y_data;

    int info = slicot_sg03ad(DICO, JOB, FACT, TRANS, UPLO, N_param,
                             A_out.data(), LDA, E_out.data(), LDE,
                             Q_out.data(), LDQ, Z_out.data(), LDZ,
                             X_out.data(), LDX,
                             &SCALE_out, &SEP_out, &FERR_out,
                             ALPHAR_out.data(), ALFAI_out.data(), BETA_out.data(),
                             1 /* row_major = true */);
    ASSERT_EQ(info, 0);
    for(size_t i=0; i < X_out.size(); ++i) {
        EXPECT_NEAR(X_out[i], X_expected_rm[i], check_tol_X) << "X_out[" << i << "]";
    }
}

TEST_F(SG03ADTestPythonB1, PythonB1ColMajor) {
    LDA = N_param; LDE = N_param; LDQ = N_param; LDZ = N_param; LDX = N_param;
    // Prepare CM inputs
    std::vector<double> A_cm(A_data.size()), E_cm(E_data.size()), Y_cm(Y_data.size());
    std::vector<double> Q_cm_in(Q_data_in.size()), Z_cm_in(Z_data_in.size());
    slicot_transpose_to_fortran_with_ld(A_data.data(),A_cm.data(),N_param,N_param,N_param,LDA,sizeof(double));
    slicot_transpose_to_fortran_with_ld(E_data.data(),E_cm.data(),N_param,N_param,N_param,LDE,sizeof(double));
    slicot_transpose_to_fortran_with_ld(Y_data.data(),Y_cm.data(),N_param,N_param,N_param,LDX,sizeof(double));
    // Q_cm_in, Z_cm_in are outputs if FACT='N', so just need buffer
    
    // Output buffers (will be CM)
    std::vector<double> X_cm_out((size_t)N_param*N_param);
    // A_cm, E_cm, Q_cm_in, Z_cm_in will be modified by the call.
    // Use copies if original CM data needs to be preserved for other checks.
    std::vector<double> A_call = A_cm, E_call = E_cm, Q_call = Q_cm_in, Z_call = Z_cm_in;
    std::vector<double> X_call = Y_cm;


    int info = slicot_sg03ad(DICO, JOB, FACT, TRANS, UPLO, N_param,
                             A_call.data(), LDA, E_call.data(), LDE,
                             Q_call.data(), LDQ, Z_call.data(), LDZ,
                             X_call.data(), LDX, // X_call contains Y_cm on input, X_cm on output
                             &SCALE_out, &SEP_out, &FERR_out,
                             ALPHAR_out.data(), ALFAI_out.data(), BETA_out.data(),
                             0 /* row_major = false */);
    ASSERT_EQ(info, 0);
    
    // Transpose X_expected_rm to CM for comparison
    std::vector<double> X_expected_cm(X_expected_rm.size());
    slicot_transpose_to_fortran_with_ld(X_expected_rm.data(), X_expected_cm.data(), N_param,N_param,N_param,LDX,sizeof(double));

    for(size_t i=0; i < X_call.size(); ++i) { // Compare against X_call which now holds the solution
        EXPECT_NEAR(X_call[i], X_expected_cm[i], check_tol_X) << "X_cm_out[" << i << "]";
    }
}


TEST_F(SG03ADTest, ParameterValidation) {
    N_param = 1; // Minimal N for some checks
    setup_common_buffers(N_param);
    LDA=1; LDE=1; LDQ=1; LDZ=1; LDX=1;

    // Test invalid DICO
    int info = slicot_sg03ad('X', JOB, FACT, TRANS, UPLO, N_param, A_out.data(),LDA,E_out.data(),LDE,Q_out.data(),LDQ,Z_out.data(),LDZ,X_out.data(),LDX, &SCALE_out,&SEP_out,&FERR_out,ALPHAR_out.data(),ALFAI_out.data(),BETA_out.data(),1);
    EXPECT_EQ(info, -1);

    // Test invalid N
    info = slicot_sg03ad(DICO, JOB, FACT, TRANS, UPLO, -1, nullptr,1,nullptr,1,nullptr,1,nullptr,1,nullptr,1, &SCALE_out,&SEP_out,&FERR_out,nullptr,nullptr,nullptr,1);
    EXPECT_EQ(info, -6);

    // Test NULL A_io with N > 0
    info = slicot_sg03ad(DICO, JOB, FACT, TRANS, UPLO, N_param, nullptr,LDA,E_out.data(),LDE,Q_out.data(),LDQ,Z_out.data(),LDZ,X_out.data(),LDX, &SCALE_out,&SEP_out,&FERR_out,ALPHAR_out.data(),ALFAI_out.data(),BETA_out.data(),1);
    EXPECT_EQ(info, -7);
    
    // Test invalid LDA with N > 0 (row_major=true, LDA is cols)
    if (N_param > 0) {
        info = slicot_sg03ad(DICO, JOB, FACT, TRANS, UPLO, N_param, A_out.data(),0,E_out.data(),LDE,Q_out.data(),LDQ,Z_out.data(),LDZ,X_out.data(),LDX, &SCALE_out,&SEP_out,&FERR_out,ALPHAR_out.data(),ALFAI_out.data(),BETA_out.data(),1);
        EXPECT_EQ(info, -8);
    }
}

TEST_F(SG03ADTest, ZeroDimensionN) {
    N_param = 0;
    setup_common_buffers(N_param); // Will set N_param to 0, resize buffers to 1 element
    LDA=1; LDE=1; LDQ=1; LDZ=1; LDX=1; // Min leading dims

    // For N=0, many pointers can be NULL as arrays are 0x0
    int info = slicot_sg03ad('C', 'X', 'N', 'N', 'U', N_param,
                             nullptr, LDA, nullptr, LDE,
                             nullptr, LDQ, nullptr, LDZ,
                             nullptr, LDX, // Y_rhs (in X) is not needed if JOB='S', but for 'X' it is. If N=0, it's 0x0.
                             &SCALE_out, &SEP_out, &FERR_out,
                             ALPHAR_out.data(), ALFAI_out.data(), BETA_out.data(), // These are N-dim, so effectively 0-dim
                             1 /* row_major = true */);
    EXPECT_EQ(info, 0);
    EXPECT_EQ(SCALE_out, 1.0); // Typically SCALE is 1.0 if no scaling needed.
    // SEP, FERR might not be well-defined or set for N=0 if JOB='X'
}

