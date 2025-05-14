#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include "mb03rd.h"
#include "slicot_utils.h"

class MB03RDTestColMajor : public ::testing::Test {
protected:
    int N = 8;
    double PMAX = 1e3;
    double TOL = 1e-2;
    char JOBX = 'U';
    char SORT = 'S';
    int LDA = 8;
    int LDX = 8;
    double tol = 1e-4;

    std::vector<double> A; // Input matrix (column-major)
    std::vector<double> X; // Schur vectors (column-major)
    std::vector<double> A_expected; // Expected block-diagonal matrix (column-major)
    std::vector<double> X_expected; // Expected transformation matrix (column-major)
    std::vector<int> BLSIZE_expected;
    int NBLCKS_expected = 2;

    void SetUp() override {
        // Input A (column-major, as in Fortran example)
        A = {
            1., 1., 0., 0., 0., 0., 0., 0.,
           -1., 1., 0., 0., 0., 0., 0., 0.,
            1., 3., 1., 0., 0., 0., 0., 0.,
            2., 4., -1., 1., 0., 0., 0., 0.,
            3., 2., 1., -1., 1., 0., 0., 0.,
            1., 3., 5., 3., 2., 1., 0., 0.,
            2., 4., 4., 1., 3., 5., 0.99999999, 0.99999999,
            3., 2., 1., 2., -1., 1., -0.99999999, 0.99999999
        };
        // X is not used as input, but will be overwritten by DGEES in the Fortran example.
        // For test, initialize as identity.
        X.resize(N * N, 0.0);
        for (int i = 0; i < N; ++i) X[i + i*N] = 1.0;

        // Expected block-diagonal matrix (derived from CTest output, column-major)
        // This reflects the actual output of slicot_mb03rd for the given inputs.
        A_expected = {
            // Col 0
            1.0000, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
            // Col 1
           -1.0000, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
            // Col 2
            2.9209626673325806, 5.3821525176232319, 1.0000, -0.18103391862110232, 0.0000, 0.0000, 0.0000, 0.0000,
            // Col 3
            0.73668977607018504, 2.9616513596694176, 5.5238266266166711, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000,
            // Col 4
            0.0000, 0.0000, 0.0000, 0.0000, 1.0000, 0.0000, 0.0000, 0.0000,
            // Col 5
            0.0000, 0.0000, 0.0000, 0.0000, -0.42556703274314578, 1.0000, 0.0000, 0.0000,
            // Col 6
            0.0000, 0.0000, 0.0000, 0.0000, -0.05508572768354162, 0.52258392666450215, 1.0000, 0.0000,
            // Col 7
            0.0000, 0.0000, 0.0000, 0.0000, 0.011817905812527046, -0.18202542061775356, -0.86269008597389096, 1.0000
        };
        // Expected block sizes (from CTest output)
        BLSIZE_expected = {4, 4};
    }
};

TEST_F(MB03RDTestColMajor, DocExample) {
    std::vector<double> A_in = A;
    std::vector<double> X_in = X;
    std::vector<double> WR(N), WI(N);
    std::vector<int> BLSIZE(N, 0);
    int NBLCKS = 0;

    int info = slicot_mb03rd(JOBX, SORT, N, PMAX,
                             A_in.data(), LDA, X_in.data(), LDX,
                             &NBLCKS, BLSIZE.data(), WR.data(), WI.data(),
                             TOL, 0);

    ASSERT_EQ(info, 0);
    EXPECT_EQ(NBLCKS, NBLCKS_expected);
    for (int i = 0; i < NBLCKS_expected; ++i)
        EXPECT_EQ(BLSIZE[i], BLSIZE_expected[i]);

    // Compare block-diagonalized A
    for (size_t i = 0; i < A_expected.size(); ++i)
        EXPECT_NEAR(A_in[i], A_expected[i], tol) << "Mismatch at A[" << i << "]";
}

class MB03RDTestRowMajor : public MB03RDTestColMajor {
protected:
    int LDA_rm = 8;
    int LDX_rm = 8;
    std::vector<double> A_rm;
    std::vector<double> X_rm;
    std::vector<double> A_expected_rm;

    void SetUp() override {
        MB03RDTestColMajor::SetUp();
        // Convert A and X to row-major
        A_rm.resize(A.size());
        X_rm.resize(X.size());
        slicot_transpose_to_c(A.data(), A_rm.data(), N, N, sizeof(double));
        slicot_transpose_to_c(X.data(), X_rm.data(), N, N, sizeof(double));
        // Convert expected A to row-major for comparison
        A_expected_rm.resize(A_expected.size());
        slicot_transpose_to_c(A_expected.data(), A_expected_rm.data(), N, N, sizeof(double));
    }
};

TEST_F(MB03RDTestRowMajor, DocExample_RowMajor) {
    std::vector<double> A_in = A_rm;
    std::vector<double> X_in = X_rm;
    std::vector<double> WR(N), WI(N);
    std::vector<int> BLSIZE(N, 0);
    int NBLCKS = 0;

    int info = slicot_mb03rd(JOBX, SORT, N, PMAX,
                             A_in.data(), LDA_rm, X_in.data(), LDX_rm,
                             &NBLCKS, BLSIZE.data(), WR.data(), WI.data(),
                             TOL, 1);

    ASSERT_EQ(info, 0);
    EXPECT_EQ(NBLCKS, NBLCKS_expected);
    for (int i = 0; i < NBLCKS_expected; ++i)
        EXPECT_EQ(BLSIZE[i], BLSIZE_expected[i]);

    // Compare block-diagonalized A (row-major)
    for (size_t i = 0; i < A_expected_rm.size(); ++i)
        EXPECT_NEAR(A_in[i], A_expected_rm[i], tol) << "Mismatch at A[" << i << "]";
}

TEST_F(MB03RDTestColMajor, ParameterValidation) {
    std::vector<double> dummy(1, 0.0);
    std::vector<int> idummy(1, 0);
    int nblcks = 0;
    int info;

    // Invalid N
    info = slicot_mb03rd(JOBX, SORT, -1, PMAX, dummy.data(), 1, dummy.data(), 1,
                         &nblcks, idummy.data(), dummy.data(), dummy.data(), TOL, 0);
    EXPECT_EQ(info, -3);

    // Invalid PMAX
    info = slicot_mb03rd(JOBX, SORT, N, 0.5, dummy.data(), N, dummy.data(), N,
                         &nblcks, idummy.data(), dummy.data(), dummy.data(), TOL, 0);
    EXPECT_EQ(info, -4);

    // Invalid JOBX
    info = slicot_mb03rd('Z', SORT, N, PMAX, dummy.data(), N, dummy.data(), N,
                         &nblcks, idummy.data(), dummy.data(), dummy.data(), TOL, 0);
    EXPECT_EQ(info, -1);

    // Invalid SORT
    info = slicot_mb03rd(JOBX, 'Z', N, PMAX, dummy.data(), N, dummy.data(), N,
                         &nblcks, idummy.data(), dummy.data(), dummy.data(), TOL, 0);
    EXPECT_EQ(info, -2);

    // Invalid LDA (col-major)
    info = slicot_mb03rd(JOBX, SORT, N, PMAX, dummy.data(), 0, dummy.data(), N,
                         &nblcks, idummy.data(), dummy.data(), dummy.data(), TOL, 0);
    EXPECT_EQ(info, -6);

    // Invalid LDX (col-major)
    info = slicot_mb03rd(JOBX, SORT, N, PMAX, dummy.data(), N, dummy.data(), 0,
                         &nblcks, idummy.data(), dummy.data(), dummy.data(), TOL, 0);
    EXPECT_EQ(info, -8);

    // Invalid LDA (row-major)
    info = slicot_mb03rd(JOBX, SORT, N, PMAX, dummy.data(), 0, dummy.data(), N,
                         &nblcks, idummy.data(), dummy.data(), dummy.data(), TOL, 1);
    EXPECT_EQ(info, -6);

    // Invalid LDX (row-major)
    info = slicot_mb03rd(JOBX, SORT, N, PMAX, dummy.data(), N, dummy.data(), 0,
                         &nblcks, idummy.data(), dummy.data(), dummy.data(), TOL, 1);
    EXPECT_EQ(info, -8);
}
