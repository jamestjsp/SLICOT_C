#include <gtest/gtest.h>
#include <vector>

#include "sb04nd.h"
#include "slicot_utils.h"

namespace {

constexpr int kN = 5;
constexpr int kM = 3;

const std::vector<double> kAExampleColMajor = {
    17.0, 23.0,  0.0,  0.0, 0.0,
    24.0,  5.0,  6.0,  0.0, 0.0,
     1.0,  7.0, 13.0, 19.0, 0.0,
     8.0, 14.0, 20.0, 21.0, 2.0,
    15.0, 16.0, 22.0,  3.0, 9.0
};

const std::vector<double> kBExampleColMajor = {
    8.0, 0.0, 0.0,
    1.0, 5.0, 9.0,
    6.0, 7.0, 2.0
};

const std::vector<double> kCExampleColMajor = {
    62.0, 59.0, 70.0, 35.0, 36.0,
   -12.0,-10.0, -6.0, 31.0,-15.0,
    26.0, 31.0,  9.0, -7.0,  7.0
};

const std::vector<double> kXExpectedColMajor = {
     0.0, 1.0, 0.0, 1.0, 2.0,
     0.0, 0.0, 1.0, 1.0,-2.0,
     1.0, 0.0, 0.0,-1.0, 1.0
};

void CheckColumnMajor(const std::vector<double>& data,
                       const std::vector<double>& expected,
                       int rows, int cols, int ld, double tol)
{
    for (int j = 0; j < cols; ++j) {
        for (int i = 0; i < rows; ++i) {
            EXPECT_NEAR(data[j * ld + i], expected[j * rows + i], tol);
        }
    }
}

void CheckRowMajor(const std::vector<double>& data,
                    const std::vector<double>& expected_col_major,
                    int rows, int cols, int ld, double tol)
{
    std::vector<double> expected_rm((size_t)rows * (size_t)cols);
    slicot_transpose_to_c_with_ld(expected_col_major.data(), expected_rm.data(),
                                  rows, cols, rows, cols, sizeof(double));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            EXPECT_NEAR(data[i * ld + j], expected_rm[i * cols + j], tol);
        }
    }
}

} // namespace

TEST(SB04ND, DocExampleColumnMajor)
{
    std::vector<double> a = kAExampleColMajor;
    std::vector<double> b = kBExampleColMajor;
    std::vector<double> c = kCExampleColMajor;
    double tol = 0.0;

    int info = slicot_sb04nd('B', 'U', 'U',
                              kN, kM,
                              a.data(), kN,
                              b.data(), kM,
                              c.data(), kN,
                              tol,
                              0);

    ASSERT_EQ(info, 0);
    CheckColumnMajor(c, kXExpectedColMajor, kN, kM, kN, 1e-10);
}

TEST(SB04ND, DocExampleRowMajor)
{
    int lda = kN;
    int ldb = kM;
    int ldc = kM;
    std::vector<double> a((size_t)kN * (size_t)lda);
    std::vector<double> b((size_t)kM * (size_t)ldb);
    std::vector<double> c((size_t)kN * (size_t)ldc);

    slicot_transpose_to_c_with_ld(kAExampleColMajor.data(), a.data(), kN, kN, kN, lda, sizeof(double));
    slicot_transpose_to_c_with_ld(kBExampleColMajor.data(), b.data(), kM, kM, kM, ldb, sizeof(double));
    slicot_transpose_to_c_with_ld(kCExampleColMajor.data(), c.data(), kN, kM, kN, ldc, sizeof(double));

    double tol = 0.0;
    int info = slicot_sb04nd('B', 'U', 'U',
                              kN, kM,
                              a.data(), lda,
                              b.data(), ldb,
                              c.data(), ldc,
                              tol,
                              1);

    ASSERT_EQ(info, 0);
    CheckRowMajor(c, kXExpectedColMajor, kN, kM, ldc, 1e-10);
}

TEST(SB04ND, ZeroDimensionN)
{
    int n = 0;
    int m = 3;
    std::vector<double> b = kBExampleColMajor;
    int info = slicot_sb04nd('S', 'U', 'U',
                              n, m,
                              nullptr, 1,
                              b.data(), m,
                              nullptr, 1,
                              0.0,
                              0);
    EXPECT_EQ(info, 0);
}

TEST(SB04ND, ZeroDimensionM)
{
    int n = 3;
    int m = 0;
    std::vector<double> a = {1.0, 0.0, 0.0,
                             0.0, 1.0, 0.0,
                             0.0, 0.0, 1.0};
    int info = slicot_sb04nd('S', 'U', 'U',
                              n, m,
                              a.data(), n,
                              nullptr, 1,
                              nullptr, n,
                              0.0,
                              0);
    EXPECT_EQ(info, 0);
}

TEST(SB04ND, ZeroDimensionBoth)
{
    int info = slicot_sb04nd('S', 'U', 'U',
                              0, 0,
                              nullptr, 1,
                              nullptr, 1,
                              nullptr, 1,
                              0.0,
                              0);
    EXPECT_EQ(info, 0);
}

TEST(SB04ND, ParameterValidation)
{
    std::vector<double> a = kAExampleColMajor;
    std::vector<double> b = kBExampleColMajor;
    std::vector<double> c = kCExampleColMajor;

    EXPECT_EQ(slicot_sb04nd('X', 'U', 'U', kN, kM,
                             a.data(), kN, b.data(), kM, c.data(), kN, 0.0, 0), -1);

    EXPECT_EQ(slicot_sb04nd('B', 'X', 'U', kN, kM,
                             a.data(), kN, b.data(), kM, c.data(), kN, 0.0, 0), -2);

    EXPECT_EQ(slicot_sb04nd('B', 'U', 'X', kN, kM,
                             a.data(), kN, b.data(), kM, c.data(), kN, 0.0, 0), -3);

    EXPECT_EQ(slicot_sb04nd('B', 'U', 'U', -1, kM,
                             a.data(), kN, b.data(), kM, c.data(), kN, 0.0, 0), -4);

    EXPECT_EQ(slicot_sb04nd('B', 'U', 'U', kN, -1,
                             a.data(), kN, b.data(), kM, c.data(), kN, 0.0, 0), -5);

    EXPECT_EQ(slicot_sb04nd('B', 'U', 'U', kN, kM,
                             a.data(), kN - 1, b.data(), kM, c.data(), kN, 0.0, 0), -7);

    EXPECT_EQ(slicot_sb04nd('B', 'U', 'U', kN, kM,
                             a.data(), kN, b.data(), kM - 1, c.data(), kN, 0.0, 0), -9);

    EXPECT_EQ(slicot_sb04nd('B', 'U', 'U', kN, kM,
                             a.data(), kN, b.data(), kM, c.data(), kN - 1, 0.0, 0), -11);
}
