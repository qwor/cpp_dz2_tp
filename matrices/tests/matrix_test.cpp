#include <iostream>

#include "gtest/gtest.h"

#include <Vector.h>
#include <Matrix.h>

TEST(MatrixTest, Size) {
  Matrix<double, 3> mat_a;
  Matrix<double, 2, 4> mat_b(3, 3);
  Matrix<double> mat_c(3, 3);
  EXPECT_EQ(mat_a.rows(), 3);
  EXPECT_EQ(mat_a.cols(), 3);
  EXPECT_EQ(mat_b.rows(), 3);
  EXPECT_EQ(mat_b.cols(), 3);
  EXPECT_EQ(mat_c.rows(), 3);
  EXPECT_EQ(mat_c.cols(), 3);
}

TEST(MatrixTest, GetValues) {
  Matrix<int, 10> mat_a;
  for (std::size_t i = 0; i < 10; i++) {
    for (std::size_t j = 1; j < 10; j++) {
      mat_a.Set(i, j, (int)(i*j));
    }
  }
  EXPECT_EQ(mat_a[1][2], 2);
  mat_a[1][2]++;
  EXPECT_EQ(mat_a[1][2], 2);  // значение не изменяется
  auto col = mat_a.GetCol(1);
  col[2]++;
  EXPECT_EQ(mat_a.GetCol(1)[2], 2); // значение не изменяется
  EXPECT_EQ(mat_a.GetDiag().Sum(), 285);
}

TEST(MatrixTest, Slice) {
  Matrix<int, 4> mat;
  mat.Set(1, 1, 3);
  auto m2 = mat.Slice(1, 3, 1, 3);
  EXPECT_EQ(m2[0][0], 3);
  EXPECT_EQ(m2.cols(), 2);
  EXPECT_EQ(m2.rows(), 2);
}

TEST(MatrixTest, Arithmetic) {
  Matrix<int, 3ULL, 3ULL> mat_a;
  auto mat_b = (((mat_a + 5) - 1) * 2) / 8;
  EXPECT_EQ(mat_b.Sum(), 9);

  mat_a += 5;
  mat_a -= 1;
  mat_a *= 2;
  mat_a /= 8;
  EXPECT_EQ(mat_a.Sum(), 9);

  auto mat_c = (((mat_a + mat_b) - mat_b) * mat_b) / mat_b;
  EXPECT_EQ(mat_c.Sum(), 9);

  mat_a += mat_b;
  mat_a -= mat_b;
  mat_a *= mat_b;
  mat_a /= mat_b;
  EXPECT_EQ(mat_a.Sum(), 9);

  Vector<int, 3ULL> vec(3);
  vec += 2;
  mat_a = mat_a.AddVecVer(vec);
  mat_a = mat_a.SubVecVer(vec);
  mat_a = mat_a.MulVecHor(vec);
  mat_a = mat_a.DivVecHor(vec);

  EXPECT_EQ(mat_a.Sum(), 9);
}

TEST(MatrixTest, Transpose) {
  Matrix<int, 2, 4> m;
  int k = 0;
  for (std::size_t i = 0; i < 2; i++) {
    for (std::size_t j = 0; j < 4; j++) {
      m(i, j) = k++;
    }
  }
  auto m2 = m.Transpose();
  EXPECT_EQ(m2.rows(), 4);
  EXPECT_EQ(m2.cols(), 2);
  EXPECT_EQ(m2(3, 1), m(1, 3));
  EXPECT_EQ(m, m2.Transpose());
}

TEST(MatrixTest, Inverse) {
  constexpr int n = 2;
  Matrix<double, n, n> m;
  m.Set(0, 0, 5);
  m.Set(0, 1, 6);
  m.Set(1, 0, 4);
  m.Set(1, 1, 5);
  auto inverse = m.Inverse();
  auto ident = m.MatrixProduct(inverse);

  // Проверка на единичную матрицу
  EXPECT_EQ(ident, m.GetIdentityMatrix());
}

TEST(MatrixTest, Det) {
  constexpr int n = 2;
  Matrix<double, n, n> m;
  m.Set(0, 0, 5);
  m.Set(0, 1, 6);
  m.Set(1, 0, 4);
  m.Set(1, 1, 5);
  EXPECT_EQ(m.Det(), 1);
}