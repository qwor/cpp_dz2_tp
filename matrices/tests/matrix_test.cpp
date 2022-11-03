#include <iostream>

#include "gtest/gtest.h"

#include <Vector.h>
#include <Matrix.h>

TEST(MatrixTest, CreateFromNumbers) {
  // матрица, заполненная тройками
  Matrix<int, 3, 3> mat_a(3);

  Matrix<int, 3, 3> mat_b;
  for (std::size_t i = 0; i < mat_b.rows(); ++i) {
    for (std::size_t j= 0; j < mat_b.cols(); ++j) {
      mat_b(i, j) = 3;
    }
  }

  EXPECT_EQ(mat_a, mat_b);
}

TEST(MatrixTest, CreateFromVector) {
  Vector<int, 3> vec(3, 3);
  Matrix<int, 3, 3> mat_a(vec);

  Matrix<int, 3, 3> mat_b;
  for (std::size_t i = 0; i < mat_b.rows(); ++i) {
    for (std::size_t j= 0; j < mat_b.cols(); ++j) {
      mat_b(i, j) = 3;
    }
  }

  EXPECT_EQ(mat_a, mat_b);
}

TEST(MatrixTest, CreateFromArray) {
  int arr[] { 1, 2, 3, 4 };

  Matrix<int, 2> mat_a(arr);
  Matrix<int, 2> mat_b;
  mat_b(0, 0) = 1;
  mat_b(0, 1) = 2;
  mat_b(1, 0) = 3;
  mat_b(1, 1) = 4;

  EXPECT_EQ(mat_a, mat_b);
}

TEST(MatrixTest, TemplateParams) {
  Matrix<double, 3, 3> mat;

  EXPECT_EQ(mat.rows(), 3);
  EXPECT_EQ(mat.cols(), 3);
  EXPECT_EQ(typeid(mat(0, 0)), typeid(double));
}

TEST(MatrixTest, HalfTemplateParams) {
  Matrix<double, 3> mat;
  Matrix<double, 3, 3> mat2;

  EXPECT_EQ(mat.rows(), 3);
  EXPECT_EQ(mat.cols(), 3);
  EXPECT_EQ(typeid(mat), typeid(mat2));
}

TEST(MatrixTest, EmptyTemplateParams) {
  Matrix<double> mat_a;
  Matrix<double, 0, 0> mat_b;

  EXPECT_EQ(mat_a.rows(), 0);
  EXPECT_EQ(mat_a.cols(), 0);
  EXPECT_EQ(typeid(mat_a), typeid(mat_b));
}

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

TEST(MatrixTest, OperatorGetSet) {
  Matrix<int, 10> mat;
  for (std::size_t i = 0; i < mat.rows(); ++i) {
    for (std::size_t j = 1; j < mat.cols(); ++j) {
      mat(i, j) = static_cast<int>(i * j);
    }
  }
  EXPECT_EQ(mat(1, 2), 2);
  mat(1, 2)++;
  EXPECT_EQ(mat(1, 2), 3);  // Значение должно измениться

  EXPECT_EQ(mat[1][2], 3);
  mat[1][2]++;
  EXPECT_EQ(mat[1][2], 3); // Значение не должно измениться
}

TEST(MatrixTest, GetSetVal) {
  Matrix<int, 10> mat;
  mat.Set(1, 2, 4);
  EXPECT_EQ(mat.Get(1, 2), 4);
}

TEST(MatrixTest, GetRow) {
  Matrix<int, 10> mat;
  for (std::size_t i = 0; i < mat.rows(); ++i) {
    for (std::size_t j = 0; j < mat.cols(); ++j) {
      mat(i, j) = static_cast<int>(i);
    }
  }
  auto row = mat.GetRow(1);
  auto row2 = mat[1];
  Vector<int, 10> example_vec;
  example_vec.Fill(1);

  EXPECT_EQ(row, example_vec);
  EXPECT_EQ(row2, example_vec);

  row += 1;
  row2 += 1;

  // Значения исходных векторов не должны меняться
  EXPECT_EQ(mat.GetRow(1), example_vec);
  EXPECT_EQ(mat[1], example_vec);
}

TEST(MatrixTest, GetCol) {
  Matrix<int, 10> mat;
  for (std::size_t i = 0; i < 10; ++i) {
    for (std::size_t j = 0; j < 10; ++j) {
      mat(i, j) = static_cast<int>(j);
    }
  }
  auto col = mat.GetCol(1);
  Vector<int, 10> example_vec;
  example_vec.Fill(1);

  EXPECT_EQ(col, example_vec);

  col += 1;

  // Значение исходного вектора не должно меняться
  EXPECT_EQ(mat.GetCol(1), example_vec);
}

TEST(MatrixTest, GetDiag) {
  Matrix<int, 10> mat;
  for (std::size_t i = 0; i < mat.rows(); ++i) {
    mat(i, i) = 1;
  }
  auto diag = mat.GetDiag();
  Vector<int, 10> example_vec;
  example_vec.Fill(1);

  EXPECT_EQ(diag, example_vec);

  diag += 1;

  // Значение исходного вектора не должно меняться
  EXPECT_EQ(mat.GetDiag(), example_vec);
}

TEST(MatrixTest, Slice) {
  Matrix<int, 4> mat;
  int k = 0;
  for (std::size_t i = 0; i < mat.rows(); ++i) {
    for (std::size_t j = 0; j < mat.cols(); ++j) {
      mat(i, j) = k++;
    }
  }

  Matrix<int> mat_slice_example(2, 2);
  mat_slice_example(0, 0) = 5;
  mat_slice_example(0, 1) = 6;
  mat_slice_example(1, 0) = 9;
  mat_slice_example(1, 1) = 10;

  auto m_slice = mat.Slice(1, 3, 1, 3);

  EXPECT_EQ(m_slice, mat_slice_example);
}

TEST(MatrixTest, Transpose) {
  Matrix<int, 2, 4> m;
  int k = 0;
  for (std::size_t i = 0; i < m.rows(); ++i) {
    for (std::size_t j = 0; j < m.cols(); ++j) {
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
  Matrix<double, 2> m;
  m(0, 0) = 5;
  m(0, 1) = 6;
  m(1, 0) = 4;
  m(1, 1) = 5;
  auto inverse = m.Inverse();
  auto ident = m.MatrixProduct(inverse);

  // Проверка на единичную матрицу
  EXPECT_EQ(ident, m.GetIdentityMatrix());
}

TEST(MatrixTest, Det) {
  Matrix<double, 2> m;
  m(0, 0) = 5;
  m(0, 1) = 6;
  m(1, 0) = 4;
  m(1, 1) = 5;

  EXPECT_EQ(m.Det(), 1);
}

TEST(MatrixTest, AddValue) {
  Matrix<int, 2> mat_a;
  auto mat_b = mat_a + 1;

  EXPECT_EQ(mat_b.Sum(), 4);
}

TEST(MatrixTest, SubValue) {
  Matrix<int, 2> mat_a(2);
  auto mat_b = mat_a - 1;

  EXPECT_EQ(mat_b.Sum(), 4);
}

TEST(MatrixTest, MulValue) {
  Matrix<int, 2> mat_a(1);
  auto mat_b = mat_a * 2;

  EXPECT_EQ(mat_b.Sum(), 8);
}

TEST(MatrixTest, DivValue) {
  Matrix<int, 2> mat_a(4);
  auto mat_b = mat_a / 4;

  EXPECT_EQ(mat_b.Sum(), 4);
}

TEST(MatrixTest, AddAssignValue) {
  Matrix<int, 2> mat;
  mat += 1;

  EXPECT_EQ(mat.Sum(), 4);
}

TEST(MatrixTest, SubAssignValue) {
  Matrix<int, 2> mat(2);
  mat -= 1;

  EXPECT_EQ(mat.Sum(), 4);
}

TEST(MatrixTest, MulAssignValue) {
  Matrix<int, 2> mat(1);
  mat *= 2;

  EXPECT_EQ(mat.Sum(), 8);
}

TEST(MatrixTest, DivAssignValue) {
  Matrix<int, 2> mat(4);
  mat /= 4;

  EXPECT_EQ(mat.Sum(), 4);
}

TEST(MatrixTest, AddColVector) {
  Matrix<int, 2> mat_a;
  Vector<int, 2> vec;
  vec[0] = 1;
  vec[1] = 2;
  auto mat_b = mat_a + vec;

  EXPECT_EQ(mat_b.GetCol(0), vec);
  EXPECT_EQ(mat_b.GetCol(1), vec);
}

TEST(MatrixTest, SubColVector) {
  Matrix<int, 2> mat_a;
  Vector<int, 2> vec;
  vec[0] = 1;
  vec[1] = 2;
  auto mat_b = mat_a - vec;

  auto minus_vec = Vector<int, 2>() - vec;

  EXPECT_EQ(mat_b.GetCol(0), minus_vec);
  EXPECT_EQ(mat_b.GetCol(1), minus_vec);
}

TEST(MatrixTest, MulColVector) {
  Matrix<int, 2> mat_a(2);
  Vector<int, 2> vec;
  vec[0] = 1;
  vec[1] = 2;
  auto mat_b = mat_a * vec;

  auto mul_vec = vec * 2;

  EXPECT_EQ(mat_b.GetCol(0), mul_vec);
  EXPECT_EQ(mat_b.GetCol(1), mul_vec);
}

TEST(MatrixTest, DivColVector) {
  Matrix<int, 2> mat_a(2);
  Vector<int, 2> vec;
  vec[0] = 1;
  vec[1] = 2;
  auto mat_b = mat_a / vec;

  auto div_vec = Vector<int, 2>(2, 2) / vec;

  EXPECT_EQ(mat_b.GetCol(0), div_vec);
  EXPECT_EQ(mat_b.GetCol(1), div_vec);
}

TEST(MatrixTest, AddAssignColVector) {
  Matrix<int, 2> mat;
  Vector<int, 2> vec;
  vec[0] = 1;
  vec[1] = 2;
  mat += vec;

  EXPECT_EQ(mat.GetCol(0), vec);
  EXPECT_EQ(mat.GetCol(1), vec);
}

TEST(MatrixTest, SubAssignColVector) {
  Matrix<int, 2> mat;
  Vector<int, 2> vec;
  vec[0] = 1;
  vec[1] = 2;
  mat -= vec;

  auto minus_vec = Vector<int, 2>() - vec;

  EXPECT_EQ(mat.GetCol(0), minus_vec);
  EXPECT_EQ(mat.GetCol(1), minus_vec);
}

TEST(MatrixTest, MulAssignColVector) {
  Matrix<int, 2> mat(2);
  Vector<int, 2> vec;
  vec[0] = 1;
  vec[1] = 2;
  mat *= vec;

  auto mul_vec = vec * 2;

  EXPECT_EQ(mat.GetCol(0), mul_vec);
  EXPECT_EQ(mat.GetCol(1), mul_vec);
}

TEST(MatrixTest, DivAssignColVector) {
  Matrix<int, 2> mat(2);
  Vector<int, 2> vec;
  vec[0] = 1;
  vec[1] = 2;
  mat /= vec;

  auto div_vec = Vector<int, 2>(2, 2) / vec;

  EXPECT_EQ(mat.GetCol(0), div_vec);
  EXPECT_EQ(mat.GetCol(1), div_vec);
}

TEST(MatrixTest, AddRowVector) {
  Matrix<int, 2> mat_a;
  Vector<int, 2> vec;
  vec[0] = 1;
  vec[1] = 2;
  auto mat_b = mat_a.AddRowVec(vec);

  EXPECT_EQ(mat_b[0], vec);
  EXPECT_EQ(mat_b[1], vec);
}

TEST(MatrixTest, SubRowVector) {
  Matrix<int, 2> mat_a;
  Vector<int, 2> vec;
  vec[0] = 1;
  vec[1] = 2;
  auto mat_b = mat_a.SubRowVec(vec);

  auto minus_vec = Vector<int, 2>() - vec;

  EXPECT_EQ(mat_b[0], minus_vec);
  EXPECT_EQ(mat_b[1], minus_vec);
}

TEST(MatrixTest, MulRowVector) {
  Matrix<int, 2> mat_a(2);
  Vector<int, 2> vec;
  vec[0] = 1;
  vec[1] = 2;
  auto mat_b = mat_a.MulRowVec(vec);

  auto mul_vec = vec * 2;

  EXPECT_EQ(mat_b[0], mul_vec);
  EXPECT_EQ(mat_b[1], mul_vec);
}

TEST(MatrixTest, DivRowVector) {
  Matrix<int, 2> mat_a(2);
  Vector<int, 2> vec;
  vec[0] = 1;
  vec[1] = 2;
  auto mat_b = mat_a.DivRowVec(vec);

  auto div_vec = Vector<int, 2>(2, 2) / vec;

  EXPECT_EQ(mat_b[0], div_vec);
  EXPECT_EQ(mat_b[1], div_vec);
}

TEST(MatrixTest, AddMatrix) {
  Matrix<int, 2> mat_a(1);
  Matrix<int, 2> mat_b(2);
  Matrix<int, 2> mat_example(3);
  auto mat_c = mat_a + mat_b;

  EXPECT_EQ(mat_c, mat_example);
}

TEST(MatrixTest, SubMatrix) {
  Matrix<int, 2> mat_a(3);
  Matrix<int, 2> mat_b(1);
  Matrix<int, 2> mat_example(2);
  auto mat_c = mat_a - mat_b;

  EXPECT_EQ(mat_c, mat_example);
}

TEST(MatrixTest, MulMatrix) {
  Matrix<int, 2> mat_a(1);
  Matrix<int, 2> mat_b(2);
  Matrix<int, 2> mat_example(2);
  auto mat_c = mat_a * mat_b;

  EXPECT_EQ(mat_c, mat_example);
}

TEST(MatrixTest, DivMatrix) {
  Matrix<int, 2> mat_a(4);
  Matrix<int, 2> mat_b(2);
  Matrix<int, 2> mat_example(2);
  auto mat_c = mat_a / mat_b;

  EXPECT_EQ(mat_c, mat_example);
}

TEST(MatrixTest, AddAssignMatrix) {
  Matrix<int, 2> mat_a(1);
  Matrix<int, 2> mat_b(2);
  Matrix<int, 2> mat_example(3);
  mat_a += mat_b;

  EXPECT_EQ(mat_a, mat_example);
}

TEST(MatrixTest, SubAssignMatrix) {
  Matrix<int, 2> mat_a(3);
  Matrix<int, 2> mat_b(1);
  Matrix<int, 2> mat_example(2);
  mat_a -= mat_b;

  EXPECT_EQ(mat_a, mat_example);
}

TEST(MatrixTest, MulAssignMatrix) {
  Matrix<int, 2> mat_a(1);
  Matrix<int, 2> mat_b(2);
  Matrix<int, 2> mat_example(2);
  mat_a *= mat_b;

  EXPECT_EQ(mat_a, mat_example);
}

TEST(MatrixTest, DivAssignMatrix) {
  Matrix<int, 2> mat_a(4);
  Matrix<int, 2> mat_b(2);
  Matrix<int, 2> mat_example(2);
  mat_a /= mat_b;

  EXPECT_EQ(mat_a, mat_example);
}
