#include <iostream>

#include "gtest/gtest.h"

#include <Matrix.h>
#include <Vector.h>

TEST(VectorTest, CreateFromNumbers) {
  // вектор, заполненный тройками
  Vector<int, 3> vec_a(3, 3);

  Vector<int, 3> vec_b;
  for (std::size_t i = 0; i < vec_b.size(); ++i) {
    vec_b[i] = 3;
  }

  EXPECT_EQ(vec_a, vec_b);
}

TEST(VectorTest, CreateFromArray) {
  int arr[] { 1, 2, 3 };

  Vector<int, 3> vec_a(arr);
  Vector<int, 3> vec_b;
  vec_b[0] = 1;
  vec_b[1] = 2;
  vec_b[2] = 3;

  EXPECT_EQ(vec_a, vec_b);
}

TEST(VectorTest, TemplateParams) {
  Vector<double, 3> vec;

  EXPECT_EQ(vec.size(), 3);
  EXPECT_EQ(typeid(vec.Get(0)), typeid(double));
}

TEST(VectorTest, EmptyTemplateParams) {
  Vector<double> vec;
  Vector<double, 0> vec_example;

  EXPECT_EQ(vec.size(), 0);
  EXPECT_EQ(typeid(vec), typeid(vec_example));
}

TEST(VectorTest, Size) {
  Vector<double, 2> vec_b(3);
  Vector<double> vec_c(3);

  EXPECT_EQ(vec_b.size(), 3);
  EXPECT_EQ(vec_c.size(), 3);
}

TEST(VectorTest, OperatorGetSet) {
  Vector<int, 10> vec;
  for (std::size_t i = 0; i < vec.size(); ++i) {
    vec[i] = static_cast<int>(i);
  }
  EXPECT_EQ(vec[1], 1);
  vec[1]++;
  EXPECT_EQ(vec[1], 2);  // Значение должно измениться
}

TEST(MatrixTest, GetVal) {
  Vector<int, 10> vec;
  for (std::size_t i = 0; i < vec.size(); ++i) {
    vec[i] = static_cast<int>(i);
  }

  EXPECT_EQ(vec.Get(1), 1);
}

TEST(VectorTest, CrossProduct) {
  Vector<int, 4> vec_a;
  Vector<int, 5> vec_b;
  for (int i = 0; i < 4; i++) {
    vec_a[i] = i + 1;
    vec_b[i] = i + 2;
  }
  vec_b[4] = 6;
  auto m = vec_a.CrossProduct(vec_b);

  EXPECT_EQ(m.rows(), vec_a.size());
  EXPECT_EQ(m.cols(), vec_b.size());
  EXPECT_EQ(m(3, 3), 20);
}

TEST(VectorTest, DotProduct) {
  Vector<int, 4> vec_a, vec_b;
  for (int i = 0; i < 4; i++) {
    vec_a[i] = i + 1;
    vec_b[i] = i + 2;
  }
  int dot = vec_a.DotProduct(vec_b);

  EXPECT_EQ(dot, 40);
}

TEST(VectorTest, Slice) {
  Vector<int, 4> vec_a;
  for (std::size_t i = 0; i < vec_a.size(); ++i) {
    vec_a[i] = static_cast<int>(i);
  }
  auto vec_slice = vec_a.Slice(1, 3);

  Vector<int> vec_slice_example(2);
  vec_slice_example[0] = 1;
  vec_slice_example[1] = 2;

  EXPECT_EQ(vec_slice, vec_slice_example);
}

TEST(VectorTest, AddValue) {
  Vector<int, 4> vec_a;
  auto vec_b = vec_a + 1;

  EXPECT_EQ(vec_b.Sum(), 4);
}

TEST(VectorTest, SubValue) {
  Vector<int, 4> vec_a;
  vec_a.Fill(2);
  auto vec_b = vec_a - 1;

  EXPECT_EQ(vec_b.Sum(), 4);
}

TEST(VectorTest, MulValue) {
  Vector<int, 4> vec_a;
  vec_a.Fill(1);
  auto vec_b = vec_a * 2;

  EXPECT_EQ(vec_b.Sum(), 8);
}

TEST(VectorTest, DivValue) {
  Vector<int, 4> vec_a;
  vec_a.Fill(8);
  auto vec_b = vec_a / 2;

  EXPECT_EQ(vec_b.Sum(), 16);
}

TEST(VectorTest, AddAssignValue) {
  Vector<int, 4> vec;
  vec += 1;

  EXPECT_EQ(vec.Sum(), 4);
}

TEST(VectorTest, SubAssignValue) {
  Vector<int, 4> vec;
  vec.Fill(2);
  vec-= 1;

  EXPECT_EQ(vec.Sum(), 4);
}

TEST(VectorTest, MulAssignValue) {
  Vector<int, 4> vec;
  vec.Fill(1);
  vec *= 2;

  EXPECT_EQ(vec.Sum(), 8);
}

TEST(VectorTest, DivAssignValue) {
  Vector<int, 4> vec_a;
  vec_a.Fill(8);
  vec_a /= 2;

  EXPECT_EQ(vec_a.Sum(), 16);
}

TEST(VectorTest, AddVector) {
  Vector<int, 4> vec_a;
  Vector<int, 4> vec_b;
  vec_b.Fill(2);
  auto vec_c = vec_a + vec_b;

  EXPECT_EQ(vec_c.Sum(), 8);
}

TEST(VectorTest, SubVector) {
  Vector<int, 4> vec_a;
  vec_a.Fill(3);
  Vector<int, 4> vec_b;
  vec_b.Fill(2);
  auto vec_c = vec_a - vec_b;

  EXPECT_EQ(vec_c.Sum(), 4);
}

TEST(VectorTest, MulVector) {
  Vector<int, 4> vec_a;
  vec_a.Fill(1);
  Vector<int, 4> vec_b;
  vec_b.Fill(2);
  auto vec_c = vec_a * vec_b;

  EXPECT_EQ(vec_c.Sum(), 8);
}

TEST(VectorTest, DivVector) {
  Vector<int, 4> vec_a;
  vec_a.Fill(8);
  Vector<int, 4> vec_b;
  vec_b.Fill(4);
  auto vec_c = vec_a / vec_b;

  EXPECT_EQ(vec_c.Sum(), 8);
}

TEST(VectorTest, AddAssignVector) {
  Vector<int, 4> vec_a;
  Vector<int, 4> vec_b;
  vec_b.Fill(2);
  vec_a += vec_b;

  EXPECT_EQ(vec_a.Sum(), 8);
}

TEST(VectorTest, SubAssignVector) {
  Vector<int, 4> vec_a;
  vec_a.Fill(3);
  Vector<int, 4> vec_b;
  vec_b.Fill(2);
  vec_a -= vec_b;

  EXPECT_EQ(vec_a.Sum(), 4);
}

TEST(VectorTest, MulAssignVector) {
  Vector<int, 4> vec_a;
  vec_a.Fill(1);
  Vector<int, 4> vec_b;
  vec_b.Fill(2);
  vec_a *= vec_b;

  EXPECT_EQ(vec_a.Sum(), 8);
}

TEST(VectorTest, DivAssignVector) {
  Vector<int, 4> vec_a;
  vec_a.Fill(8);
  Vector<int, 4> vec_b;
  vec_b.Fill(4);
  vec_a /= vec_b;

  EXPECT_EQ(vec_a.Sum(), 8);
}
