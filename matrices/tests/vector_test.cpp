#include <iostream>

#include "gtest/gtest.h"

#include <Matrix.h>
#include <Vector.h>

TEST(VectorTest, Size) {
  Vector<double, 3> vec_a;
  Vector<double, 2> vec_b(3);
  Vector<double> vec_c(3);
  EXPECT_EQ(vec_a.size(), vec_b.size());
  EXPECT_EQ(vec_b.size(), vec_c.size());
}

TEST(VectorTest, CrossDotProducts) {
  Vector<int, 4> vec_a, vec_b;
  for (int i = 0; i < 4; i++) {
    vec_a[i] = i + 1;
    vec_b[i] = i + 2;
  }
  auto m = vec_a.CrossProduct(vec_b);
  EXPECT_EQ(m(3, 3), 20);
  int dot = vec_a.DotProduct(vec_b);
  EXPECT_EQ(dot, 40);
}

TEST(VectorTest, Slice) {
  Vector<int, 4> vec_a;
  vec_a[1] = 3;
  auto vec_b = vec_a.Slice(1, 3);
  EXPECT_EQ(vec_b[0], 3);
  EXPECT_EQ(vec_b.size(), 2);
}

TEST(VectorTest, Arithmetic) {
  Vector<int, 4> vec_a;
  auto vec_b = (((vec_a + 5) - 1) * 2) / 8;
  EXPECT_EQ(vec_b.Sum(), 4);

  vec_a += 5;
  vec_a -= 1;
  vec_a *= 2;
  vec_a /= 8;
  EXPECT_EQ(vec_a.Sum(), 4);

  auto vec_c = (((vec_a + vec_b) - vec_b) * vec_a) / vec_b;
  EXPECT_EQ(vec_c.Sum(), 4);

  vec_c += vec_a;
  vec_c -= vec_a;
  vec_c *= vec_a;
  vec_c /= vec_a;
  EXPECT_EQ(vec_c.Sum(), 4);
}
