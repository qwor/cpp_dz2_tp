#include <iostream>

#include <Vector.h>
#include <Matrix.h>

int main() {
  Vector<int, 3> vec_a;
  for (std::size_t i = 0; i < vec_a.size(); ++i) {
    vec_a[i] = static_cast<int>(i);
  }
  Matrix<int, 4, 3> mat;
  int k = 0;
  for (std::size_t i = 0; i < mat.rows(); ++i) {
    for (std::size_t j = 0; j < mat.cols(); ++j) {
      mat(i, j) = k++;
    }
  }
  auto vec_b = vec_a * mat;
  std::cout << vec_b;
}
