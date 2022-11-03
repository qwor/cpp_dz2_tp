#include <iostream>

#include <Vector.h>
#include <Matrix.h>

int main() {
  constexpr int n = 2;
  Matrix<double, n, n> m(1);
  Vector<double, 3> v;
  std::cout << v;
  std::cout << m.Det();

  Matrix<double, 2, 4> m2;
  std::cout << m2.Transpose();
}
