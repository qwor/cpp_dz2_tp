#include <iostream>
#include <array>

#include <Vector.h>
#include <Matrix.h>

int main() {
  constexpr int n = 2;
  Matrix<double, n, n> m;
  for (std::size_t i = 0; i < n; i++) {
    m.Set(i, i, 1);
  }
  Vector<double, 3> v;
  std::cout << m.Det();
  //auto inverse = m.Inverse();
  //auto ident = m.MatrixProduct(inverse);

}
