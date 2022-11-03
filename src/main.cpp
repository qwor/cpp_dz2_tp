#include <iostream>

#include <Vector.h>
#include <Matrix.h>

int main() {
  constexpr int n = 2;
  Matrix<double, n, n> m;
  for (std::size_t i = 0; i < n; i++) {
    m.Set(i, i, 1);
  }
  Vector<double, 3> v;
  std::cout << v;
  std::cout << m;
  std::cout << m.Det();
}
