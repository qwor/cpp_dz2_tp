#include <iostream>

#include <Vector.h>
#include <Matrix.h>

int main() {
  Matrix<double, 2> m1(1);
  Matrix<double, 2> m2(2);
  auto m3 = m1 + m2;
  std::cout << m1 << std::endl << m2 << std::endl << m3;
}
