#include <iostream>
#include <array>

#include " matrix_lib/Vector.h"
#include " matrix_lib/Matrix.h"

int main() {
  constexpr int n = 4;
  Matrix<double, n, n> m;
  int k = 0;
  for (std::size_t i = 0; i < 4; i++) {
    for (std::size_t j = 0; j < 4; j++) {
      m.Set(i, j, k++);
    }
  }
  EXPECT_EQ(m.Det(), );
  return 0;
}
