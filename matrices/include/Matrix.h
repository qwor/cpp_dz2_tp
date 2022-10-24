//
// Created by Mipku on 20.10.2022.
//

#ifndef DZ2_MATRIX_H
#define DZ2_MATRIX_H


#include <utility>
#include <iostream>
#include <array>
#include <cassert>

#include "Vector.h"

template <typename T, std::size_t Rows = 0, std::size_t Cols = Rows>
class Matrix {
 public:
  enum VectorType {
    kCol,
    kRow
  };

  Matrix();
  Matrix(std::size_t rows, std::size_t cols);
  // конструктор создания матрицы из массива чисел
  explicit Matrix(T* a) : a_(a), rows_(Rows), cols_(Cols) {}
  Matrix(T* a, std::size_t rows, std::size_t cols) : a_(a), rows_(rows), cols_(cols) {}
  // конструктор создания матрицы из массива векторов-столбцов
  explicit Matrix(Vector<T>* vecs);
  Matrix(Vector<T>* vecs, std::size_t rows, std::size_t cols);

  Matrix(Matrix const&);
  Matrix& operator=(Matrix const& ref) {*this = Matrix(ref); return *this; }

  std::size_t rows() const { return rows_; }
  std::size_t cols() const { return cols_; }
  auto& UnderLyingArray() { return a_; }

  Vector<T, Cols> GetRow(std::size_t row) const;
  Vector<T, Cols> operator[](std::size_t row) const { return GetRow(row); }
  Vector<T, Rows> GetCol(std::size_t col) const;
  Vector<T, Rows> GetDiag() const;
  T Get(std::size_t row, std::size_t col) const;
  void Set(std::size_t row, std::size_t col, T value);

  template<std::size_t Cols2>
  Matrix<T, Rows, Cols2> MatrixProduct(const Matrix<T, Cols, Cols2>& other);
  Matrix<T, Cols, Rows> Transpose() const;
  // Только для квадратных матриц
  T Trace() { return GetDiag().Sum(); }
  Matrix<T, Rows, Cols> Cofactor(std::size_t p, std::size_t q) const;
  static T Det(Matrix<T, Rows, Cols> mat, std::size_t n);
  T Det() const;
  Matrix<T, Rows, Cols> Adjugate() const;
  Matrix<double, Rows, Cols> Inverse() const;

  T Sum() const;
  void Print() const;

  // Каждый ряд (или столбец) матрицы и вектор
  Matrix<T, Rows, Cols> AddVector(const Vector<T>& v, VectorType type = VectorType::kCol);
  Matrix<T, Rows, Cols> SubVector(const Vector<T>& v, VectorType type = VectorType::kCol);
  Matrix<T, Rows, Cols> MulVector(const Vector<T>& v, VectorType type = VectorType::kCol);
  Matrix<T, Rows, Cols> DivVector(const Vector<T>& v, VectorType type = VectorType::kCol);

  Matrix<T, Rows, Cols> operator+(T value);
  Matrix<T, Rows, Cols> operator-(T value);
  Matrix<T, Rows, Cols> operator*(T value);
  Matrix<T, Rows, Cols> operator/(T value);
  Matrix<T, Rows, Cols> operator+(const Vector<T>& v) { return AddVector(v); }
  Matrix<T, Rows, Cols> operator-(const Vector<T>& v) { return SubVector(v); }
  Matrix<T, Rows, Cols> operator*(const Vector<T>& v) { return MulVector(v); }
  Matrix<T, Rows, Cols> operator/(const Vector<T>& v) { return DivVector(v); }
  Matrix<T, Rows, Cols> operator+(const Matrix<T, Rows, Cols>& other);
  Matrix<T, Rows, Cols> operator-(const Matrix<T, Rows, Cols>& other);
  Matrix<T, Rows, Cols> operator*(const Matrix<T, Rows, Cols>& other);
  Matrix<T, Rows, Cols> operator/(const Matrix<T, Rows, Cols>& other);

  Matrix<T, Rows, Cols> operator+=(T value);
  Matrix<T, Rows, Cols> operator-=(T value);
  Matrix<T, Rows, Cols> operator*=(T value);
  Matrix<T, Rows, Cols> operator/=(T value);
  Matrix<T, Rows, Cols> operator+=(const Vector<T>& v) { *this = AddVector(v); return *this; }
  Matrix<T, Rows, Cols> operator-=(const Vector<T>& v) { *this = SubVector(v); return *this; }
  Matrix<T, Rows, Cols> operator*=(const Vector<T>& v) { *this = MulVector(v); return *this; }
  Matrix<T, Rows, Cols> operator/=(const Vector<T>& v) { *this = DivVector(v); return *this; }
  Matrix<T, Rows, Cols> operator+=(const Matrix<T, Rows, Cols>& other);
  Matrix<T, Rows, Cols> operator-=(const Matrix<T, Rows, Cols>& other);
  Matrix<T, Rows, Cols> operator*=(const Matrix<T, Rows, Cols>& other);
  Matrix<T, Rows, Cols> operator/=(const Matrix<T, Rows, Cols>& other);

  Matrix<T> Slice(std::size_t start_rows, std::size_t end_rows,
                              std::size_t start_cols, std::size_t end_cols) const;

 private:
  T* a_;
  std::size_t rows_;
  std::size_t cols_;
};


template<typename T, std::size_t Rows, std::size_t Cols>
Matrix<T, Rows, Cols>::Matrix() : rows_(Rows), cols_(Cols) {
  a_ = new T[rows_ * cols_];
  for (std::size_t i = 0; i < rows_ * cols_; i++) {
    a_[i] = 0;
  }
}

template<typename T, std::size_t Rows, std::size_t Cols>
Matrix<T, Rows, Cols>::Matrix(std::size_t rows, std::size_t cols) : rows_(rows), cols_(cols) {
  a_ = new T[rows_ * cols_];
  for (std::size_t i = 0; i < rows_ * cols_; i++) {
    a_[i] = 0;
  }
}

template<typename T, std::size_t Rows, std::size_t Cols>
Matrix<T, Rows, Cols>::Matrix(Matrix<T, Rows, Cols> const& other) : rows_(other.rows_), cols_(other.cols_) {
  a_ = new T[rows_ * cols_];
  for (std::size_t i = 0; i < rows_ * cols_; i++) {
    a_[i] = other.a_[i];
  }
}

template<typename T, std::size_t Rows, std::size_t Cols>
Matrix<T, Rows, Cols>::Matrix(Vector<T>* vecs, std::size_t rows, std::size_t cols) : rows_(rows), cols_(cols) {
  for (std::size_t col = 0; col < cols_; col++) {
    for (std::size_t row = 0; row < rows_; row++) {
      a_[row * Rows + col] = vecs[col][row];
    }
  }
}

template<typename T, std::size_t Rows, std::size_t Cols>
Matrix<T, Rows, Cols>::Matrix(Vector<T>* vecs) : rows_(Rows), cols_(Cols) {
  for (std::size_t col = 0; col < cols_; col++) {
    for (std::size_t row = 0; row < rows_; row++) {
      a_[row * Rows + col] = vecs[col][row];
    }
  }
}

template<typename T, std::size_t Rows, std::size_t Cols>
Vector<T, Cols> Matrix<T, Rows, Cols>::GetRow(std::size_t row) const{
  Vector<T, Cols> res(cols_);
  for (std::size_t i = 0; i < cols_; i++) {
    res[i] = a_[row * rows_ + i];
  }
  return res;
}

template<typename T, std::size_t Rows, std::size_t Cols>
Vector<T, Rows>  Matrix<T, Rows, Cols>::GetCol(std::size_t col) const {
  Vector<T, Rows> res(rows_);
  for (std::size_t i = 0; i < rows_; i++) {
    res[i] = a_[i * rows_ + col];
  }
  return res;
}

template<typename T, std::size_t Rows, std::size_t Cols>
Vector<T, Rows> Matrix<T, Rows, Cols>::GetDiag() const {
  Vector<T, Rows> res(rows_);
  for (std::size_t i = 0; i < rows_; i++) {
    res[i] = a_[i * rows_ + i];
  }
  return res;
}

template <typename T, std::size_t Rows, std::size_t Cols>
T Matrix<T, Rows, Cols>::Get(std::size_t row, std::size_t col) const {
  assert(row < rows_ and col < cols_);
  return a_[row * rows_ + col];
}

template<typename T, std::size_t Rows, std::size_t Cols>
void Matrix<T, Rows, Cols>::Set(std::size_t row, std::size_t col, T value) {
  assert(row < rows_ and col < cols_);
  a_[row * rows_ + col] = value;
}

template<typename T, std::size_t Rows, std::size_t Cols>
Matrix<T, Cols, Rows> Matrix<T, Rows, Cols>::Transpose() const{
  Matrix<T, Cols, Rows> t(cols_, rows_);
  for (std::size_t i = 0; i < rows_; i++) {
    for (std::size_t j = 0; j < cols_; j++) {
      t.Set(j, i, Get(i, j));
    }
  }
  return t;
}

template<typename T, std::size_t Rows, std::size_t Cols>
Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::Cofactor(std::size_t p, std::size_t q) const{
  assert(cols_ == rows_);
  std::size_t n = cols_;
  Matrix<T, Rows, Cols> res(cols_, cols_);
  std::size_t i = 0, j = 0;

  // Цикл проходит по всем элементам в матрице
  for (std::size_t row = 0; row < n; row++) {
    for (std::size_t col = 0; col < n; col++)
    {
      //  Копирует в res только те элементы, которые не находятся в данном ряде или столбце
      if (row != p && col != q)
      {
        res.Set(i, j++, Get(row, col));

        // Ряд заполнен, поэтому инкрементируем индекс ряда и обнуляем индекс столбца
        if (j == n - 1)
        {
          j = 0;
          i++;
        }
      }
    }
  }
  return res;
}

template<typename T, std::size_t Rows, std::size_t Cols>
T Matrix<T, Rows, Cols>::Det(Matrix<T, Rows, Cols> mat, std::size_t n) {
  int sign = 1;
  T det = 0;  // Результат

  //  Если в матрице один элемент
  if (n == 1)
    return mat.a_[0];


  for (std::size_t f = 0; f < n; f++) {
    // Вычисляем матрицу дополнений для a_[0][f]
    Matrix<T, Rows, Cols> temp = mat.Cofactor(0, f);
    det += sign * mat.a_[f] * Det(temp, n - 1);
    sign = -sign;
  }

  return det;
}

template<typename T, std::size_t Rows, std::size_t Cols>
T Matrix<T, Rows, Cols>::Det() const {
  assert(cols_ == rows_);
  return Det(*this, cols_);
}

template<typename T, std::size_t Rows, std::size_t Cols>
Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::Adjugate() const {
  assert(cols_ == rows_);
  Matrix<T, Rows, Cols> res(cols_, cols_);
  std::size_t n = cols_;
  if (n == 1) {
    res.Set(0, 0, 1);
    return res;
  }

  int sign = 1;

  for (std::size_t i = 0; i < n; i++) {
    for (std::size_t j = 0; j < n; j++) {
       // temp нужен для хранянения матрицы дополнений
      Matrix<T, Rows, Cols> temp = Cofactor(i, j);

      // Знак положительный если сумма ряда и столбца чётна
      sign = ((i + j) % 2 == 0) ? 1 : -1;

      //  Меняем ряды и столбцы, чтобы получить транспонированную матрицу дополнений
      T value = (sign) * Det(temp, n - 1);
      res.Set(j, i, value);
    }
  }
  return res;
}

template<typename T, std::size_t Rows, std::size_t Cols>
Matrix<double, Rows, Cols> Matrix<T, Rows, Cols>::Inverse() const {
  T det = Det();
  assert(det != 0);

  // Находим матрицу дополнений
  auto adj = Adjugate();

  std::size_t n = cols_;
  Matrix<double, Cols, Cols> inverse(n, n);
  for (std::size_t i = 0; i < n; i++) {
    for (std::size_t j = 0; j < n; j++) {
      inverse.Set(i, j, adj[i][j] / det);
    }
  }

  return inverse;
}

template<typename T, std::size_t Rows, std::size_t Cols>
void Matrix<T, Rows, Cols>::Print() const {
  for (std::size_t i = 0; i < rows_; i++) {
    for (std::size_t j = 0; j < cols_; j++) {
      std::cout << Get(i, j) << '\t';
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

template<typename T, std::size_t Rows, std::size_t Cols>
Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::AddVector(const Vector<T>& v, Matrix::VectorType type) {
  Matrix<T, Rows, Cols> res(rows_, cols_);
  for (std::size_t row = 0; row < rows_; row++) {
    for (std::size_t col = 0; col < cols_; col++) {
      if (type == VectorType::kCol) {
        assert(rows_ == v.size());
        res.Set(row, col, a_[row*rows_ + col] + v[row]);
      } else if (type == VectorType::kRow) {
        assert(rows_ == v.size());
        res.Set(row, col, a_[row*rows_ + col] + v[col]);
      }
    }
  }
  return res;
}

template<typename T, std::size_t Rows, std::size_t Cols>
Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::SubVector(const Vector<T>& v, Matrix::VectorType type) {
  Matrix<T, Rows, Cols> res(rows_, cols_);
  for (std::size_t row = 0; row < rows_; row++) {
    for (std::size_t col = 0; col < cols_; col++) {
      if (type == VectorType::kCol) {
        assert(rows_ == v.size());
        res.Set(row, col, a_[row*rows_ + col] - v[row]);
      } else if (type == VectorType::kRow) {
        assert(rows_ == v.size());
        res.Set(row, col, a_[row*rows_ + col] - v[col]);
      }
    }
  }
  return res;
}

template<typename T, std::size_t Rows, std::size_t Cols>
Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::MulVector(const Vector<T>& v, Matrix::VectorType type) {
  Matrix<T, Rows, Cols> res(rows_, cols_);
  for (std::size_t row = 0; row < rows_; row++) {
    for (std::size_t col = 0; col < cols_; col++) {
      if (type == VectorType::kCol) {
        assert(rows_ == v.size());
        res.Set(row, col, a_[row*rows_ + col] * v[row]);
      } else if (type == VectorType::kRow) {
        assert(rows_ == v.size());
        res.Set(row, col, a_[row*rows_ + col] * v[col]);
      }
    }
  }
  return res;
}

template<typename T, std::size_t Rows, std::size_t Cols>
Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::DivVector(const Vector<T>& v, Matrix::VectorType type) {
  Matrix<T, Rows, Cols> res(rows_, cols_);
  for (std::size_t row = 0; row < rows_; row++) {
    for (std::size_t col = 0; col < cols_; col++) {
      if (type == VectorType::kCol) {
        assert(rows_ == v.size());
        res.Set(row, col, a_[row*rows_ + col] / v[row]);
      } else if (type == VectorType::kRow) {
        assert(rows_ == v.size());
        res.Set(row, col, a_[row*rows_ + col] / v[col]);
      }
    }
  }
  return res;
}

template<typename T, std::size_t Rows, std::size_t Cols>
Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::operator+(T value) {
  Matrix<T, Rows, Cols> res(rows_, cols_);
  for (std::size_t row = 0; row < rows_; row++) {
    for (std::size_t col = 0; col < cols_; col++) {
      res.Set(row, col, a_[row*rows_ + col] + value);
    }
  }
  return res;
}

template<typename T, std::size_t Rows, std::size_t Cols>
Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::operator-(T value) {
  Matrix<T, Rows, Cols> res(rows_, cols_);
  for (std::size_t row = 0; row < rows_; row++) {
    for (std::size_t col = 0; col < cols_; col++) {
      res.Set(row, col, a_[row*rows_ + col] - value);
    }
  }
  return res;
}

template<typename T, std::size_t Rows, std::size_t Cols>
Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::operator*(T value) {
  Matrix<T, Rows, Cols> res(rows_, cols_);
  for (std::size_t row = 0; row < rows_; row++) {
    for (std::size_t col = 0; col < cols_; col++) {
      res.Set(row, col, a_[row*rows_ + col] * value);
    }
  }
  return res;
}

template<typename T, std::size_t Rows, std::size_t Cols>
Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::operator/(T value) {
  Matrix<T, Rows, Cols> res(rows_, cols_);
  for (std::size_t row = 0; row < rows_; row++) {
    for (std::size_t col = 0; col < cols_; col++) {
      res.Set(row, col, a_[row*rows_ + col] / value);
    }
  }
  return res;
}

template<typename T, std::size_t Rows, std::size_t Cols>
Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::operator+(const Matrix<T, Rows, Cols> &other) {
  assert(rows_ == other.rows_ && cols_ == other.cols_);
  Matrix<T, Rows, Cols> res(rows_, cols_);
  for (std::size_t row = 0; row < rows_; row++) {
    for (std::size_t col = 0; col < cols_; col++) {
      res.Set(row, col, a_[row*rows_ + col] + other.a_[row*rows_ + col]);
    }
  }
  return res;
}

template<typename T, std::size_t Rows, std::size_t Cols>
Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::operator-(const Matrix<T, Rows, Cols> &other) {
  assert(rows_ == other.rows_ && cols_ == other.cols_);
  Matrix<T, Rows, Cols> res(rows_, cols_);
  for (std::size_t row = 0; row < rows_; row++) {
    for (std::size_t col = 0; col < cols_; col++) {
      res.Set(row, col, a_[row*rows_ + col] - other.a_[row*rows_ + col]);
    }
  }
  return res;
}

template<typename T, std::size_t Rows, std::size_t Cols>
Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::operator*(const Matrix<T, Rows, Cols> &other) {
  assert(rows_ == other.rows_ && cols_ == other.cols_);
  Matrix<T, Rows, Cols> res(rows_, cols_);
  for (std::size_t row = 0; row < rows_; row++) {
    for (std::size_t col = 0; col < cols_; col++) {
      res.Set(row, col, a_[row*rows_ + col] * other.a_[row*rows_ + col]);
    }
  }
  return res;
}

template<typename T, std::size_t Rows, std::size_t Cols>
Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::operator/(const Matrix<T, Rows, Cols> &other) {
  assert(rows_ == other.rows_ && cols_ == other.cols_);
  Matrix<T, Rows, Cols> res(rows_, cols_);
  for (std::size_t row = 0; row < rows_; row++) {
    for (std::size_t col = 0; col < cols_; col++) {
      res.Set(row, col, a_[row*rows_ + col] / other.a_[row*rows_ + col]);
    }
  }
  return res;
}

template<typename T, std::size_t Rows, std::size_t Cols>
Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::operator+=(T value) {
  for (std::size_t row = 0; row < rows_; row++) {
    for (std::size_t col = 0; col < cols_; col++) {
      Set(row, col, a_[row*rows_ + col] + value);
    }
  }
  return *this;
}

template<typename T, std::size_t Rows, std::size_t Cols>
Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::operator-=(T value) {
  for (std::size_t row = 0; row < rows_; row++) {
    for (std::size_t col = 0; col < cols_; col++) {
      Set(row, col, a_[row*rows_ + col] - value);
    }
  }
  return *this;
}

template<typename T, std::size_t Rows, std::size_t Cols>
Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::operator*=(T value) {
  for (std::size_t row = 0; row < rows_; row++) {
    for (std::size_t col = 0; col < cols_; col++) {
      Set(row, col, a_[row*rows_ + col] * value);
    }
  }
  return *this;
}

template<typename T, std::size_t Rows, std::size_t Cols>
Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::operator/=(T value) {
  for (std::size_t row = 0; row < rows_; row++) {
    for (std::size_t col = 0; col < cols_; col++) {
      Set(row, col, a_[row*rows_ + col] / value);
    }
  }
  return *this;
}

template<typename T, std::size_t Rows, std::size_t Cols>
Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::operator+=(const Matrix<T, Rows, Cols> &other) {
  assert(rows_ == other.rows_ && cols_ == other.cols_);
  for (std::size_t row = 0; row < rows_; row++) {
    for (std::size_t col = 0; col < cols_; col++) {
      Set(row, col, a_[row*rows_ + col] + other.a_[row*rows_ + col]);
    }
  }
  return *this;
}

template<typename T, std::size_t Rows, std::size_t Cols>
Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::operator-=(const Matrix<T, Rows, Cols> &other) {
  assert(rows_ == other.rows_ && cols_ == other.cols_);
  for (std::size_t row = 0; row < rows_; row++) {
    for (std::size_t col = 0; col < cols_; col++) {
      Set(row, col, a_[row*rows_ + col] - other.a_[row*rows_ + col]);
    }
  }
  return *this;
}

template<typename T, std::size_t Rows, std::size_t Cols>
Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::operator*=(const Matrix<T, Rows, Cols> &other) {
  assert(rows_ == other.rows_ && cols_ == other.cols_);
  for (std::size_t row = 0; row < rows_; row++) {
    for (std::size_t col = 0; col < cols_; col++) {
      Set(row, col, a_[row*rows_ + col] * other.a_[row*rows_ + col]);
    }
  }
  return *this;
}

template<typename T, std::size_t Rows, std::size_t Cols>
Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::operator/=(const Matrix<T, Rows, Cols> &other) {
  assert(rows_ == other.rows_ && cols_ == other.cols_);
  for (std::size_t row = 0; row < rows_; row++) {
    for (std::size_t col = 0; col < cols_; col++) {
      Set(row, col, a_[row*rows_ + col] / other.a_[row*rows_ + col]);
    }
  }
  return *this;
}

template<typename T, std::size_t Rows, std::size_t Cols>
T Matrix<T, Rows, Cols>::Sum() const {
  T sum = 0;
  for (std::size_t i = 0; i < rows_*cols_; i++) {
    sum += a_[i];
  }
  return sum;
}

template<typename T, std::size_t Rows, std::size_t Cols>
Matrix<T> Matrix<T, Rows, Cols>::Slice(std::size_t start_rows, std::size_t end_rows,
                            std::size_t start_cols, std::size_t end_cols) const {
  Matrix<T> res(end_rows - start_rows, end_cols - start_cols);
  for (std::size_t row = 0; row < (end_rows - start_rows); row++) {
    for (std::size_t col = 0; col < (end_cols - start_cols); col++) {
      T temp = Get(row + start_rows, col + start_cols);
      res.Set(row, col, temp);
    }
  }
  return res;
}

template<typename T, std::size_t Rows, std::size_t Cols>
template<std::size_t Cols2>
Matrix<T, Rows, Cols2> Matrix<T, Rows, Cols>::MatrixProduct(const Matrix<T, Cols, Cols2> &other) {
  assert(cols_ == other.rows());
  Matrix<T, Rows, Cols2> res(rows_, other.cols());

  for (std::size_t i = 0; i < rows_; i++) {
    for (std::size_t j = 0; j < other.cols(); j++) {
      T sum = T();
      for (std::size_t k = 0; k < cols_; k++) {
        sum += Get(i, k) * other.Get(k, j);
      }
      res.Set(i, j, sum);
    }
  }
  return res;
}

#endif //DZ2_MATRIX_H
