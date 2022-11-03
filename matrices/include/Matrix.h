//
// Created by Ripku on 20.10.2022.
//

#pragma once

#include <utility>
#include <iostream>
#include <array>
#include <stdexcept>

#include "Vector.h"

template <typename T, std::size_t N, std::size_t M>
class Matrix;

template <typename T, std::size_t N, std::size_t M>
std::ostream& operator<<(std::ostream &os, const Matrix<T, N, M>& m);

template <typename T, std::size_t N = 0, std::size_t M = N>
class Matrix {
 public:
  // конструктор создания матрицы из чисел
  Matrix(T value, std::size_t rows, std::size_t cols);
  explicit Matrix(T value) : Matrix(value, N, M) {}
  Matrix(std::size_t rows, std::size_t cols) : Matrix(0, rows, cols) {}
  Matrix() : Matrix(N, M) {}

  // конструктор создания матрицы из векторов-столбцов
  Matrix(const Vector<T, N>& value, std::size_t rows, std::size_t cols);
  explicit Matrix(const Vector<T, N>& value) : Matrix(value, rows, cols) {}

  // Конструктор копирования
  Matrix(Matrix const&);
  Matrix& operator=(Matrix const& other);
  ~Matrix() { delete[] a_; }

  std::size_t rows() const { return rows_; }
  std::size_t cols() const { return cols_; }
  std::size_t size() const { return rows_ * cols_; }
  auto& UnderLyingArray() { return a_; }

  Vector<T, M> GetRow(std::size_t row) const;
  Vector<T, M> operator[](std::size_t row) const { return GetRow(row); }
  Vector<T, N> GetCol(std::size_t col) const;
  Vector<T, N> GetDiag() const;
  T Get(std::size_t row, std::size_t col) const;
  T operator()(std::size_t row, std::size_t col) const { return Get(row, col); }
  T& operator()(std::size_t row, std::size_t col);
  void Set(std::size_t row, std::size_t col, T value);

  template<std::size_t M2>
  Matrix<T, N, M2> MatrixProduct(const Matrix<T, M, M2>& other);
  Matrix<T, M, N> Transpose() const;
  // Только для квадратных матриц
  T Trace() { return GetDiag().Sum(); }
  Matrix<T, N, M> Rinor(std::size_t p, std::size_t q) const;
  T Det(std::size_t n) const;
  T Det() const;
  Matrix<T, N, M> Adjugate() const;
  Matrix<double, N, M> Inverse() const;

  T Sum() const;

  static Matrix<T, N, M> GetIdentityMatrix();

  // Каждый ряд (или столбец) матрицы и вектор
  Matrix<T, N, M> AddVectorHor(const Vector<T, M>& v) { return OperateVector(Op::kAdd, v, true); }
  Matrix<T, N, M> AddVectorVert(const Vector<T, N>& v) { return OperateVector(Op::kAdd, v); }
  Matrix<T, N, M> SubVectorHor(const Vector<T, M>& v) { return OperateVector(Op::kSub, v, true); }
  Matrix<T, N, M> SubVectorVert(const Vector<T, N>& v) { return OperateVector(Op::kSub, v); }
  Matrix<T, N, M> MulVectorHor(const Vector<T, M>& v) { return OperateVector(Op::kMul, v, true); }
  Matrix<T, N, M> MulVectorVert(const Vector<T, N>& v) { return OperateVector(Op::kMul, v); }
  Matrix<T, N, M> DivVectorHor(const Vector<T, M>& v) { return OperateVector(Op::kDiv, v, true); }
  Matrix<T, N, M> DivVectorVert(const Vector<T, M>& v) { return OperateVector(Op::kDiv, v); }

  Matrix<T, N, M> operator+(T value) { return OperateValue(Op::kAdd, value); }
  Matrix<T, N, M> operator-(T value) { return OperateValue(Op::kSub, value); }
  Matrix<T, N, M> operator*(T value) { return OperateValue(Op::kMul, value); }
  Matrix<T, N, M> operator/(T value) { return OperateValue(Op::kDiv, value); }
  Matrix<T, N, M> operator+(const Vector<T, N>& v) { return AddVectorVert(v); }
  Matrix<T, N, M> operator-(const Vector<T, N>& v) { return SubVectorVert(v); }
  Matrix<T, N, M> operator*(const Vector<T, N>& v) { return MulVectorVert(v); }
  Matrix<T, N, M> operator/(const Vector<T, N>& v) { return DivVectorVert(v); }
  Matrix<T, N, M> operator+(const Matrix<T, N, M>& other) { return OperateMatrix(Op::kAdd, other); }
  Matrix<T, N, M> operator-(const Matrix<T, N, M>& other) { return OperateMatrix(Op::kSub, other); }
  Matrix<T, N, M> operator*(const Matrix<T, N, M>& other) { return OperateMatrix(Op::kMul, other); }
  Matrix<T, N, M> operator/(const Matrix<T, N, M>& other) { return OperateMatrix(Op::kDiv, other); }

  Matrix<T, N, M> operator+=(T value) { *this = OperateValue(Op::kAdd, value); return *this; }
  Matrix<T, N, M> operator-=(T value) { *this = OperateValue(Op::kSub, value); return *this; }
  Matrix<T, N, M> operator*=(T value) { *this = OperateValue(Op::kMul, value); return *this; }
  Matrix<T, N, M> operator/=(T value) { *this = OperateValue(Op::kDiv, value); return *this; }
  Matrix<T, N, M> operator+=(const Vector<T>& v) { *this = AddVectorVert(v); return *this; }
  Matrix<T, N, M> operator-=(const Vector<T>& v) { *this = SubVectorVert(v); return *this; }
  Matrix<T, N, M> operator*=(const Vector<T>& v) { *this = MulVectorVert(v); return *this; }
  Matrix<T, N, M> operator/=(const Vector<T>& v) { *this = DivVectorVert(v); return *this; }
  Matrix<T, N, M> operator+=(const Matrix<T, N, M>& other) { *this = OperateMatrix(Op::kAdd, other); return *this; }
  Matrix<T, N, M> operator-=(const Matrix<T, N, M>& other) { *this = OperateMatrix(Op::kSub, other); return *this; }
  Matrix<T, N, M> operator*=(const Matrix<T, N, M>& other) { *this = OperateMatrix(Op::kMul, other); return *this; }
  Matrix<T, N, M> operator/=(const Matrix<T, N, M>& other) { *this = OperateMatrix(Op::kDiv, other); return *this; }

  bool operator==(const Matrix<T, N, M>& other) const;

  Matrix<T> Slice(std::size_t start_rows, std::size_t end_rows, std::size_t start_cols, std::size_t end_cols) const;

  friend std::ostream & operator<< <> (std::ostream &os, const Matrix<T, N, M>& m);

 private:
  T* a_;
  std::size_t rows_;
  std::size_t cols_;

  enum Op {
    kAdd,
    kSub,
    kMul,
    kDiv
  };

  Matrix<T, N, M> OperateValue(Op op, T value);
  template<std::size_t R>
  Matrix<T, N, M> OperateVector(Op op, const Vector<T, R>& v, bool horizontal=false);
  Matrix<T, N, M> OperateMatrix(Op op, const Matrix<T, N, M>& other);

  void CheckVectorRowSize(const Vector<T, N>& v) const;
  void CheckVectorColSize(const Vector<T, M>& v) const;
  void CheckMatrixSize(const Matrix<T, N, M>& other) const;
  void CheckMatrixSquare() const;
  void CheckIndex(std::size_t row, std::size_t col) const;
};

template<typename T, std::size_t N, std::size_t M>
Matrix<T, N, M>::Matrix(T value, std::size_t rows, std::size_t cols) : rows_(rows), cols_(cols) {
  a_ = new T[size()];
  for (std::size_t i = 0; i < size(); ++i) {
    a_[i] = value;
  }
}

template<typename T, std::size_t N, std::size_t M>
Matrix<T, N, M>::Matrix(const Vector<T, N>& vec, std::size_t rows, std::size_t cols)
    : rows_(rows), cols_(cols) {
  a_ = new T[size()];
  for (std::size_t col = 0; col < cols_; ++col) {
    for (std::size_t row = 0; row < rows_; ++row) {
      a_[row * rows_ + col] = vec[row];
    }
  }
}

template<typename T, std::size_t N, std::size_t M>
Matrix<T, N, M>::Matrix(Matrix<T, N, M> const& other) : rows_(other.rows_), cols_(other.cols_) {
  a_ = new T[size()];
  for (std::size_t i = 0; i < rows_ * cols_; ++i) {
    a_[i] = other.a_[i];
  }
}

template<typename T, std::size_t N, std::size_t M>
Matrix<T, N, M>& Matrix<T, N, M>::operator=(Matrix const& other) {
  if (this == &other) {
    return *this;
  }

  rows_ = other.rows_;
  cols_ = other.cols_;
  a_ = new T[size()];
  for (std::size_t i = 0; i < rows_ * cols_; ++i) {
    a_[i] = other.a_[i];
  }
  return *this;
}

template<typename T, std::size_t N, std::size_t M>
Vector<T, M> Matrix<T, N, M>::GetRow(std::size_t row) const{
  Vector<T, M> res(cols_);
  for (std::size_t i = 0; i < cols_; ++i) {
    res[i] = a_[row * rows_ + i];
  }
  return res;
}

template<typename T, std::size_t N, std::size_t M>
Vector<T, N>  Matrix<T, N, M>::GetCol(std::size_t col) const {
  Vector<T, N> res(rows_);
  for (std::size_t i = 0; i < rows_; ++i) {
    res[i] = a_[i * rows_ + col];
  }
  return res;
}

template<typename T, std::size_t N, std::size_t M>
Vector<T, N> Matrix<T, N, M>::GetDiag() const {
  Vector<T, N> res(rows_);
  for (std::size_t i = 0; i < rows_; ++i) {
    res[i] = a_[i * rows_ + i];
  }
  return res;
}

template <typename T, std::size_t N, std::size_t M>
T Matrix<T, N, M>::Get(std::size_t row, std::size_t col) const {
  CheckIndex(row, col);
  return a_[row * rows_ + col];
}

template <typename T, std::size_t N, std::size_t M>
T& Matrix<T, N, M>::operator()(std::size_t row, std::size_t col) {
  if (row >= rows_) {
    throw std::out_of_range("Matrix row index out of range");
  }
  if (col >= cols_) {
    throw std::out_of_range("Matrix col index out of range");
  }
  return a_[row * rows_ + col];
}

template<typename T, std::size_t N, std::size_t M>
void Matrix<T, N, M>::Set(std::size_t row, std::size_t col, T value) {
  CheckIndex(row, col);
  a_[row * rows_ + col] = value;
}

template<typename T, std::size_t N, std::size_t M>
Matrix<T, M, N> Matrix<T, N, M>::Transpose() const{
  Matrix<T, M, N> t(cols_, rows_);
  for (std::size_t row = 0; row < rows_; ++row) {
    for (std::size_t col = 0; col < cols_; ++col) {
      t(col, row) = Get(row, col);
    }
  }
  return t;
}

template<typename T, std::size_t N, std::size_t M>
Matrix<T, N, M> Matrix<T, N, M>::Rinor(std::size_t p, std::size_t q) const{
  CheckMatrixSquare();
  std::size_t n = cols_;

  Matrix<T, N, M> res(n-1, n-1);
  std::size_t i = 0, j = 0;

  // Цикл проходит по всем элементам в матрице
  for (std::size_t row = 0; row < n; ++row) {
    for (std::size_t col = 0; col < n; ++col)
    {
      //  Копирует в res только те элементы, которые не находятся в данном ряде или столбце
      if (row != p && col != q)
      {
        res.Set(i, j++, Get(row, col));

        // Ряд заполнен, поэтому инкрементируем индекс ряда и обнуляем индекс столбца
        if (j == n - 1)
        {
          j = 0;
          ++i;
        }
      }
    }
  }
  return res;
}

template<typename T, std::size_t N, std::size_t M>
T Matrix<T, N, M>::Det(std::size_t n) const {
  if (n == 0) {
    return T();
  } else if (n == 1) {
    return a_[0];
  }

  int sign = 1;
  T det = 0;  // Результат
  for (std::size_t i = 0; i < n; ++i) {
    // Вычисляем матрицу дополнений для a_[0][i]
    auto minor = Rinor(0, i);
    det += sign * Get(0, i) * minor.Det(n-1);
    sign = -sign;
  }
  return det;
}

template<typename T, std::size_t N, std::size_t M>
T Matrix<T, N, M>::Det() const {
  CheckMatrixSquare();
  return Det(cols_);
}

template<typename T, std::size_t N, std::size_t M>
Matrix<T, N, M> Matrix<T, N, M>::Adjugate() const {
  CheckMatrixSquare();
  Matrix<T, N, M> res(cols_, cols_);
  std::size_t n = cols_;
  if (n == 1) {
    res.Set(0, 0, 1);
    return res;
  }

  int sign = 1;

  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
       // temp нужен для хранянения матрицы дополнений
      auto temp = Rinor(i, j);

      // Знак положительный если сумма ряда и столбца чётна
      sign = ((i + j) % 2 == 0) ? 1 : -1;

      //  Меняем ряды и столбцы, чтобы получить транспонированную матрицу дополнений
      T value = (sign) * temp.Det(n - 1);
      res.Set(j, i, value);
    }
  }
  return res;
}

template<typename T, std::size_t N, std::size_t M>
Matrix<double, N, M> Matrix<T, N, M>::Inverse() const {
  T det = Det();
  if (det == 0) {
    throw std::logic_error("Determinant must not be 0 if you want to find the inverse of a matrix");
  }

  // Находим матрицу дополнений
  auto adj = Adjugate();

  std::size_t n = cols_;
  Matrix<double, M, M> inverse(n, n);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      inverse.Set(i, j, adj[i][j] / det);
    }
  }

  return inverse;
}

template<typename T, std::size_t N, std::size_t M>
T Matrix<T, N, M>::Sum() const {
  T sum = 0;
  for (std::size_t i = 0; i < rows_*cols_; ++i) {
    sum += a_[i];
  }
  return sum;
}

template<typename T, std::size_t N, std::size_t M>
Matrix<T> Matrix<T, N, M>::Slice(std::size_t start_rows, std::size_t end_rows,
                            std::size_t start_cols, std::size_t end_cols) const {
  Matrix<T> res(end_rows - start_rows, end_cols - start_cols);
  for (std::size_t row = 0; row < (end_rows - start_rows); ++row) {
    for (std::size_t col = 0; col < (end_cols - start_cols); ++col) {
      T temp = Get(row + start_rows, col + start_cols);
      res.Set(row, col, temp);
    }
  }
  return res;
}

template<typename T, std::size_t N, std::size_t M>
template<std::size_t M2>
Matrix<T, N, M2> Matrix<T, N, M>::MatrixProduct(const Matrix<T, M, M2> &other) {
  if (cols_ != other.rows()) {
    throw std::invalid_argument("Invalid matrix argument proportions, cannot compute product");
  }
  Matrix<T, N, M2> res(rows_, other.cols());

  for (std::size_t i = 0; i < rows_; ++i) {
    for (std::size_t j = 0; j < other.cols(); ++j) {
      T sum = T();
      for (std::size_t k = 0; k < cols_; ++k) {
        sum += Get(i, k) * other.Get(k, j);
      }
      res.Set(i, j, sum);
    }
  }
  return res;
}

template<typename T, std::size_t N, std::size_t M>
Matrix<T, N, M> Matrix<T, N, M>::OperateValue(Op op, T value) {
  Matrix<T, N, M> res(rows_, cols_);
  for (std::size_t row = 0; row < rows_; ++row) {
    for (std::size_t col = 0; col < cols_; ++col) {
      switch (op) {
        case Op::kAdd:
          res(row, col) = a_[row*rows_ + col] + value;
          break;
        case Op::kSub:
          res(row, col) = a_[row*rows_ + col] - value;
          break;
        case Op::kMul:
          res(row, col) = a_[row*rows_ + col] * value;
          break;
        case Op::kDiv:
          res(row, col) = a_[row*rows_ + col] / value;
          break;
      }
    }
  }
  return res;
}

template<typename T, std::size_t N, std::size_t M>
template<std::size_t R>
Matrix<T, N, M> Matrix<T, N, M>::OperateVector(Matrix::Op op, const Vector<T, R> &v, bool horizontal) {
  if (horizontal) {
    CheckVectorColSize(v);
  } else {
    CheckVectorColSize(v);
  }
  Matrix<T, N, M> res(rows_, cols_);
  for (std::size_t row = 0; row < rows_; ++row) {
    for (std::size_t col = 0; col < cols_; ++col) {
      auto value = v[row];
      if (horizontal) {
        value = v[col];
      }
      switch (op) {
        case Op::kAdd:
          res(row, col) = a_[row*rows_ + col] + value;
          break;
        case Op::kSub:
          res(row, col) = a_[row*rows_ + col] - value;
          break;
        case Op::kMul:
          res(row, col) = a_[row*rows_ + col] * value;
          break;
        case Op::kDiv:
          res(row, col) = a_[row*rows_ + col] / value;
          break;
      }
    }
  }
  return res;
}

template<typename T, std::size_t N, std::size_t M>
Matrix<T, N, M> Matrix<T, N, M>::OperateMatrix(Op op, const Matrix<T, N, M>& other) {
  CheckMatrixSize(other);
  Matrix<T, N, M> res(rows_, cols_);
  for (std::size_t row = 0; row < rows_; ++row) {
    for (std::size_t col = 0; col < cols_; ++col) {
      switch (op) {
        case Op::kAdd:
          res(row, col) = a_[row * rows_ + col] + other.a_[row * rows_ + col];
          break;
        case Op::kSub:
          res(row, col) = a_[row * rows_ + col] - other.a_[row * rows_ + col];
          break;
        case Op::kMul:
          res(row, col) = a_[row * rows_ + col] * other.a_[row * rows_ + col];
          break;
        case Op::kDiv:
          res(row, col) = a_[row * rows_ + col] / other.a_[row * rows_ + col];
          break;
      }
    }
  }
  return res;
}

template<typename T, std::size_t N, std::size_t M>
void Matrix<T, N, M>::CheckVectorRowSize(const Vector<T, N> &v) const {
  if (rows_ != v.size()) {
    throw std::invalid_argument("Vector has wrong size, it is not equal to the number of rows in the matrix");
  }
}

template<typename T, std::size_t N, std::size_t M>
void Matrix<T, N, M>::CheckVectorColSize(const Vector<T, M> &v) const {
  if (cols_ != v.size()) {
    throw std::invalid_argument("Vector has wrong size, it is not equal to the number of cols in the matrix");
  }
}

template<typename T, std::size_t N, std::size_t M>
void Matrix<T, N, M>::CheckMatrixSize(const Matrix<T, N, M> &other) const {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::invalid_argument("Matrix dimensions must agree");
  }
}

template<typename T, std::size_t N, std::size_t M>
void Matrix<T, N, M>::CheckMatrixSquare() const {
  if (rows_ != cols_ || rows_ == 0) {
    throw std::logic_error("Matrix must be square to perform this operation");
  }
}

template<typename T, std::size_t N, std::size_t M>
void Matrix<T, N, M>::CheckIndex(std::size_t row, std::size_t col) const {
  if (row >= rows_) {
    throw std::out_of_range("Matrix row index out of range");
  }
  if (col >= cols_) {
    throw std::out_of_range("Matrix col index out of range");
  }
}

template<typename T, std::size_t N, std::size_t M>
Matrix<T, N, M> Matrix<T, N, M>::GetIdentityMatrix() {
  Matrix<T, N, M> res;
  res.CheckMatrixSquare();
  for (std::size_t i = 0; i < res.cols_; ++i) {
    res(i, i) = static_cast<T>(1);
  }
  return res;
}

template<typename T, std::size_t N, std::size_t M>
bool Matrix<T, N, M>::operator==(const Matrix<T, N, M> &other) const {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    return false;
  }
  for (std::size_t row = 0; row < rows_; ++row) {
    for (std::size_t col = 0; col < cols_; ++col) {
      if (a_[row * rows_ + col] != other.a_[row * rows_ + col]) {
        return false;
      }
    }
  }
  return true;
}

template<typename T, std::size_t N, std::size_t M>
std::ostream & operator<<(std::ostream &os, const Matrix<T, N, M>& m) {
  for (std::size_t row = 0; row < m.rows_; ++row) {
    for (std::size_t col = 0; col < m.cols_; ++col) {
      os << m(row, col) << '\t';
    }
    os << std::endl;
  }
  return os;
}
