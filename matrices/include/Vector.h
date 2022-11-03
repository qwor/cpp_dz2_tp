//
// Created by Mipku on 20.10.2022.
//

#pragma once

#include <iostream>
#include <stdexcept>
#include <cmath>

template <typename T, std::size_t Rows, std::size_t Cols>
class Matrix;

template <typename T, std::size_t N>
class Vector;

template <typename T, std::size_t N>
std::ostream & operator<<(std::ostream &os, const Vector<T, N>& v);

template <typename T, std::size_t N = 0>
class Vector {
 public:
  // конструктор создания вектора из числа
  Vector(T value, std::size_t size);
  explicit Vector(std::size_t size) : Vector(T(), size) {}
  Vector() : Vector(N) {}

  // конструктор создания вектора из массива
  Vector(T* a, std::size_t size);
  explicit Vector(T* a) : Vector(a, N) {}

  Vector(Vector const& other);
  Vector& operator=(Vector const& other);
  ~Vector() { delete[] a_; }

  std::size_t size() const { return size_; }

  template <std::size_t M>
  Matrix<T, N, M> CrossProduct(const Vector<T, M>& other) const;
  T DotProduct(const Vector<T, N>& other) const;
  Vector<T> Slice(std::size_t start, std::size_t end) const;
  Vector<T> Slice(std::size_t start) { return Slice(start, size()); };
  T* UnderlyingArray() { return a_; }
  double Magnitude() const;
  void Normalize();
  T Sum() const;
  void Fill(T value);

  T Get(std::size_t i) const { return a_[i]; }
  T operator[](std::size_t i) const { return a_[i]; };
  T& operator[](std::size_t i) { return a_[i]; }

  Vector<T, N> operator+=(T value) { return OperateAssignValue(Op::kAdd, value); }
  Vector<T, N> operator-=(T value) { return OperateAssignValue(Op::kSub, value); }
  Vector<T, N> operator*=(T value) { return OperateAssignValue(Op::kMul, value); }
  Vector<T, N> operator/=(T value) { return OperateAssignValue(Op::kDiv, value); }
  Vector<T, N> operator+=(const Vector<T, N>& other) { return OperateAssignVector(Op::kAdd, other);  }
  Vector<T, N> operator-=(const Vector<T, N>& other) { return OperateAssignVector(Op::kSub, other); }
  Vector<T, N> operator*=(const Vector<T, N>& other) { return OperateAssignVector(Op::kMul, other); }
  Vector<T, N> operator/=(const Vector<T, N>& other) { return OperateAssignVector(Op::kDiv, other); }

  Vector<T, N> operator+(T value) { Vector<T, N> res(*this); return res += value; }
  Vector<T, N> operator-(T value) { Vector<T, N> res(*this); return res -= value; }
  Vector<T, N> operator*(T value) { Vector<T, N> res(*this); return res *= value; }
  Vector<T, N> operator/(T value) { Vector<T, N> res(*this); return res /= value; }
  Vector<T, N> operator+(const Vector<T, N>& other) { Vector<T, N> res(*this); return res += other; }
  Vector<T, N> operator-(const Vector<T, N>& other) { Vector<T, N> res(*this); return res -= other; }
  Vector<T, N> operator*(const Vector<T, N>& other) { Vector<T, N> res(*this); return res *= other; }
  Vector<T, N> operator/(const Vector<T, N>& other) { Vector<T, N> res(*this); return res /= other; }

  bool operator==(const Vector<T, N>& other) const;
  friend std::ostream& operator<< <> (std::ostream& os, const Vector<T, N>& v);
 private:
  std::size_t size_;
  T* a_;

  enum Op {
    kAdd,
    kSub,
    kMul,
    kDiv
  };

  Vector<T, N> OperateAssignValue(Op op, T value);
  Vector<T, N> OperateAssignVector(Op op, const Vector<T, N>& other);

  void CheckSize(const Vector<T, N>& other) const;
};

template<typename T, std::size_t N>
Vector<T, N>::Vector(T value, std::size_t size) : size_(size) {
  a_ = new T[size];
  for (std::size_t i = 0; i < size; ++i) {
    a_[i] = value;
  }
}

template<typename T, std::size_t N>
Vector<T, N>::Vector(T* a, std::size_t size) : size_(size) {
  a_ = new T[size_];
  std::copy(a, a + size_, a_);
}

template<typename T, std::size_t N>
Vector<T, N>::Vector(Vector<T, N> const& other) : size_(other.size_) {
  a_ = new T[size_];
  for (std::size_t i = 0; i < size_; ++i) {
    a_[i] = other.a_[i];
  }
}

template<typename T, std::size_t N>
Vector<T, N>& Vector<T, N>::operator=(Vector<T, N> const& other) {
  if (this == &other) {
    return *this;
  }
  size_ = other.size_;
  a_ = new T[size_];
  for (std::size_t i = 0; i < size_; ++i) {
    a_[i] = other.a_[i];
  }
  return *this;
}

template <typename T, std::size_t N>
double Vector<T, N>::Magnitude() const{
  double res = 0;
  for (std::size_t i = 0; i < size_; ++i) {
    res += static_cast<double>(a_[i] * a_[i]);
  }
  return sqrt(res);
}

template<typename T, std::size_t N>
void Vector<T, N>::Normalize() {
  double magnitude = Magnitude();
  for (std::size_t i = 0; i < size_; ++i) {
    a_[i] /= magnitude;
  }
}

template<typename T, std::size_t N>
template<std::size_t M>
Matrix<T, N, M> Vector<T, N>::CrossProduct(const Vector<T, M> &other) const {
  Matrix<T, N, M> res(size_, other.size());
  for (std::size_t i = 0; i < size_; ++i) {
    for (std::size_t j = 0; j < other.size(); ++j) {
      res(i, j) = a_[i] * other[j];
    }
  }
  return res;
}

template<typename T, std::size_t N>
T Vector<T, N>::DotProduct(const Vector<T, N> &other) const {
  CheckSize(other);
  T res = 0;
  for (std::size_t i = 0; i < size_; ++i) {
    res += a_[i] * other.a_[i];
  }
  return res;
}

template<typename T, std::size_t N>
Vector<T> Vector<T, N>::Slice(std::size_t start, std::size_t end) const {
  Vector<T> res(end - start);
  for (std::size_t i = 0; i < (end - start); ++i) {
    res[i] = a_[i + start];
  }
  return res;
}

template<typename T, std::size_t N>
T Vector<T, N>::Sum() const {
  T sum = T();
  for (std::size_t i = 0; i < size_; ++i) {
    sum += a_[i];
  }
  return sum;
}

template<typename T, std::size_t N>
void Vector<T, N>::Fill(T value) {
  for (std::size_t i = 0; i < size_; i++) {
    a_[i] = value;
  }
}

template<typename T, std::size_t N>
Vector<T, N> Vector<T, N>::OperateAssignValue(Op op, T value) {
  for (std::size_t i = 0; i < size_; ++i) {
    switch (op) {
      case Op::kAdd:
        a_[i] += value;
        break;
      case Op::kSub:
        a_[i] -= value;
        break;
      case Op::kMul:
        a_[i] *= value;
        break;
      case Op::kDiv:
        a_[i] /= value;
        break;
    }
  }
  return *this;
}

template<typename T, std::size_t N>
Vector<T, N> Vector<T, N>::OperateAssignVector(Op op, const Vector<T, N>& other) {
  CheckSize(other);
  for (std::size_t i = 0; i < size_; ++i) {
    switch (op) {
      case Op::kAdd:
        a_[i] += other.a_[i];
        break;
      case Op::kSub:
        a_[i] -= other.a_[i];
        break;
      case Op::kMul:
        a_[i] *= other.a_[i];
        break;
      case Op::kDiv:
        a_[i] /= other.a_[i];
        break;
    }
  }
  return *this;
}

template<typename T, std::size_t N>
void Vector<T, N>::CheckSize(const Vector<T, N> &other) const {
  if (size_ != other.size_) {
    throw std::invalid_argument("Vectors must have the same size");
  }
}

template<typename T, std::size_t N>
bool Vector<T, N>::operator==(const Vector<T, N>& other) const {
  if (size_ != other.size_) {
    return false;
  }
  for (std::size_t i = 0; i < size_; ++i) {
    if (a_[i] != other.a_[i]) {
      return false;
    }
  }
  return true;
}

template<typename T, std::size_t N>
std::ostream &operator<<(std::ostream &os, const Vector<T, N> &v) {
  for (std::size_t i = 0; i < v.size_; i++) {
    os << v[i] << '\t';
  }
  return os << std::endl;
}

// TODO: Vector * Matrix
