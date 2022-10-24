//
// Created by Mipku on 20.10.2022.
//

#ifndef DZ2_VECTOR_H
#define DZ2_VECTOR_H

#include <array>
#include <iostream>
#include <cassert>
#include <cmath>

template <typename T, std::size_t Rows, std::size_t Cols>
class Matrix;

template <typename T, std::size_t N = 0>
class Vector {
 public:
  Vector();
  Vector(Vector const&);
  explicit Vector(T* a) : a_(a), size_(N) {}
  explicit Vector(std::size_t size);
  Vector(T* a, std::size_t size) : a_(a), size_(size) {}

  std::size_t size() const { return size_; }

  Matrix<T, 0, 0> CrossProduct(const Vector<T, N>& other) const;
  T DotProduct(const Vector<T, N>& other) const;
  Vector<T> Slice(std::size_t start, std::size_t end) const;
  Vector<T> Slice(std::size_t start) { return Slice(start, size()); };
  std::array<T, N> UnderlyingArray() { return a_; }
  double Magnitude() const;
  void Normalize();
  T Sum() const;
  void Print() const;

  T Get(std::size_t i) const { return a_[i]; }
  T operator[](std::size_t i) const { return a_[i]; };
  T& operator[](std::size_t i) { return a_[i]; }

  Vector<T, N> operator+(T value);
  Vector<T, N> operator-(T value);
  Vector<T, N> operator*(T value);
  Vector<T, N> operator/(T value);
  Vector<T, N> operator+(const Vector<T, N>& other);
  Vector<T, N> operator-(const Vector<T, N>& other);
  Vector<T, N> operator*(const Vector<T, N>& other);
  Vector<T, N> operator/(const Vector<T, N>& other);

  Vector<T, N> operator+=(T value);
  Vector<T, N> operator-=(T value);
  Vector<T, N> operator*=(T value);
  Vector<T, N> operator/=(T value);
  Vector<T, N> operator+=(const Vector<T, N>& other);
  Vector<T, N> operator-=(const Vector<T, N>& other);
  Vector<T, N> operator*=(const Vector<T, N>& other);
  Vector<T, N> operator/=(const Vector<T, N>& other);

 private:
  std::size_t size_;
  T* a_;
};

template<typename T, std::size_t N>
Vector<T, N>::Vector() : size_(N) {
  a_ = new T[size_];
  for (std::size_t i = 0; i < size_; i++)
    a_[i] = T();
}

template<typename T, std::size_t N>
Vector<T, N>::Vector(Vector<T, N> const& other) : size_(other.size_) {
  a_ = new T[size_];
  for (std::size_t i = 0; i < size_; i++)
    a_[i] = other.a_[i];
}

template<typename T, std::size_t N>
Vector<T, N>::Vector(std::size_t size) : size_(size) {
  a_ = new T[size];
  for (std::size_t i = 0; i < size; i++) {
    a_[i] = T();
  }
}

template <typename T, std::size_t N>
double Vector<T, N>::Magnitude() const{
  double res = 0;
  for (std::size_t i = 0; i < size_; i++) {
    res += a_[i] * a_[i];
  }
  return sqrt(res);
}

template<typename T, std::size_t N>
void Vector<T, N>::Normalize() {
  double magnitude = Magnitude();
  for (std::size_t i = 0; i < size_; i++) {
    a_[i] /= magnitude;
  }
}

template<typename T, std::size_t N>
Matrix<T, 0, 0> Vector<T, N>::CrossProduct(const Vector<T, N> &other) const {
  Matrix<T, 0, 0> res(size_, other.size_);
  for (std::size_t i = 0; i < size_; i++) {
    for (std::size_t j = 0; j < other.size_; j++) {
      T product = a_[i] * other.a_[j];
      res.Set(i, j, product);
    }
  }
  return res;
}

template<typename T, std::size_t N>
T Vector<T, N>::DotProduct(const Vector<T, N> &other) const {
  assert(size_ == other.size_);
  T res = 0;
  for (std::size_t i = 0; i < size_; i++) {
    res += a_[i] * other.a_[i];
  }
  return res;
}

template<typename T, std::size_t N>
Vector<T> Vector<T, N>::Slice(std::size_t start, std::size_t end) const {
  Vector<T> res(end - start);
  for (std::size_t i = 0; i < (end - start); i++) {
    res[i] = a_[i + start];
  }
  return res;
}

template<typename T, std::size_t N>
T Vector<T, N>::Sum() const {
  T sum = T();
  for (std::size_t i = 0; i < size_; i++) {
    sum += a_[i];
  }
  return sum;
}

template<typename T, std::size_t N>
void Vector<T, N>::Print() const {
  for (std::size_t i = 0; i < size_; i++) {
    std::cout << a_[i] << ' ';
  }
  std::cout << std::endl;
}

template<typename T, std::size_t N>
Vector<T, N> Vector<T, N>::operator+(const T value) {
  Vector<T, N> res(size_);
  for (std::size_t i = 0; i < size_; i++) {
    res[i] = a_[i] + value;
  }
  return res;
}

template<typename T, std::size_t N>
Vector<T, N> Vector<T, N>::operator-(const T value) {
  Vector<T, N> res(size_);
  for (std::size_t i = 0; i < size_; i++) {
    res[i] = a_[i] - value;
  }
  return res;
}

template<typename T, std::size_t N>
Vector<T, N> Vector<T, N>::operator*(T value) {
  Vector<T, N> res(size_);
  for (std::size_t i = 0; i < size_; i++) {
    res[i] = a_[i] * value;
  }
  return res;
}

template<typename T, std::size_t N>
Vector<T, N> Vector<T, N>::operator/(T value) {
  Vector<T, N> res(size_);
  for (std::size_t i = 0; i < size_; i++) {
    res[i] = a_[i] / value;
  }
  return res;
}

template<typename T, std::size_t N>
Vector<T, N> Vector<T, N>::operator+(const Vector<T, N> &other) {
  Vector<T, N> res(size_);
  for (std::size_t i = 0; i < size_; i++) {
    res[i] = a_[i] + other.a_[i];
  }
  return res;
}

template<typename T, std::size_t N>
Vector<T, N> Vector<T, N>::operator-(const Vector<T, N> &other) {
  Vector<T, N> res(size_);
  for (std::size_t i = 0; i < size_; i++) {
    res[i] = a_[i] - other.a_[i];
  }
  return res;
}

template<typename T, std::size_t N>
Vector<T, N> Vector<T, N>::operator*(const Vector<T, N> &other) {
  Vector<T, N> res(size_);
  for (std::size_t i = 0; i < size_; i++) {
    res[i] = a_[i] * other.a_[i];
  }
  return res;
}

template<typename T, std::size_t N>
Vector<T, N> Vector<T, N>::operator/(const Vector<T, N> &other) {
  Vector<T, N> res(size_);
  for (std::size_t i = 0; i < size_; i++) {
    res[i] = a_[i] / other.a_[i];
  }
  return res;
}

template<typename T, std::size_t N>
Vector<T, N> Vector<T, N>::operator+=(T value) {
  for (std::size_t i = 0; i < size_; i++) {
    a_[i] += value;
  }
  return *this;
}

template<typename T, std::size_t N>
Vector<T, N> Vector<T, N>::operator-=(T value) {
  for (std::size_t i = 0; i < size_; i++) {
    a_[i] -= value;
  }
  return *this;
}

template<typename T, std::size_t N>
Vector<T, N> Vector<T, N>::operator*=(T value) {
  for (std::size_t i = 0; i < size_; i++) {
    a_[i] *= value;
  }
  return *this;
}

template<typename T, std::size_t N>
Vector<T, N> Vector<T, N>::operator/=(T value) {
  for (std::size_t i = 0; i < size_; i++) {
    a_[i] /= value;
  }
  return *this;
}

template<typename T, std::size_t N>
Vector<T, N> Vector<T, N>::operator+=(const Vector<T, N> &other) {
  for (std::size_t i = 0; i < size_; i++) {
    a_[i] += other.a_[i];
  }
  return *this;
}

template<typename T, std::size_t N>
Vector<T, N> Vector<T, N>::operator-=(const Vector<T, N> &other) {
  for (std::size_t i = 0; i < size_; i++) {
    a_[i] -= other.a_[i];
  }
  return *this;
}

template<typename T, std::size_t N>
Vector<T, N> Vector<T, N>::operator*=(const Vector<T, N> &other) {
  for (std::size_t i = 0; i < size_; i++) {
    a_[i] *= other.a_[i];
  }
  return *this;
}

template<typename T, std::size_t N>
Vector<T, N> Vector<T, N>::operator/=(const Vector<T, N> &other) {
  for (std::size_t i = 0; i < size_; i++) {
    a_[i] /= other.a_[i];
  }
  return *this;
}

#endif //DZ2_VECTOR_H
