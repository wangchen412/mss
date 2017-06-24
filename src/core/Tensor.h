// ----------------------------------------------------------------------
//
// Copyright © 2017 mss authors.
//
// mss is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// mss is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ----------------------------------------------------------------------

// Tensor class templates.

#ifndef MSS_TENSOR_H
#define MSS_TENSOR_H

#include <fstream>
#include <iostream>
#include "Math.h"

namespace mss {

template <typename T>
struct Vector {
  explicit Vector(const T& x = 0, const T& y = 0) : x(x), y(y) {}
  Vector(const Vector& other) : x(other.x), y(other.y) {}
  ~Vector() {}

  Vector& operator+=(const Vector& other);
  Vector& operator-=(const Vector& other);
  Vector operator+(const Vector& other) const;
  Vector operator-(const Vector& other) const;
  bool operator==(const Vector& other) const;

  Vector& RotateInPlace(const double& angle);
  Vector Rotate(const double& angle) const;

  double Length() const;
  double Angle() const;
  double Angle(const double& L) const;
  Vector Polar() const;

  T x, y;
};

template <typename T>
struct Tensor {
  explicit Tensor(const T& xx = 0, const T& yy = 0, const T& xy = 0)
      : xx(xx), yy(yy), xy(xy) {}
  Tensor(const Tensor& other) : xx(other.xx), yy(other.yy), xy(other.xy) {}
  ~Tensor() {}

  Tensor& operator+=(const Tensor& other);
  Tensor& operator-=(const Tensor& other);
  Tensor operator+(const Tensor& other) const;
  Tensor operator-(const Tensor& other) const;
  bool operator==(const Tensor& other) const;

  Tensor& RotateInPlace(const double& angle);
  Tensor Rotate(const double& angle) const;

  T xx, yy, xy;
};

// Typedefs for antiplane (AP) and in-plane (IP) problems.
typedef dcomp DispAP;
typedef Vector<dcomp> StressAP;
typedef Vector<dcomp> DispIP;
typedef Tensor<dcomp> StressIP;
typedef Vector<double> PosiVect;

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T>
inline Vector<T>& Vector<T>::operator+=(const Vector<T>& other) {
  x += other.x;
  y += other.y;
  return *this;
}
template <typename T>
inline Vector<T> Vector<T>::operator+(const Vector<T>& other) const {
  return Vector<T>(x + other.x, y + other.y);
}
template <typename T>
inline Vector<T>& Vector<T>::operator-=(const Vector<T>& other) {
  x -= other.x;
  y -= other.y;
  return *this;
}
template <typename T>
inline Vector<T> Vector<T>::operator-(const Vector<T>& other) const {
  return Vector<T>(x - other.x, y - other.y);
}
template <typename T>
inline bool Vector<T>::operator==(const Vector<T>& other) const {
  if(x == other.x && y == other.y) return true;
  return Lp<2>({ x - other.x, y - other.y }) /
    std::max(Lp<2>({ x, y }), Lp<2>({ other.x, other.y })) < epsilon;
}
template <typename T>
inline std::ostream& operator<<(std::ostream& os, const Vector<T>& v) {
  return os << v.x << "\t" << v.y;
}

template <typename T>
inline std::istream& operator>>(std::istream& is, Vector<T>& v) {
  return is >> v.x >> v.y;
}

template <typename T>
inline Vector<T>& Vector<T>::RotateInPlace(const double& a) {
  // Change its components into the ones in a rotated CS.

  double c = cos(a), s = sin(a);
  T t1 = x, &t2 = y;
  x = t1 * c + t2 * s;
  y = -t1 * s + t2 * c;
  return *this;
}
template <typename T>
inline Vector<T> Vector<T>::Rotate(const double& a) const {
  // Return a new vector with components in a rotated C.S

  Vector<T> rst(*this);
  rst.RotateInPlace(a);
  return rst;
}
template <typename T>
inline double Vector<T>::Length() const {
  // Return the length of the real vector.

  assert(typeid(T) == typeid(double));
  return Lp<2>({x, y});
}
template <typename T>
inline double Vector<T>::Angle() const {
  // Return the angle the real vector without knowing its length.

  return Angle(Length());
}
template <typename T>
inline double Vector<T>::Angle(const double& L) const {
  // Return the angle of the real vector with knowing its length.

  if (L < epsilon) return 0;
  assert(typeid(T) == typeid(double));
  return x > 0 ? asin(y / L) : pi - asin(y / L);
}
template <typename T>
inline Vector<T> Vector<T>::Polar() const {
  // Return the components of the position vector as if it is in a polar
  // coordinate system, of which the polar axis is aligned along the x axis.
  // The first value is the radius, and the second one is the polar angle.

  assert(typeid(T) == typeid(double));
  double r = Length();
  return Vector<T>(r, Angle(r));
}

// ----------------------------------------------------------------------

template <typename T>
inline Tensor<T>& Tensor<T>::operator+=(const Tensor<T>& other) {
  xx += other.xx;
  yy += other.yy;
  xy += other.xy;
  return *this;
}
template <typename T>
inline Tensor<T> Tensor<T>::operator+(const Tensor<T>& other) const {
  return Tensor<T>(xx + other.xx, yy + other.yy, xy + other.xy);
}
template <typename T>
inline Tensor<T>& Tensor<T>::operator-=(const Tensor<T>& other) {
  xx -= other.xx;
  yy -= other.yy;
  xy -= other.xy;
  return *this;
}
template <typename T>
inline Tensor<T> Tensor<T>::operator-(const Tensor<T>& other) const {
  return Tensor<T>(xx - other.xx, yy - other.yy, xy - other.xy);
}
template <typename T>
bool Tensor<T>::operator==(const Tensor<T>& other) const {
  if(xx == other.xx && yy == other.yy && xy == other.xy) return true;
  return Lp<2>({xx - other.xx, yy - other.yy, xy - other.xy}) /
    std::max(Lp<2>({ xx, yy, xy }), Lp<2>({ other.xx, other.yy, other.xy })) < epsilon;
}
template <typename T>
inline std::ostream& operator<<(std::ostream& os, const Tensor<T>& t) {
  return os << t.xx << "\t" << t.yy << "\t" << t.xy;
}
template <typename T>
inline std::istream& operator>>(std::istream& is, Tensor<T>& t) {
  return is >> t.xx >> t.yy >> t.xy;
}
template <typename T>
inline Tensor<T>& Tensor<T>::RotateInPlace(const double& a) {
  // Change its components into the ones in a rotated CS.

  double c2 = cos(a * 2), s2 = sin(a * 2);
  double cc = pow(cos(a), 2), ss = pow(sin(a), 2);
  T t1 = xx, t2 = yy, &t3 = xy;
  xx = t1 * cc + t2 * ss + t3 * s2;
  yy = t1 * ss + t2 * cc - t3 * s2;
  xy = (t2 - t1) * 0.5 * s2 + t3 * c2;
  return *this;
}
template <typename T>
inline Tensor<T> Tensor<T>::Rotate(const double& a) const {
  // Return a new tensor with components in a rotated CS.

  Tensor<T> rst(*this);
  rst.RotateInPlace(a);
  return rst;
}

}  // namespace mss

#endif
