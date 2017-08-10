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

// Basic mathematical typedef and functions.

#ifndef MSS_MATH_H
#define MSS_MATH_H

#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <typeindex>
#include <vector>

namespace mss {

typedef std::complex<double> dcomp;

const dcomp ii(0.0, 1.0);
const double ee(2.71828182845904523536);
const double pi(3.14159265358979323846);
const double pi2(pi * 2);
const double pi_2(pi / 2);
const double epsilon(1e-14);
const auto setMaxPrecision =
    std::setprecision(std::numeric_limits<double>::digits10 + 1);

inline dcomp operator*(const dcomp& lhs, int rhs) {
  return lhs * double(rhs);
}

inline dcomp operator*(int lhs, const dcomp& rhs) {
  return double(lhs) * rhs;
}

inline dcomp Jn(int n, const double& x) {
  return jn(n, x);
}

inline dcomp Hn(int n, const double& x) {
  return jn(n, x) + ii * yn(n, x);
}

inline bool AngEqu(const double& a, const double& b) {
  double t = (a - b) / pi / 2;
  return std::abs(t - (long long)t) * 2 * pi < epsilon;
}

template <int p, typename T>
inline double Lp(std::initializer_list<T> l) {
  assert(p > 0);
  double sum = 0;
  for (T i : l) sum += pow(std::abs(i), p);
  return pow(sum, 1.0 / p);
}

template <typename T>
inline double RelativeDiff(const T& a, const T& b) {
  return std::abs(a - b) / std::max(std::abs(a), std::abs(b));
}

// Check if the two values are approximately equal with relative error.
template <typename T>
inline bool ApproxRv(const T& a, const T& b, const double& re = epsilon) {
  return RelativeDiff(a, b) < re;
}

// Check if the elements of two vectors are approximately equal with relative
// error.
template <typename T>
inline bool ApproxVectRv(const Eigen::Matrix<T, Eigen::Dynamic, 1>& a,
                         const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
                         const double& re = epsilon, int k = 0) {
  assert(a.size() == b.size());
  for (long i = k; i < a.size() - k; i++)
    if (!ApproxRv(a(i), b(i), re)) return false;
  return true;
}

inline Eigen::VectorXcd IntVect(const size_t& p, int n) {
  Eigen::VectorXcd g(p), h(p);
  double d    = pi2 / p;
  long long P = p;
  if (n)
    for (long long i = 0; i < P; i++) {
      g(i) = exp(-(i + 1) * n * d * ii) / (n * n * d) *
             (exp(n * d * ii) - n * d * ii * exp(n * d * ii) - 1.0);
      h(i) = exp(-(i + 1) * n * d * ii) / (n * n * d) *
             (-exp(n * d * ii) + n * d * ii + 1.0);
    }
  else {
    g.setConstant(d / 2);
    h.setConstant(d / 2);
  }

  Eigen::VectorXcd rst(p);
  rst(0) = g(0) + h(p - 1);
  for (size_t i = 1; i < p; i++) rst(i)= h(i - 1) + g(i);
  rst /= pi2;

  return rst;
}

inline Eigen::MatrixXcd IntMat(int m, const size_t& p, int n) {
  Eigen::VectorXcd v = IntVect(p, n);
  Eigen::MatrixXcd rst(m, m * p);
  rst.setZero();
  for (int i = 0; i < m; i++)
    for (size_t j = 0; j < p; j++) rst(i, i + m * j) = v(j);
  return rst;
}

}  // namespace mss

#endif
