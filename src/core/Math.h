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
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <tuple>
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

using Eigen::MatrixXcd;
using Eigen::VectorXcd;

template <typename T>
using MatrixNcd = Eigen::Matrix<dcomp, T::NumBv, T::NumBv>;

enum SolveMethod { COLLOCATION, DFT };
enum Tessellation { RECTANGULAR, HEXAGONAL };

inline dcomp operator*(const dcomp& lhs, int rhs) {
  return lhs * double(rhs);
}

inline dcomp operator*(int lhs, const dcomp& rhs) {
  return double(lhs) * rhs;
}

inline dcomp operator/(const dcomp& lhs, int rhs) {
  return lhs / double(rhs);
}

inline dcomp operator/(int lhs, const dcomp& rhs) {
  return double(lhs) / rhs;
}

inline dcomp Jn(int n, double x) {
  return jn(n, x);
}

inline dcomp Hn(int n, double x) {
  return jn(n, x) + ii * yn(n, x);
}

// The functor f returns the value of f(x) / f'(x).
template <typename Func>
inline double Newton(const Func& f, double x0, double e = 1e-16) {
  double x = x0;
  for (double dx = f(x); std::abs(dx) > e; dx = f(x)) x -= dx;
  return x -= f(x);
}

// Return a tuple of Pn(x) and P'n(x).
inline std::pair<double, double> Legendre(int N, double x) {
  // assert(N > 2);
  double pn = x, pn_1 = x, pn_2 = 1;
  for (int i = 2; i <= N; i++) {
    pn   = ((2 * i - 1) * x * pn_1 - (i - 1) * pn_2) / i;
    pn_2 = pn_1;
    pn_1 = pn;
  }
  return std::make_pair(pn, N * (x * pn - pn_2) / (x * x - 1));
}

template <int N>
class LegendreRoot {
 public:
  LegendreRoot() : root_(N), weight_(N) {
    for (int i = 0; i < N; i++) {
      double x   = Newton(dx, cos(pi * (i + 0.75) / (N + 0.5)));
      double d   = Legendre(N, x).second;
      root_[i]   = x;
      weight_[i] = 2 / ((1 - x * x) * d * d);
    }
  }

  double root(int n) const { return root_[n]; }
  double weight(int n) const { return weight_[n]; }

 private:
  std::vector<double> root_;
  std::vector<double> weight_;

  static double dx(double x) {
    auto p = Legendre(N, x);
    return p.first / p.second;
  }
};

inline bool AngEqu(double a, double b) {
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
inline bool ApproxRv(const T& a, const T& b, double re = epsilon) {
  return RelativeDiff(a, b) < re;
}

// Check if the elements of two vectors are approximately equal with relative
// error.
template <typename T>
inline bool ApproxVectRv(const Eigen::Matrix<T, Eigen::Dynamic, 1>& a,
                         const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
                         double re = epsilon, int k = 0) {
  assert(a.size() == b.size());
  for (long i = k; i < a.size() - k; i++)
    if (!ApproxRv(a(i), b(i), re)) return false;
  return true;
}

}  // namespace mss

#include "Integrators.h"

#endif
