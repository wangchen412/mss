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

using Eigen::Matrix2cd;
using Eigen::Matrix4cd;
using Eigen::MatrixXcd;
using Eigen::Vector2cd;
using Eigen::Vector4cd;
using Eigen::VectorXcd;

template <typename T>
using MatrixNcd = Eigen::Matrix<dcomp, T::NumBv, T::NumBv>;

template <typename T>
using VectorNcd = Eigen::Matrix<dcomp, T::NumBv, 1>;

enum SolveMethod { COLLOCATION, DFT };
enum BoundaryShape { RECTANGULAR, HEXAGONAL };

inline dcomp operator+(const dcomp& lhs, int rhs) {
  return lhs + double(rhs);
}

inline dcomp operator+(int lhs, const dcomp& rhs) {
  return double(lhs) + rhs;
}

inline dcomp operator-(const dcomp& lhs, int rhs) {
  return lhs - double(rhs);
}

inline dcomp operator-(int lhs, const dcomp& rhs) {
  return double(lhs) - rhs;
}

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
double Newton(const Func& f, double x0, double e = 1e-16) {
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
double Lp(std::initializer_list<T> l) {
  assert(p > 0);
  double sum = 0;
  for (T i : l) sum += pow(std::abs(i), p);
  return pow(sum, 1.0 / p);
}

template <typename T>
double RelativeDiff(const T& a, const T& b) {
  return std::abs(a - b) / std::max(std::abs(a), std::abs(b));
}

// Check if the two values are approximately equal with relative error.
template <typename T>
bool ApproxRv(const T& a, const T& b, double re = epsilon, bool v = false) {
  if (v)
    std::cout << a << "\t" << b << "\t" << RelativeDiff(a, b) << std::endl;
  return RelativeDiff(a, b) < re;
}

// Check if the elements of two vectors are approximately equal with relative
// error.
template <typename T>
bool ApproxVectRv(const Eigen::Matrix<T, Eigen::Dynamic, 1>& a,
                  const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
                  double re = epsilon, int k = 0, bool v = false) {
  assert(a.size() == b.size());

  bool rst = true;
  for (long i = k; i < a.size() - k; i++)
    if (!ApproxRv(a(i), b(i), re, v)) {
      if (v)
        rst = false;
      else
        return false;
    }
  return rst;
}

template <typename T>
std::vector<T> DeterminantFactor(const Eigen::Matrix<T, -1, -1>& m) {
  auto lu = m.partialPivLu();

  Eigen::Matrix<T, -1, 1> diag = lu.matrixLU().diagonal();

  T p = 0;
  for (int i = 0; i < diag.size(); i++) p += std::log10(diag(i));
  p /= diag.size();

  std::vector<T> v;
  for (int i = 0; i < diag.size(); i++) v.push_back(diag(i) / pow(10, p));
  std::sort(v.begin(), v.end(),
            [](const T& a, const T& b) { return std::abs(a) < std::abs(b); });

  v.push_back(p);
  v.push_back(lu.permutationP().determinant());

  return v;
}

template <typename T>
std::tuple<T, T> Determinant(const Eigen::Matrix<T, -1, -1>& m) {
  auto v = DeterminantFactor(m);
  T rst  = v.back();
  v.pop_back();

  T ave_p = v.back();
  v.pop_back();

  size_t i = 0, j = v.size() - 1;
  while (i < j) rst *= v[i++] * v[j--];
  if (i == j) rst *= v[i];

  return std::make_tuple(rst, ave_p);
}

double DeterExpon(const MatrixXcd& m) {
  VectorXcd diag = m.partialPivLu().matrixLU().diagonal();

  dcomp p = 0;
  for (long i = 0; i < diag.size(); i++) p += std::log10(diag(i));
  p /= diag.size();

  return p.real();
}

}  // namespace mss

#include "Integrators.h"

#endif
