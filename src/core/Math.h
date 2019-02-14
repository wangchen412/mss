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

#include <complex_bessel.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
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
enum BoundaryShape { RECTANGULAR, HEXAGONAL, CIRCULAR };

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

inline dcomp Jn(int n, const dcomp& z) {
  return sp_bessel::besselJ(n, z);
}

inline dcomp Yn(int n, double x) {
  return yn(n, x);
}

inline dcomp Yn(int n, const dcomp& z) {
  return sp_bessel::besselY(n, z);
}

inline dcomp Hn(int n, double x) {
  return jn(n, x) + ii * yn(n, x);
}

inline dcomp Hn(int n, const dcomp& z) {
  return sp_bessel::hankelH1(n, z);
}

inline dcomp H2n(int n, double x) {
  return jn(n, x) - ii * yn(n, x);
}

inline dcomp H2n(int n, const dcomp& z) {
  return sp_bessel::hankelH2(n, z);
}

// The functor f returns the value of f(x) / f'(x).
template <typename Func>
double Newton(const Func& f, double x0, double e = 1e-16) {
  double x = x0;
  for (double dx = f(x); std::abs(dx) > e; dx = f(x)) x -= dx;
  return x -= f(x);
}

// The function/functor f takes Eigen::VectorXd and returns a double.
template <typename Func, int n>
Eigen::Matrix<double, n, 1> Gradient(const Func& f,
                                     const Eigen::Matrix<double, n, 1>& x,
                                     double d = 1e-2) {
  double fx = f(x);
  std::cout << "Residue: " << fx << std::endl;
  Eigen::VectorXd dy(n);
  Eigen::MatrixXd dx = Eigen::MatrixXd::Identity(n, n) * d;
  for (int i = 0; i < n; i++) dy(i) = (f(x + dx.col(i)) - fx);
  return dy / d;
}

template <typename Func, int n>
Eigen::VectorXd GradientDescent(const Func& f,
                                const Eigen::Matrix<double, n, 1>& x0,
                                double d = 5e-3, double e = 1e-4,
                                size_t max_iter = 1e3) {
  auto x(x0), g(Gradient(f, x0));
  for (size_t i = 0; i < max_iter && g.norm() > e; i++, g = Gradient(f, x)) {
    x -= g * d;
    std::cout << "Gradient: " << g.transpose() << std::endl;
    std::cout << "Properties: " << x.transpose() << std::endl;
  }
  if (g.norm() > e) std::cout << "Not converged." << std::endl;
  return x;
}

// Return a tuple of Pn(x) and P'n(x).
inline std::pair<double, double> Legendre(int N, double x) {
  // assert(N > 2);
  double pn = x, pn_1 = x, pn_2 = 1;
  for (int i = 2; i <= N; i++) {
    pn = ((2 * i - 1) * x * pn_1 - (i - 1) * pn_2) / i;
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
      double x = Newton(dx, cos(pi * (i + 0.75) / (N + 0.5)));
      double d = Legendre(N, x).second;
      root_[i] = x;
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
  // If one of a and b is zero, then return the norm of the other one as the
  // difference.

  double aa(std::abs(a)), ab(std::abs(b));
  double ac(std::max(aa, ab));
  if (aa * ab == 0) return ac;
  return std::abs(a - b) / ac;
}

// Check if the two values are approximately equal with relative error.
template <typename T>
bool ApproxRv(const T& a, const T& b, double re = epsilon, bool v = false,
              std::ostream& os = std::cout) {
  if (v)
    if (RelativeDiff(a, b) > re)
      os << a << "\t" << b << "\t" << RelativeDiff(a, b) << std::endl;
  return RelativeDiff(a, b) < re;
}

// Check if the elements of two vectors are approximately equal with relative
// error.
template <typename T>
bool ApproxVectRv(const Eigen::Matrix<T, Eigen::Dynamic, 1>& a,
                  const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
                  double re = epsilon, int k = 0, bool v = false,
                  std::ostream& os = std::cout) {
  assert(a.size() == b.size());

  bool rst = true;
  for (long i = k; i < a.size() - k; i++)
    if (!ApproxRv(a(i), b(i), re, v, os)) {
      if (v)
        rst = false;
      else
        return false;
    }
  return rst;
}

template <typename T>
bool ApproxMatRv(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& a,
                 const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& b,
                 double re = epsilon, int k = 0, bool v = false,
                 std::ostream& os = std::cout) {
  assert(a.size() == b.size());

  bool rst = true;
  for (long i = 0; i < a.cols(); i++)
    if (!ApproxVectRv<T>(a.col(i), b.col(i), re, k, v, os)) {
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
  T rst = v.back();
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

MatrixXcd PseudoInverse(const MatrixXcd& m) {
  auto svd = m.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
  auto& sv = svd.singularValues();
  double tol = std::numeric_limits<double>::epsilon() * m.rows() *
               sv.array().abs().maxCoeff();
  MatrixXcd svi(m.cols(), m.cols());
  svi.setZero();
  for (long i = 0; i < m.cols(); i++)
    svi(i, i) = std::abs(sv(i)) > tol ? 1 / sv(i) : 0;
  return svd.matrixV() * svi * svd.matrixU().adjoint();
}

double ConditionNum(const MatrixXcd& m) {
  auto svd = m.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
  auto& sv = svd.singularValues();
  return sv(0) / sv(sv.size() - 1);
}

template <typename T>
typename T::Scalar GeometricMean(const Eigen::ArrayBase<T>& a) {
  typename T::Scalar t = 0;
  for (long i = 0; i < a.size(); i++) t += log(a(i));
  return exp(t / double(a.size()));
}

double ApproxIdentity(const MatrixXcd& m) {
  MatrixXcd I(m.rows(), m.cols());
  I.setIdentity();
  return (m - I).norm();
}

std::vector<long> FindUnitEigenvalue(const VectorXcd& v,
                                     double tol = epsilon) {
  std::vector<long> rst;
  for (long i = 0; i < v.size(); i++) {
    double norm = sqrt(v(i).real() * v(i).real() + v(i).imag() * v(i).imag());
    if (norm > 1 - tol && norm < 1 + tol) rst.push_back(i);
  }
  return rst;
}

}  // namespace mss

#include "Integrators.h"

#endif
