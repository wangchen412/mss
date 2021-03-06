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
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <random>
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
enum BoundaryShape { INCL_RECT, RECTANGULAR, HEXAGONAL, CIRCULAR, INPUT };

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

// The functor f takes Eigen::VectorXd and returns a double.
template <typename Func>
Eigen::VectorXd Gradient(const Func& f, const Eigen::VectorXd& x,
                         std::ostream* os = nullptr, double d = 1e-4) {
  double fx = f(x);
  if (os != nullptr) *os << fx << "\t";
  int n = x.rows();
  Eigen::VectorXd dy(n);
  Eigen::MatrixXd dx = Eigen::MatrixXd::Identity(n, n) * d;
  for (int i = 0; i < n; i++) dy(i) = f(x + dx.col(i)) - fx;
  return dy / d;
}

template <typename Func>
Eigen::VectorXd GradientDescent(
    const Func& f, std::ostream* os = nullptr,
    const Eigen::Vector4d& x0 = Eigen::Vector4d::Ones(), double d = 1e-3,
    double e = 1e-4, size_t max_iter = 1e4) {
  Eigen::VectorXd x(x0), g(Gradient(f, x0));
  for (size_t i = 0; i < max_iter && g.norm() > e; i++) {
    x -= g * d;
    g = Gradient(f, x, os);
    if (os != nullptr)
      *os << x.transpose() << "\t" << g.transpose() << std::endl;
  }
  // if (g.norm() > e) std::cout << "Not converged." << std::endl;
  return x;
}

template <typename T>
Eigen::VectorXi sort_index(Eigen::Matrix<T, -1, 1>& v) {
  Eigen::VectorXi x(v.size());
  x.setLinSpaced(0, x.size());
  std::sort(x.data(), x.data() + x.size(),
            [&v](size_t i, size_t j) { return v(i) < v(j); });
  std::sort(v.data(), v.data() + v.size());
  return x;
}
Eigen::MatrixXd Permutation(const Eigen::VectorXi& index) {
  Eigen::MatrixXd p = Eigen::MatrixXd::Zero(index.size(), index.size());
  for (long i = 0; i < index.size(); i++) p(index(i), i) = 1;
  return p;
}
template <typename Func>
Eigen::VectorXd NelderMead(const Func& f, const Eigen::VectorXd& x0,
                           std::ostream* os = nullptr, double e = 1e-4,
                           size_t max_iter = 1e3, double* fval = nullptr,
                           bool* conv = nullptr) {
  long N = x0.size();
  bool converge = false;

  Eigen::MatrixXd x = x0 * Eigen::VectorXd::Ones(N + 1).transpose();
  Eigen::VectorXd y(N + 1);

  for (long i = 1; i < N + 1; i++)
    x(i - 1, i) =
        std::abs(x(i - 1, i)) < epsilon ? 0.00025 : x(i - 1, i) * 1.05;
  for (long i = 0; i < N + 1; i++) y(i) = f(x.col(i));
  for (size_t i = 0; i < max_iter; i++) {
    x *= Permutation(sort_index(y));
    if (os != nullptr)
      *os << setMaxPrecision << y(0) << "\t" << x.col(0).transpose()
          << std::endl;
    if (x.col(N).norm() > 100 * x0.norm()) break;
    if ((x.col(N) - x.col(0)).norm() < e) {
      converge = true;
      break;
    }
    Eigen::VectorXd m = x.block(0, 0, N, N).rowwise().mean();
    Eigen::VectorXd r = 2 * m - x.col(N);
    double fr = f(r);
    if (fr < y(0)) {
      Eigen::VectorXd c = m + 2 * (m - x.col(N));
      double fc = f(c);
      x.col(N) = fc < fr ? c : r;
      y(N) = std::min(fc, fr);
    } else if (fr < y(N - 1)) {
      x.col(N) = r;
      y(N) = fr;
    } else {
      Eigen::VectorXd c = m + (fr < y(N) ? r : x.col(N) - m) / 2;
      double fc = f(c);
      if (fc < std::min(fr, y(N))) {
        x.col(N) = c;
        y(N) = fc;
      } else {
        for (long j = 1; j < N + 1; j++) {
          x.col(j) = x.col(0) + (x.col(j) - x.col(0)) / 2;
          y(j) = f(x.col(j));
        }
      }
    }
  }
  if (fval != nullptr) *fval = y(0);
  if (conv != nullptr) *conv = converge;
  return x.col(0);
}
template <typename T>
T perturb(const T& x, const T& step_size) {
  std::random_device rd;
  std::mt19937 mt(rd());

  T xx(x);
  for (long i = 0; i < x.size(); i++) {
    std::uniform_real_distribution<double> dist(-step_size(i), step_size(i));
    // std::normal_distribution<double> dist(0, step_size(i));
    xx(i) += dist(mt);
  }

  return xx;
}
template <typename Func>
Eigen::VectorXd BasinHopping(size_t num_iter, size_t num_hop, const Func& f,
                             const Eigen::VectorXd& x0,
                             std::ostream* os = nullptr, double e = 1e-4,
                             size_t max_iter = 1e3) {
  bool conv;
  double y;
  Eigen::VectorXd x = x0, s = x0;
  s.setOnes();

  x = NelderMead(f, x0, os, e, max_iter, &y, &conv);

  if (!conv)
    for (size_t i = 0; i < 400; i++) {
      std::cout << "Not converged, try: " << i << std::endl;
      Eigen::VectorXd s0 = s;
      x = NelderMead(f, perturb(x0, s0), os, e, max_iter, &y, &conv);
      if (conv)
        if ((x - x0).array().abs().maxCoeff() < 0.4) break;
    }

  if (!conv) {
    std::cout << "[mss]: Basin Hopping: First minimum not found."
              << std::endl;
    exit(1);
  }

  std::ofstream file("BasinHopping.dat");
  file << "First min: " << y << "  " << x.transpose() << std::endl;

  for (size_t i = 0; i < num_iter; i++) {
    size_t na = 0;
    size_t nb = 0;
    file << i << "  step size: " << s.transpose() << std::endl;

    for (size_t j = 0; j < num_hop; j++) {
      double yy;
      Eigen::VectorXd xx0 = perturb(x, s);
      Eigen::VectorXd xx = NelderMead(f, xx0, os, e, max_iter, &yy, &conv);
      file << "  " << j << " Start from: " << xx0.transpose();
      if (conv) {
        file << "\t local min: " << yy << "\t" << xx.transpose()
             << "\t Dist: " << (x - xx).norm() << std::endl;
        if (yy < y && (x - xx).norm() > e * 10) {
          x = xx;
          y = yy;
          s.setOnes();
          s *= 0.6;  // Newly found min, start from 0.75 (0.6 * 1.25).
          break;
        }
      } else {
        file << "\t not converged: " << yy << "\t" << xx.transpose()
             << std::endl;
      }
    }
    s *= 1.25;
  }

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
  if (a.size() != b.size()) {
    if (v)
      os << "[mss]: Vector sizes inconsistent: " << std::to_string(a.size())
         << ", " << std::to_string(b.size()) << std::endl;
    return false;
  }

  bool rst = true;
  for (long i = k; i < a.size() - k; i++)
    if (!ApproxRv(a(i), b(i), re, v, os)) {
      if (v)
        rst = false;
      else
        return false;
    }
  return rst;
}  // namespace mss

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

VectorXcd InverseIteration(const MatrixXcd& A, const dcomp& e,
                           size_t max_iter = 1e2, double tol = 1e-8) {
  VectorXcd v1 = VectorXcd::Random(A.cols());
  v1 /= v1.norm();
  VectorXcd v2(v1);

  std::ofstream file("iter.dat");

  size_t i = 0;
  while (i++ < max_iter) {
    v2 = (A - MatrixXcd::Identity(A.rows(), A.cols()) * e).lu().solve(v1);
    v2 /= v2.norm();
    file << v2.transpose() << std::endl;
    if ((v2 - v1).norm() < tol) {
      std::cout << i << std::endl;
      return v2;
    } else
      v1 = v2;
  }
  std::cout << "[mss:] Inverse Iteration didn't converge." << std::endl;
  return v2;
}

VectorXcd InverseIteration(const MatrixXcd& A, const MatrixXcd& B,
                           const dcomp& e, size_t max_iter = 1e2,
                           double tol = 1e-8) {
  // Av = eBv
  // (A - eB) v_n = B v_n_1

  // MatrixXcd inv = PseudoInverse(A - B * e) * B;

  VectorXcd v1 = VectorXcd::Random(A.cols());
  v1 /= v1.norm();
  VectorXcd v2(v1);

  size_t i = 0;
  while (i++ < max_iter) {
    v2 = (A - B * e).jacobiSvd(40).solve(B * v1);
    v2 /= v2.norm();
    dcomp theta = v1.dot(v2);
    if ((v2 - theta * v1).norm() < tol * std::abs(theta)) {
      std::cout << i << std::endl;
      return v2;
    } else
      v1 = v2;
  }
  std::cout << "[mss:] Inverse Iteration didn't converge." << std::endl;
  return v2;
}

MatrixXcd PhaseShift(const dcomp& ka, double angle, long m) {
  MatrixXcd rst(m, m);
  rst.setZero();
  double c = cos(angle), s = sin(angle);
  for (long i = 0; i < m / 2; i++) {
    rst(i, i) = exp(ii * ka * c);
    rst(i + m / 2, i + m / 2) = exp(ii * ka * s);
  }
  return rst;
}

MatrixXcd PhaseShiftDiff(const dcomp& ka, double angle, long m) {
  MatrixXcd rst(m, m);
  rst.setZero();
  double c = cos(angle), s = sin(angle);
  for (long i = 0; i < m / 2; i++) {
    rst(i, i) = ii * c * exp(ii * ka * c);
    rst(i + m / 2, i + m / 2) = ii * s * exp(ii * ka * s);
  }
  return rst;
}

double MinDet(const MatrixXcd& A, const MatrixXcd& B, double angle,
              int num_det = 5000) {
  // Return KL for minimum |A - SB| of Av = S Bv

  Eigen::JacobiSVD<MatrixXcd> svd = B.jacobiSvd(40);
  MatrixXcd u = svd.matrixU();
  Eigen::VectorXd det(num_det);
#ifdef NDEBUG
#pragma omp parallel for
#endif
  for (int i = 0; i < num_det; ++i)
    det(i) = std::abs(
        (u.adjoint() * A -
         u.adjoint() * PhaseShift(pi * i / num_det, angle, B.rows()) * B)
            .determinant());
  long index;
  det.minCoeff(&index);
  return double(index) / double(num_det);
}

VectorXcd NewtonEigen(const MatrixXcd& A, const MatrixXcd& Ap,
                      size_t max_iter = 1e4, double tol = 1e-10) {
  // VectorXcd v1 = VectorXcd::Random(A.cols());
  VectorXcd v1(A.cols());
  v1.setOnes();
  v1 /= v1.norm();
  VectorXcd v2(v1);

  size_t i = 0;
  while (i++ < max_iter) {
    v2 = A.jacobiSvd(40).solve(Ap * v1);
    v2 /= v2.norm();
    dcomp theta = v1.dot(v2);
    if ((v2 - theta * v1).norm() < tol * std::abs(theta)) {
      std::cout << i << std::endl;
      return v2;
    } else
      v1 = v2;
  }
  std::cout << "[mss:] Newton Iteration didn't converge." << std::endl;
  return v2;
}

// VectorXcd NewtonEigen(const MatrixXcd& A, const MatrixXcd& B, dcomp& kl,
//                       double angle, size_t max_iter = 1e2,
//                       double tol = 1e-10) {
//   // VectorXcd v1 = VectorXcd::Random(A.cols());

//   VectorXcd v1(A.cols());
//   v1.setOnes();
//   v1 /= v1.norm();
//   VectorXcd v2(v1);

//   size_t i = 0;
//   while (i++ < max_iter) {
//     MatrixXcd M(A - PhaseShift(exp(ii * kl * pi), angle, A.rows()) * B);
//     v2 = M.jacobiSvd(40).solve(M.transpose() * v1);
//     v2 /= v2.norm();
//     dcomp theta = v1.dot(v2);
//     if ((v2 - theta * v1).norm() < tol * std::abs(theta)) {
//       std::cout << i << std::endl;
//       return v2;
//     } else {
//       v1 = v2;
//       kl -= 1/theta;
//     }
//   }
//   std::cout << "[mss:] Inverse Iteration didn't converge." << std::endl;
//   return v2;
// }

}  // namespace mss

#include "Integrators.h"

#endif
