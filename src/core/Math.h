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
#include <cassert>
#include <complex>
#include <cstdlib>
#include <limits>
#include <string>

namespace mss {

typedef std::complex<double> dcomp;

const dcomp ii(0.0, 1.0);
const double ee(2.71828182845904523536);
const double pi(3.14159265358979323846);
const double epsilon(1e-14);

dcomp Jn(int n, const double& x) {
  return jn(n, x);
}

dcomp Hn(int n, const double& x) {
  return jn(n, x) + ii * yn(n, x);
}

// class Bessel {
//  public:
//   Bessel(BesselType f, int n, const double& k) : f(f), n(n), k(k) {}

//   double N() const { return n; }
//   dcomp operator()(const double& x) const { return f(n, x); }

//   // dcomp d(int j, const double& x, int i) const { return i > 0 ?
//   //   (d(j - 1, x, i - 1) - d(j + 1, x, i - 1)) / 2.0 : f(j, x);
//   //   }
//   //
//   // dcomp d(const double& x, int i) const { return d(n, x, i); }

//   dcomp d(const double& x) const { return 0.5 * (f(n - 1, x) - f(n + 1,
//   x)); } dcomp dd(const double& x) const {
//     return 0.25 * (f(n - 2, x) - 2.0 * f(n, x) + f(n + 2, x));
//   }

//   dcomp dr(const double& r) const { return k * d(k * r); }
//   dcomp ddr(const double& r) const { return k * k * dd(k * r); }
//   dcomp dt(const double& r) const { return dcomp(0, n) * f(n, k * r); }
//   dcomp ddt(const double& r) const { return dcomp(-n * n, 0) * f(n, k * r);
//   } dcomp drdt(const double& r) const { return dcomp(0, n) * k * d(k * r);
//   } dcomp Laplacian(const double& r) const {
//     return ddr(r) + dr(r) / r + ddt(r) / r / r;
//   }

//  private:
//   BesselType f;
//   int n;
//   double k;
// };

inline bool angEqu(const double& a, const double& b) {
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

}  // namespace mss

#endif
