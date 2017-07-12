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
#include <map>
#include <string>

namespace mss {

typedef std::complex<double> dcomp;

const dcomp ii(0.0, 1.0);
const double ee(2.71828182845904523536);
const double pi(3.14159265358979323846);
const double pi2(pi * 2);
const double pi_2(pi / 2);
const double epsilon(1e-14);

dcomp Jn(int n, const double& x) {
  return jn(n, x);
}

dcomp Hn(int n, const double& x) {
  return jn(n, x) + ii * yn(n, x);
}

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
