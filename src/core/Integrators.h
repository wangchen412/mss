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

#ifndef MSS_INTEGRATORS_H
#define MSS_INTEGRATORS_H

#include "Math.h"

namespace mss {

inline VectorXcd DFT_v(size_t p, int n) {
  VectorXcd g(p), h(p);
  double d = pi2 / p;
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

  VectorXcd rst(p);
  rst(0) = g(0) + h(p - 1);
  for (size_t i = 1; i < p; i++) rst(i) = h(i - 1) + g(i);
  rst /= pi2;

  return rst;
}

inline MatrixXcd DFT_m(int m, size_t p, int n) {
  VectorXcd v = DFT_v(p, n);
  MatrixXcd rst(m, m * p);
  rst.setZero();
  for (int i = 0; i < m; i++)
    for (size_t j = 0; j < p; j++) rst(i, i + m * j) = v(j);
  return rst;
}

// template <int N>
// class GaussQuad {
//  public:
//   template <typename F>
//   auto operator()(const F& f, double a = -1.0, double b = 1.0,
//                               decltype(f(1.0)) rst = 0) {
//     double c = (b - a) / 2, d = (b + a) / 2;
//     for (int i = 0; i < N; i++) rst += l_.weight(i) * f(c * l_.root(i) +
//     d); return rst * c;
//   }

//  private:
//   static LegendreRoot<N> l_;
// };

// template <int N>
// LegendreRoot<N> GaussQuad<N>::l_;

}  // namespace mss

#endif
