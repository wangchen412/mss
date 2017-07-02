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

#ifndef MSS_FUNCTORS_H
#define MSS_FUNCTORS_H

#include "Tensor.h"

namespace mss {

typedef decltype((Hn)) BesselFunc;

class EigenFunctor {
  /// The functors for the eigenfunction of scattering problem in cylindrical
  /// CS. Eigenfunction: Bessel(n, k * r) * exp(ii * n * t).

 public:
  EigenFunctor(BesselFunc f, int n, const double& k) : f(f), n(n), k(k) {}

  dcomp d(const double& x) const { return 0.5 * (f(n - 1, x) - f(n + 1, x)); }
  dcomp dd(const double& x) const {
    return 0.25 * (f(n - 2, x) - 2.0 * f(n, x) + f(n + 2, x));
  }

  dcomp operator()(const PosiVect& p) const {
    return f(n, k * p.x) * exp(ii * n * p.y);
  }
  dcomp dr(const PosiVect& p) const {
    return k * d(k * p.x) * exp(ii * n * p.y);
  }
  dcomp ddr(const PosiVect& p) const {
    return k * k * dd(k * p.x) * exp(ii * n * p.y);
  }
  dcomp dt(const PosiVect& p) const { return ii * n * (*this)(p); }
  dcomp ddt(const PosiVect& p) const { return ii * n * dt(p); }
  dcomp drdt(const PosiVect& p) const { return ii * n * dr(p); }

 private:
  BesselFunc f;
  double n, k;  // The order and the wave number.
};

}  // namespace mss

#endif
