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

#include "State.h"

namespace mss {

// typedef decltype((Hn)) BesselFunc;
typedef dcomp (*BesselFunc)(int, const double&);

class BesselFunctor {
  /// The functor for the Bessel functions: Bessel(n, k * r).

 public:
  BesselFunctor(const BesselFunc f, int n, const double& k, const double& R)
      : f(f), n(n), k(k), norm(f(n, k * R)) {}

  dcomp operator()(const double& r) const { return f(n, k * r) / norm; }
  dcomp dr(const double& r) const { return k * d(k * r) / norm; }
  dcomp ddr(const double& r) const { return k * k * dd(k * r) / norm; }

 private:
  const BesselFunc f;
  const int n;
  const double k;
  const dcomp norm;

  // The derivatives for dBessel(x) / dx:
  dcomp d(const double& x) const { return 0.5 * (f(n - 1, x) - f(n + 1, x)); }
  dcomp dd(const double& x) const {
    return 0.25 * (f(n - 2, x) - 2.0 * f(n, x) + f(n + 2, x));
  }
};

class ExpFunctor {
  /// The functor for the exponential function: exp(ii * n * t).

 public:
  ExpFunctor(int n) : n(n) {}

  dcomp operator()(const double& r) const { return exp(ii * n * r); }
  dcomp dt(const double& r) const { return ii * n * (*this)(r); }
  dcomp ddt(const double& r) const { return ii * n * dt(r); }

 private:
  double n;
};

class EigenFunctor {
  /// The functor for the eigenfunction of scattering problem in cylindrical
  /// CS. Eigenfunction: Bessel(n, k * r) * exp(ii * n * t).

 public:
  EigenFunctor(const BesselFunc f, int n, const double& k, const double& R)
      : f(f, n, k, R), g(n), n(n), k(k) {}

  dcomp operator()(const PosiVect& p) const { return f(p.x) * g(p.y); }
  dcomp dr(const PosiVect& p) const { return f.dr(p.x) * g(p.y); }
  dcomp dt(const PosiVect& p) const { return f(p.x) * g.dt(p.y); }
  dcomp ddr(const PosiVect& p) const { return f.ddr(p.x) * g(p.y); }
  dcomp ddt(const PosiVect& p) const { return f(p.x) * g.ddt(p.y); }
  dcomp drdt(const PosiVect& p) const { return f.dr(p.x) * g.dt(p.y); }

  int N() const { return n; }
  const double& K() const { return k; }

 private:
  const BesselFunctor f;
  const ExpFunctor g;

  const int n;
  const double k;
};

}  // namespace mss

#endif
