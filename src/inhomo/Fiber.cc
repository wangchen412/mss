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

#include "Fiber.h"

namespace mss {

template <>
StateAP Fiber<StateAP>::_modeT(const CS* objCS, EigenFunctor& f,
                               const Material& m) const {
  CS cs(objCS->in(LocalCS()));
  DispAP w = f(cs.Position());
  StressAP t = m.C(Geo(f, cs.Position()));
  dcomp norm = f(Config()->CharLength());

  return StateAP(w, t, &cs).in(objCS) / norm;
}
template <>
StateIP Fiber<StateIP>::_modeL(const CS* objCS, EigenFunctor& f,
                               const Material& m) const {
  CS cs(objCS->in(LocalCS()));

  auto uf = f.dr();
  auto vf = f.dt();

  dcomp ur = f.dr(r);
  dcomp ut = f.dt(r) / r;
  dcomp grr = f.ddr(r);
  dcomp gtt = ur / r + f.ddt(r) / r / r;
  dcomp grt = f.drdt(r) / r - ut / r;
  dcomp norm = exp(ii * f.N() * t) / f(Config()->CharLength());

  return StateIP({ur, ut}, m.C({grr, gtt, grt}), &cs).in(objCS)*norm;
}
template <>
StateIP Fiber<StateIP>::_modeT(const CS* objCS, const Bessel& f,
                               const Material& m) const {
  CS cs(objCS->in(LocalCS()));
  const double& r = cs.Position().Polar().x;
  const double& t = cs.Position().Polar().y;

  dcomp ur = f.dt(r) / r;
  dcomp ut = -f.dr(r);

  dcomp grr = f.drdt(r) / r;

  dcomp trr = 2 * m.mu * (f.drdt(r) - f.dt(r) / r);
  dcomp ttt = 2 * m.mu * (f.dt(r) / r - f.drdt(r));
  dcomp trt = m.mu * (f.dr(r) / r - f.ddr(r) - f.ddt(r));
  dcomp norm = f(Config()->CharLength());

  return StateIP(ur, ut, trr, ttt, trt, &cs).in(objCS)*exp(ii * f.N() * t) /
         norm;
}

template <>
StateAP Fiber<StateAP>::ScatterModeT(const CS* objCS, int n) const {
  return _modeT(objCS, Bessel(Hn, n, kt_m), Config()->Matrix()->Material());
}
template <>
StateAP Fiber<StateAP>::InnerModeT(const CS* objCS, int n) const {
  return _modeT(objCS, Bessel(Jn, n, kt_f), Config()->Material());
}

}  // namespace mss
