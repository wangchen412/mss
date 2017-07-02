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

// ---------------------------------------------------------------------------
// In-plane modes:

template <>
StateIP Fiber<StateIP>::modeL(const CS* objCS, const EigenFunctor& f,
                              const Material& m) const {
  /// Return the effect of the longitude mode in in-plane problems.

  // Position in the local CS, which is seen as a polar CS:
  CS cs(objCS->in(LocalCS()));
  PosiVect p = cs.Position().Polar();
  const double& r = p.x;

  // Displacement in the local CS:
  dcomp ur = f.dr(p);
  dcomp ut = f.dt(p) / r;

  // Stress in the local CS:
  dcomp grr = f.ddr(p);
  dcomp gtt = ur / r + f.ddt(p) / r / r;
  dcomp grt = 2.0 * (f.drdt(r) / r - ut / r);
  StressIP t = m.C(grr, gtt, grt);

  // Normalized state in the objective CS.
  return StateIP(ur, ut, t, &cs).in(objCS) / f(Config()->CharLength());
}
template <>
StateIP Fiber<StateIP>::modeT(const CS* objCS, const EigenFunctor& f,
                              const Material& m) const {
  /// Return the effect of the transverse mode in in-plane problems.

  // Position in the local CS, which is seen as a polar CS:
  CS cs(objCS->in(LocalCS()));
  PosiVect p = cs.Position().Polar();
  const double& r = p.x;

  // Displacement in the local CS:
  dcomp ur = f.dt(r) / r;
  dcomp ut = -f.dr(r);

  // Stress in the local CS:
  dcomp grr = f.drdt(r) / r;
  dcomp gtt = ur / r - f.drdt(r) / r;
  dcomp grt = f.ddt(r) / r / r - f.ddr(r) - ut / r;
  StressIP t = m.C(grr, gtt, grt);

  // Normalized state in the objective CS.
  return StateIP(ur, ut, t, &cs).in(objCS) / f(Config()->CharLength());
}

// ---------------------------------------------------------------------------
// Antiplane modes:

template <>
StateAP Fiber<StateAP>::modeL(const CS* objCS, const EigenFunctor&,
                              const Material&) const {
  /// Return zero state for the longitude mode in antiplane problems.

  return StateAP(0, 0, 0, objCS);
}
template <>
StateAP Fiber<StateAP>::modeT(const CS* objCS, const EigenFunctor& f,
                              const Material& m) const {
  /// Return the effect of the transverse mode in antiplane problems.

  // Position in the local CS, which is seen as a polar CS:
  CS cs(objCS->in(LocalCS()));
  PosiVect p = cs.Position().Polar();

  // Displacement in the local CS:
  DispAP w = f(p);

  // Stress in the local CS:
  dcomp gzr = f.dr(p);
  dcomp gzt = f.dt(p);
  StressAP t = m.C(gzr, gzt);

  // Normalized state in the objective CS.
  return StateAP(w, t, &cs).in(objCS) / f(Config()->CharLength());
}

}  // namespace mss
