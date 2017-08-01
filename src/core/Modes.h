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

#ifndef MSS_MODES_H
#define MSS_MODES_H

#include "CS.h"
#include "Functors.h"
#include "Material.h"

namespace mss {

template <typename T>
inline T ModeL(const CS*, const CS*, const EigenFunctor&, const Material& m);

template <typename T>
inline T ModeT(const CS*, const CS*, const EigenFunctor&, const Material& m);

template <>
inline StateIP ModeL<StateIP>(const CS* localCS, const CS* objCS,
                              const EigenFunctor& f, const Material& m) {
  /// Return the effect of the longitude mode in in-plane problems.

  // Position in the local CS, which is seen as a polar CS:
  const PosiVect pc = objCS->PositionIn(localCS);
  const PosiVect p  = pc.Polar();
  const double& r   = p.x;
  const CS cs(pc, p.y, localCS);

  // TODO when r < epsilon

  // Displacement in the local CS:
  dcomp ur = f.dr(p);
  dcomp ut = f.dt(p) / r;

  // Stress in the local CS:
  dcomp grr  = f.ddr(p);
  dcomp gtt  = ur / r + f.ddt(p) / r / r;
  dcomp grt  = 2.0 * (f.drdt(r) / r - ut / r);
  StressIP t = m.C(grr, gtt, grt);

  // Normalized state in the objective CS.
  return StateIP(ur, ut, t, &cs).in(objCS);
}

template <>
inline StateIP ModeT<StateIP>(const CS* localCS, const CS* objCS,
                              const EigenFunctor& f,
                              const class Material& m) {
  /// Return the effect of the transverse mode in in-plane problems.

  // Position in the local CS, which is seen as a polar CS:
  const PosiVect pc = objCS->PositionIn(localCS);
  const PosiVect p  = pc.Polar();
  const double& r   = p.x;
  const CS cs(pc, p.y, localCS);

  // TODO when r < epsilon

  // Displacement in the local CS:
  dcomp ur = f.dt(r) / r;
  dcomp ut = -f.dr(r);

  // Stress in the local CS:
  dcomp grr  = f.drdt(r) / r;
  dcomp gtt  = ur / r - f.drdt(r) / r;
  dcomp grt  = f.ddt(r) / r / r - f.ddr(r) - ut / r;
  StressIP t = m.C(grr, gtt, grt);

  // Normalized state in the objective CS.
  return StateIP(ur, ut, t, &cs).in(objCS);
}

template <>
inline StateAP ModeT<StateAP>(const CS* localCS, const CS* objCS,
                              const EigenFunctor& f,
                              const class Material& m) {
  /// Return the effect of the transverse mode in antiplane problems.

  // Position in the local CS, which is seen as a polar CS:
  const PosiVect pc = objCS->PositionIn(localCS);
  const PosiVect p  = pc.Polar();
  const double& r   = p.x;
  const CS cs(pc, p.y, localCS);

  // Displacement in the local CS:
  DispAP w = f(p);

  // Stress in the local CS:
  dcomp gzr = f.dr(p), gzt;
  if (r > epsilon)
    gzt = f.dt(p) / r;
  else
    gzt = std::abs(f.N()) == 1 ? 0.5 * ii * f.K() * exp(f.N() * p.y * ii) : 0;
  StressAP t = m.C(gzr, gzt);

  // Normalized state in the objective CS.
  return StateAP(w, t, &cs).in(objCS);
}

}  // namespace mss

#endif
