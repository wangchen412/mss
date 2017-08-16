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

#include "Functors.h"
#include "Material.h"

namespace mss {

template <typename T>
inline T ModeL(const CS*, const CS*, const EigenFunctor&, const Material&);
template <typename T>
inline T ModeT(const CS*, const CS*, const EigenFunctor&, const Material&);

template <>
inline StateIP ModeL<StateIP>(const CS* localCS, const CS* objCS,
                              const EigenFunctor& f, const Material& m) {
  /// Return the effect of the longitude mode in in-plane problems.

  // Position in the local CS, which is seen as a polar CS:
  const PosiVect pc = objCS->PositionIn(localCS);
  const PosiVect p  = pc.Polar();
  const double& r   = p.x;
  const CS cs(pc, p.y, localCS);

  // Displacement in the local CS:
  dcomp ur = f.dr(p);
  dcomp ut = f.dt_r(p);

  // Stress in the local CS:
  dcomp grr  = f.ddr(p);
  dcomp gtt  = f.dr_r(p) + f.ddt_rr(p);
  dcomp grt  = 2.0 * (f.drdt_r(r) - f.dt_rr(p));
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

  // Displacement in the local CS:
  dcomp ur = f.dt_r(r);
  dcomp ut = -f.dr(r);

  // Stress in the local CS:
  dcomp grr  = f.drdt_r(r);
  dcomp gtt  = f.dt_rr(r) - f.drdt_r(r);
  dcomp grt  = f.ddt_rr(r) - f.ddr(r) + f.dr_r(r);
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
  dcomp gzr  = f.dr(p);
  dcomp gzt  = f.dt_r(p);
  StressAP t = m.C(gzr, gzt);

  // Normalized state in the objective CS.
  return StateAP(w, t, &cs).in(objCS);
}

}  // namespace mss

#endif
