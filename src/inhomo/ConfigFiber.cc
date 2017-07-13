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

#include "ConfigFiber.h"

namespace mss {

template <>
StateIP ConfigFiber<StateIP>::ModeL(const CS* localCS, const CS* objCS,
                                    const EigenFunctor& f,
                                    const class Material& m) const {
  /// Return the effect of the longitude mode in in-plane problems.

  // Position in the local CS, which is seen as a polar CS:
  CS cs(objCS->in(localCS));
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
  return StateIP(ur, ut, t, &cs).in(objCS) / f(CharLength());
}
template <>
StateIP ConfigFiber<StateIP>::ModeT(const CS* localCS, const CS* objCS,
                                    const EigenFunctor& f,
                                    const class Material& m) const {
  /// Return the effect of the transverse mode in in-plane problems.

  // Position in the local CS, which is seen as a polar CS:
  CS cs(objCS->in(localCS));
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
  return StateIP(ur, ut, t, &cs).in(objCS) / f(CharLength());
}
template <>
StateAP ConfigFiber<StateAP>::ModeT(const CS* localCS, const CS* objCS,
                                    const EigenFunctor& f,
                                    const class Material& m) const {
  /// Return the effect of the transverse mode in antiplane problems.

  // Position in the local CS, which is seen as a polar CS:
  CS cs(objCS->in(localCS));
  PosiVect p = cs.Position().Polar();

  // Displacement in the local CS:
  DispAP w = f(p);

  // Stress in the local CS:
  dcomp gzr = f.dr(p);
  dcomp gzt = f.dt(p);
  StressAP t = m.C(gzr, gzt);

  // Normalized state in the objective CS.
  return StateAP(w, t, &cs).in(objCS) / f(CharLength());
}

template <>
dcomp ConfigFiber<StateAP>::TT(int n) const {
  BesselFunctor Jf(Jn, n, KT()), Jm(Jn, n, KT_m()), Hm(Hn, n, KT_m());
  dcomp muJR_f = Jf.dr(R_) * Material().Mu();
  dcomp muJR_m = Jm.dr(R_) * Matrix()->Material().Mu();
  dcomp muHR_m = Hm.dr(R_) * Matrix()->Material().Mu();
  return (Jf(R_) - Hm(R_)) * muJR_m / (muJR_f - muHR_m) / Jm(R_);
}
template <>
void ConfigFiber<StateAP>::compute_MatrixQ() {
  for (int n = -N_; n <= N_; n++) {
    dcomp tn = TT(n);
    EigenFunctor J(Jn, n, KT()), H(Hn, n, KT_m());
    for (size_t i = 0; i < P_; i++) {
      StateAP s = ModeT(nullptr, &node_[i], J, Material()) * tn -
                  ModeT(nullptr, &node_[i], H, Material_m());
      Q_(i * 2, n + N_) = s.Displacement().x;
      Q_(i * 2 + 1, n + N_) = s.Stress().x;
    }
  }
}
template <>
dcomp ConfigFiber<StateIP>::TL(int) const {
  return 0;  // TODO: tL of in-plane problem
}
template <>
dcomp ConfigFiber<StateIP>::TT(int) const {
  return 0;  // TODO: tT of in-plane problem
}
template <>
void ConfigFiber<StateIP>::compute_MatrixQ() {
  // TODO: tT of in-plane problem
}

}  // namespace mss
