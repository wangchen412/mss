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
#include "Matrix.h"

namespace mss {

// Cylindrical modes
// TODO Rename to CylinL and CylinT
template <typename T>
T ModeL(const CS*, const CS*, const EigenFunctor&, const Material&);
template <typename T>
T ModeT(const CS*, const CS*, const EigenFunctor&, const Material&);

template <>
StateIP ModeL<IP>(const CS* localCS, const CS* objCS, const EigenFunctor& f,
                  const Material& m) {
  /// Return the effect of the longitude mode in in-plane problems.

  // Position in the local CS, which is seen as a polar CS:
  const PosiVect pc = objCS->PositionIn(localCS);
  const PosiVect p = pc.Polar();
  double r = p.x;
  const CS cs(pc, p.y, localCS);

  // Displacement in the local CS:
  dcomp ur = f.dr(p);
  dcomp ut = f.dt_r(p);

  // Stress in the local CS:
  dcomp grr = f.ddr(p);
  dcomp gtt = f.dr_r(p) + f.ddt_rr(p);
  dcomp grt = 2.0 * (f.drdt_r(r) - f.dt_rr(p));
  StressIP t = m.C(grr, gtt, grt);

  // Normalized state in the objective CS.
  return StateIP(ur, ut, t, &cs).in(objCS);
}
template <>
StateIP ModeT<IP>(const CS* localCS, const CS* objCS, const EigenFunctor& f,
                  const class Material& m) {
  /// Return the effect of the transverse mode in in-plane problems.

  // Position in the local CS, which is seen as a polar CS:
  const PosiVect pc = objCS->PositionIn(localCS);
  const PosiVect p = pc.Polar();
  double r = p.x;
  const CS cs(pc, p.y, localCS);

  // Displacement in the local CS:
  dcomp ur = f.dt_r(r);
  dcomp ut = -f.dr(r);

  // Stress in the local CS:
  dcomp grr = f.drdt_r(r);
  dcomp gtt = f.dt_rr(r) - f.drdt_r(r);
  dcomp grt = f.ddt_rr(r) - f.ddr(r) + f.dr_r(r);
  StressIP t = m.C(grr, gtt, grt);

  // Normalized state in the objective CS.
  return StateIP(ur, ut, t, &cs).in(objCS);
}
template <>
StateAP ModeT<AP>(const CS* localCS, const CS* objCS, const EigenFunctor& f,
                  const class Material& m) {
  /// Return the effect of the transverse mode in antiplane problems.

  // Position in the local CS, which is seen as a polar CS:
  const PosiVect pc = objCS->PositionIn(localCS);
  const PosiVect p = pc.Polar();
  double r = p.x;
  const CS cs(pc, p.y, localCS);

  // Displacement in the local CS:
  DispAP w = f(p);

  // Stress in the local CS:
  dcomp gzr = f.dr(p);
  dcomp gzt = f.dt_r(p);
  StressAP t = m.C(gzr, gzt);

  // Normalized state in the objective CS.
  return StateAP(w, t, &cs).in(objCS);
}

// Plane modes
template <typename T>
T PlaneL(const Vector<double>&, const CS*, const Material&);
template <typename T>
T PlaneT(const Vector<double>&, const CS*, const Material&);

template <>
StateIP PlaneL<IP>(const Vector<double>& k, const CS* objCS,
                   const Material& m) {
  // The "local CS" in this case is the global CS, so the information of the
  // "source" is angle, which is carried by wave vector.
  double n = k.Length();
  dcomp e = exp(k * objCS->PositionGLB() * ii);

  // Displacement in the global CS:
  dcomp u = k.x * e / n;
  dcomp v = k.y * e / n;

  // Stress in the global CS:
  dcomp gxx = ii * k.x * u;
  dcomp gyy = ii * k.y * v;
  dcomp gxy = ii * (k.y * u + k.x * v);
  StressIP t = m.C(gxx, gyy, gxy);

  return StateIP(u, v, t).in(objCS);
}
template <>
StateIP PlaneT<IP>(const Vector<double>& k, const CS* objCS,
                   const Material& m) {
  // The "local CS" in this case is the global CS, and the information of the
  // "source" is angle, which is carried by wave vector.
  double n = k.Length();
  dcomp e = exp(k * objCS->PositionGLB() * ii);

  // Displacement in the global CS:
  dcomp u = -k.y * e / n;
  dcomp v = k.x * e / n;

  // Stress in the global CS:
  dcomp gxx = ii * k.x * u;
  dcomp gyy = ii * k.y * v;
  dcomp gxy = ii * (k.y * u + k.x * v);
  StressIP t = m.C(gxx, gyy, gxy);

  return StateIP(u, v, t).in(objCS);
}
template <>
StateAP PlaneT<AP>(const Vector<double>& k, const CS* objCS,
                   const Material& m) {
  // The "local CS" in this case is the global CS, and the information of the
  // "source" is angle, which is carried by wave vector.

  // Displacement in the global CS:
  dcomp w = exp(k * objCS->PositionGLB() * ii);

  // Stress in the global CS:
  dcomp gzx = ii * k.x * w;
  dcomp gzy = ii * k.y * w;
  StressAP t = m.C(gzx, gzy);

  return StateAP(w, t).in(objCS);
}

// The return should be a matrix of which size is 4x4 for in-plane and 2x2 for
// antiplane.
template <typename T>
MatrixNcd<T> GreenT(const CS* localCS, const CS* objCS, const Matrix* matrix);

// The matrix transform the boundary value vector of the source point to the
// one of the field point. TODO Need renaming.
template <>
Matrix2cd GreenT<AP>(const CS* localCS, const CS* objCS,
                     const Matrix* matrix) {
  Matrix2cd rst;

  CS X = objCS->inGLB(), Y = localCS->inGLB();
  // distance
  double r = (X.Position() - Y.Position()).Length();
  // n_j(y), n_l(x)
  double nxy = cos(Y.Angle()), nyy = sin(Y.Angle());
  double nxx = cos(X.Angle()), nyx = sin(X.Angle());
  // r_{,j} n_j(y)
  double rny = (Y.Position().x - X.Position().x) / r * nxy +
               (Y.Position().y - X.Position().y) / r * nyy;
  // r_{,j} n_j(x)
  double rnx = (X.Position().x - Y.Position().x) / r * nxx +
               (X.Position().y - Y.Position().y) / r * nyx;

  dcomp mu = matrix->Material().Mu_comp();
  dcomp k = matrix->KT_comp();

  dcomp H = Hn(0, k * r), Hd = -k * Hn(1, k * r), H2r = Hn(2, k * r);

  rst(0, 0) = -ii / 4 * Hd * rny;
  rst(0, 1) = ii / 4 * H / mu;
  rst(1, 0) = ii / 4 / r * mu * Hd * (nxx * nxy + nyx * nyy) -
              ii / 4 * k * k * mu * H2r * rnx * rny;
  rst(1, 1) = ii / 4 * Hd * rnx;

  return rst;
}

template <typename T>
Eigen::Matrix<dcomp, T::NumV, T::NumBv> GreenStateT(const CS* localCS,
                                                    const CS* objCS,
                                                    const Matrix* matrix);

// The matrix transform the boundary value vector of the source point to the
// vector with all the components of state of the field point.
// The state components are arranged as (w, tx, ty).
template <>
Eigen::Matrix<dcomp, 3, 2> GreenStateT<AP>(const CS* localCS, const CS* objCS,
                                           const Matrix* matrix) {
  Eigen::Matrix<dcomp, 3, 2> rst;

  CS X = objCS->inGLB(), Y = localCS->inGLB();
  // distance
  double r = (X.Position() - Y.Position()).Length();
  // n_j(y), n_l(x)
  double nxy = cos(Y.Angle()), nyy = sin(Y.Angle());
  double nxx = cos(X.Angle()), nyx = sin(X.Angle());
  double nxx2 = cos(X.Angle() + pi_2), nyx2 = sin(X.Angle() + pi_2);
  // r_{,j} n_j(y)
  double rny = (Y.Position().x - X.Position().x) / r * nxy +
               (Y.Position().y - X.Position().y) / r * nyy;
  // r_{,j} n_j(x)
  double rnx = (X.Position().x - Y.Position().x) / r * nxx +
               (X.Position().y - Y.Position().y) / r * nyx;
  double rnx2 = (X.Position().x - Y.Position().x) / r * nxx2 +
                (X.Position().y - Y.Position().y) / r * nyx2;

  dcomp mu = matrix->Material().Mu_comp();
  dcomp k = matrix->KT_comp();

  dcomp H = Hn(0, k * r), Hd = -k * Hn(1, k * r), H2r = Hn(2, k * r);

  rst(0, 0) = -ii / 4 * Hd * rny;
  rst(0, 1) = ii / 4 * H / mu;
  rst(1, 0) = ii / 4 / r * mu * Hd * (nxx * nxy + nyx * nyy) -
              ii / 4 * k * k * mu * H2r * rnx * rny;
  rst(1, 1) = ii / 4 * Hd * rnx;
  rst(2, 0) = ii / 4 / r * mu * Hd * (nxx2 * nxy + nyx2 * nyy) -
              ii / 4 * k * k * mu * H2r * rnx2 * rny;
  rst(2, 1) = ii / 4 * Hd * rnx2;

  return rst;
}

StateAP _planeWaveAP(const CS* objCS, double a, const Matrix* matrix) {
  const PosiVect p = objCS->PositionGLB();

  dcomp w = exp(ii * matrix->KT() * cos(a) * p.x +
                ii * matrix->KT() * sin(a) * p.y);

  dcomp gzx = ii * matrix->KT() * cos(a) * w;
  dcomp gzy = ii * matrix->KT() * sin(a) * w;
  StressAP t = matrix->Material().C(gzx, gzy);

  return StateAP(w, t).in(objCS);
}

}  // namespace mss

#endif
