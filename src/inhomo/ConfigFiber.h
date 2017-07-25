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

#ifndef MSS_CONFIGFIBER_H
#define MSS_CONFIGFIBER_H

#include "../core/Matrix.h"
#include "../core/State.h"
#include "../incident/Incident.h"
#include "../pre/Input.h"

namespace mss {

template <typename T>
class ConfigFiber {
 public:
  ConfigFiber(const std::string& ID, const size_t& N_max, const size_t& P,
              const double& R, const Material& material, const Matrix* matrix)
      : ID_(ID),
        N_(N_max),
        NoC_((2 * N_ + 1) * T::NoBV / 2),
        P_(P),
        R_(R),
        material_(material),
        kl_(matrix->Frequency() / material_.CL()),
        kt_(matrix->Frequency() / material_.CT()),
        matrix_(matrix),
        Q_(NoE(), NoC()) {
    add_node();
    compute_MatrixQ();
  }

  ConfigFiber(const input::ConfigFiber& input, const Matrix* matrix)
      : ConfigFiber(input.ID, input.N_max, input.P, input.radius,
                    *input.material, matrix) {}

  virtual ~ConfigFiber() { delete_node(); }

  const Eigen::MatrixXcd& TransMatrix() const { return Q_; }
  const double& CharLength() const { return R_; }
  size_t NoN() const { return P_; }
  size_t NoE() const { return P_ * T::NoBV; }
  size_t NoC() const { return NoC_; }
  int TopOrder() const { return N_; }

  const std::string& ID() const { return ID_; }
  const double& Radius() const { return R_; }
  const class Material& Material() const { return material_; }
  const double& KL() const { return kl_; }
  const double& KT() const { return kt_; }
  const class Matrix* Matrix() const { return matrix_; }
  const class Material& Material_m() const { return matrix_->Material(); }
  const double& KL_m() const { return matrix_->KL(); }
  const double& KT_m() const { return matrix_->KT(); }
  const CSCPtrs& Node() const { return node_; }

  T ModeT(const CS* localCS, const CS* objCS, const EigenFunctor& f,
          const class Material& mat) const;
  T ModeL(const CS* localCS, const CS* objCS, const EigenFunctor& f,
          const class Material& mat) const;

  // Functions for the factor T_n: B_n = T_n A_n.
  // tL is for longitude modes and tT is for transverse modes.
  dcomp TL(int n) const;  // TODO: T-matrix for in-plane problem.
  dcomp TT(int n) const;

  Eigen::VectorXcd InciVect(const InciPtrs<T>& incident) const;
  Eigen::VectorXcd Solve(const Eigen::VectorXcd& v) const;
  Eigen::VectorXcd Solve(const InciPtrs<T>& incident) const;

 protected:
  const std::string ID_;           // The ID.
  const int N_;                    // The top order of the series.
  const int NoC_;                  // Number of the scattering coefficients.
  const size_t P_;                 // Number of the collocation points.
  const double R_;                 // Radius of the fiber.
  const class Material material_;  // Material of the fiber.
  const double kl_, kt_;           // Wave numbers of the fiber.
  const class Matrix* matrix_;     // The matrix.
  CSCPtrs node_;                   // Collocation points.
  Eigen::MatrixXcd Q_;             // Transform matrix.

  void add_node();
  void delete_node();
  void compute_MatrixQ();
};

template <typename T>
using ConfigFiberPtrs = std::vector<ConfigFiber<T>*>;

template <typename T>
using ConfigFiberCPtrs = std::vector<const ConfigFiber<T>*>;

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T>
inline void ConfigFiber<T>::add_node() {
  node_.reserve(P_);
  for (size_t i = 0; i < P_; i++)
    node_.push_back(
        new CS(PosiVect(R_, i * pi2 / P_).Cartesian(), i * pi2 / P_));
}
template <typename T>
inline void ConfigFiber<T>::delete_node() {
  for (auto i : node_) delete i;
}
template <>
StateIP ConfigFiber<StateIP>::ModeL(const CS* localCS, const CS* objCS,
                                    const EigenFunctor& f,
                                    const class Material& m) const {
  /// Return the effect of the longitude mode in in-plane problems.

  // Position in the local CS, which is seen as a polar CS:
  CS cs(objCS->in(localCS));
  PosiVect p      = cs.Position().Polar();
  const double& r = p.x;

  // Displacement in the local CS:
  dcomp ur = f.dr(p);
  dcomp ut = f.dt(p) / r;

  // Stress in the local CS:
  dcomp grr  = f.ddr(p);
  dcomp gtt  = ur / r + f.ddt(p) / r / r;
  dcomp grt  = 2.0 * (f.drdt(r) / r - ut / r);
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
  PosiVect p      = cs.Position().Polar();
  const double& r = p.x;

  // Displacement in the local CS:
  dcomp ur = f.dt(r) / r;
  dcomp ut = -f.dr(r);

  // Stress in the local CS:
  dcomp grr  = f.drdt(r) / r;
  dcomp gtt  = ur / r - f.drdt(r) / r;
  dcomp grt  = f.ddt(r) / r / r - f.ddr(r) - ut / r;
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
  dcomp gzr  = f.dr(p);
  dcomp gzt  = f.dt(p);
  StressAP t = m.C(gzr, gzt);

  // Normalized state in the objective CS.
  return StateAP(w, t, &cs).in(objCS) / f(CharLength());
}
template <>
dcomp ConfigFiber<StateAP>::TT(int n) const {
  BesselFunctor Jf(Jn, n, KT()), Jm(Jn, n, KT_m()), Hm(Hn, n, KT_m());
  dcomp mJf = Jf.dr(R_) * Material().Mu();
  dcomp mJm = Jm.dr(R_) * Matrix()->Material().Mu();
  dcomp mHm = Hm.dr(R_) * Matrix()->Material().Mu();

  return (mJm - mHm * (Jm(R_) / Hm(R_))) / (mJm - mJf * (Jm(R_) / Jf(R_)));
}
template <>
void ConfigFiber<StateAP>::compute_MatrixQ() {
  for (int n = -N_; n <= N_; n++) {
    dcomp tn = TT(n);
    EigenFunctor J(Jn, n, KT()), H(Hn, n, KT_m());
    for (size_t i = 0; i < P_; i++) {
      StateAP s = ModeT(nullptr, node_[i], J, Material()) * tn -
                  ModeT(nullptr, node_[i], H, Material_m());
      Q_(i * 2, n + N_)     = s.Displacement().x;
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
template <typename T>
Eigen::VectorXcd ConfigFiber<T>::InciVect(const InciPtrs<T>& incident) const {
  Eigen::VectorXcd rst(NoE());
  rst.setZero();
  for (auto& i : incident) rst += i->EffectBV(Node());
  return rst;
}
template <typename T>
Eigen::VectorXcd ConfigFiber<T>::Solve(const Eigen::VectorXcd& v) const {
  return Q_.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(v);
}
template <typename T>
Eigen::VectorXcd ConfigFiber<T>::Solve(const InciPtrs<T>& incident) const {
  return Solve(InciVect(incident));
}

}  // namespace mss

#endif
