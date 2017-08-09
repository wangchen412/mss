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

#ifndef MSS_FIBERCONFIG_H
#define MSS_FIBERCONFIG_H

#include "../core/Matrix.h"
#include "../core/Modes.h"
#include "../incident/Incident.h"
#include "../pre/Input.h"

namespace mss {

template <typename T>
class FiberConfig {
 public:
  FiberConfig(const std::string& ID, const size_t& N_max, const size_t& P,
              const double& R, const Material& material, const Matrix* matrix)
      : ID_(ID),
        N_(N_max),
        NumCoeff_((2 * N_ + 1) * T::NumBv / 2),
        P_(P),
        R_(R),
        material_(material),
        kl_(matrix->Frequency() / material_.CL()),
        kt_(matrix->Frequency() / material_.CT()),
        matrix_(matrix),
        CQ_(NumBv(), NumCoeff()) {
    add_node();
    com_CQ();
  }

  FiberConfig(const input::FiberConfig& input, const Matrix* matrix)
      : FiberConfig(input.ID, input.N_max, input.P, input.radius,
                    *input.material, matrix) {}

  virtual ~FiberConfig() { del_node(); }

  const Eigen::MatrixXcd& ColloMat() const { return CQ_; }
  size_t NumNode() const { return P_; }
  size_t NumBv() const { return P_ * T::NumBv; }
  size_t NumCoeff() const { return NumCoeff_; }
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
  const int NumCoeff_;             // Number of the scattering coefficients.
  const size_t P_;                 // Number of the collocation points.
  const double R_;                 // Radius of the fiber.
  const class Material material_;  // Material of the fiber.
  const double kl_, kt_;           // Wave numbers of the fiber.
  const class Matrix* matrix_;     // The matrix.
  CSCPtrs node_;                   // Collocation points.
  Eigen::MatrixXcd CQ_;            // Collocation matrix.

  void add_node();
  void del_node();
  void com_CQ();
};

template <typename T>
using FiberConfigPtrs = std::vector<FiberConfig<T>*>;

template <typename T>
using FiberConfigCPtrs = std::vector<const FiberConfig<T>*>;

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T>
inline void FiberConfig<T>::add_node() {
  node_.reserve(P_);
  for (size_t i = 0; i < P_; i++)
    node_.push_back(
        new CS(PosiVect(R_, i * pi2 / P_).Cartesian(), i * pi2 / P_));
}
template <typename T>
inline void FiberConfig<T>::del_node() {
  for (auto i : node_) delete i;
}
template <>
dcomp FiberConfig<StateAP>::TT(int n) const {
  BesselFunctor Jf(Jn, n, KT(), R_), Jm(Jn, n, KT_m(), R_);
  BesselFunctor Hm(Hn, n, KT_m(), R_);
  dcomp mJf = Jf.dr(R_) * Material().Mu();
  dcomp mJm = Jm.dr(R_) * Matrix()->Material().Mu();
  dcomp mHm = Hm.dr(R_) * Matrix()->Material().Mu();

  return (mJm * Hm(R_) - mHm * Jm(R_)) / (mJm * Jf(R_) - mJf * Jm(R_));
}
template <>
void FiberConfig<StateAP>::com_CQ() {
#ifdef NDEBUG
#pragma omp parallel for
#endif

  for (int n = -N_; n <= N_; n++) {
    dcomp tn = TT(n);
    EigenFunctor J(Jn, n, KT(), R_), H(Hn, n, KT_m(), R_);
    for (size_t i = 0; i < P_; i++) {
      StateAP s = ModeT<StateAP>(nullptr, node_[i], J, Material()) * tn -
                  ModeT<StateAP>(nullptr, node_[i], H, Material_m());
      CQ_.block<2, 1>(2 * i, n + N_) = s.BV();
    }
  }
}
template <>
dcomp FiberConfig<StateIP>::TL(int) const {
  return 0;  // TODO: tL of in-plane problem
}
template <>
dcomp FiberConfig<StateIP>::TT(int) const {
  return 0;  // TODO: tT of in-plane problem
}
template <>
void FiberConfig<StateIP>::com_CQ() {
  // TODO: tT of in-plane problem
}
template <typename T>
Eigen::VectorXcd FiberConfig<T>::InciVect(const InciPtrs<T>& incident) const {
  Eigen::VectorXcd rst(NumBv());
  rst.setZero();
  for (auto& i : incident) rst += i->EffectBV(Node());
  return rst;
}
template <typename T>
Eigen::VectorXcd FiberConfig<T>::Solve(const Eigen::VectorXcd& v) const {
  return CQ_.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(v);
}
template <typename T>
Eigen::VectorXcd FiberConfig<T>::Solve(const InciPtrs<T>& incident) const {
  return Solve(InciVect(incident));
}

}  // namespace mss

#endif
