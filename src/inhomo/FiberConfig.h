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
              const double& R, const Material& material, const Matrix* matrix,
              SolveMethod method)
      : ID_(ID),
        N_(N_max),
        NumCoeff_((2 * N_ + 1) * T::NumBv / 2),
        P_(P),
        r_(R),
        material_(material),
        kl_(matrix->Frequency() / material_.CL()),
        kt_(matrix->Frequency() / material_.CT()),
        matrix_(matrix),
        method_(method) {
    add_node();
    com_mat();
  }

  FiberConfig(const input::FiberConfig& input, const Matrix* matrix,
              SolveMethod method)
      : FiberConfig(input.ID, input.N_max, input.P, input.radius,
                    *input.material, matrix, method) {}

  virtual ~FiberConfig() { del_node(); }

  const Eigen::MatrixXcd& ColloMat() const;
  const Eigen::MatrixXcd& TransMat() const;
  size_t NumNode() const { return P_; }
  size_t NumBv() const { return P_ * T::NumBv; }
  size_t NumCoeff() const { return NumCoeff_; }
  int TopOrder() const { return N_; }

  const std::string& ID() const { return ID_; }
  const double& Radius() const { return r_; }
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

  Eigen::VectorXcd InciVect(const InciCPtrs<T>& inc) const;
  Eigen::VectorXcd Solve(const InciCPtrs<T>& inc) const;
  Eigen::VectorXcd CSolve(const InciCPtrs<T>& inc) const;
  Eigen::VectorXcd DSolve(const InciCPtrs<T>& inc) const;

 protected:
  const std::string ID_;           // The ID.
  const int N_;                    // The top order of the series.
  const int NumCoeff_;             // Number of the scattering coefficients.
  const size_t P_;                 // Number of the collocation points.
  const double r_;                 // Radius of the fiber.
  const class Material material_;  // Material of the fiber.
  const double kl_, kt_;           // Wave numbers of the fiber.
  const class Matrix* matrix_;     // The matrix.
  CSCPtrs node_;                   // Nodes.
  const SolveMethod method_;       // Solving method.
  Eigen::MatrixXcd Q_;             // Transform matrix.
  Eigen::MatrixXcd R_;             // Inner transform matrix.
  Eigen::MatrixXcd CQ_;            // Collocation matrix.
  bool QR_nc{true}, CQ_nc{true};

  void add_node();
  void del_node();
  void com_mat();
  void com_QR();
  void com_CQ();
  auto tm(int n) const;
};

template <typename T>
using FiberConfigPtrs = std::vector<FiberConfig<T>*>;

template <typename T>
using FiberConfigCPtrs = std::vector<const FiberConfig<T>*>;

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T>
const Eigen::MatrixXcd& FiberConfig<T>::ColloMat() const {
  assert(method_ == COLLOCATION);
  return CQ_;
}
template <typename T>
const Eigen::MatrixXcd& FiberConfig<T>::TransMat() const {
  assert(method_ == DFT);
  return Q_;
}

template <typename T>
inline void FiberConfig<T>::add_node() {
  node_.reserve(P_);
  for (size_t i = 0; i < P_; i++)
    node_.push_back(
        new CS(PosiVect(r_, i * pi2 / P_).Cartesian(), i * pi2 / P_));
}
template <typename T>
inline void FiberConfig<T>::del_node() {
  for (auto i : node_) delete i;
}

template <typename T>
inline void FiberConfig<T>::com_mat() {
  switch (method_) {
    case COLLOCATION:
      com_CQ();
      break;
    case DFT:
      com_QR();
      break;
    default:
      exit_error_msg({"Unknown method."});
  }
}

template <>
auto FiberConfig<StateIP>::tm(int n) const {
  const double& lm = Matrix()->Material().Lambda();
  const double& mm = Matrix()->Material().Mu();
  const double& lf = Material().Lambda();
  const double& mf = Material().Mu();
  const double lmm = lm + 2 * mm;
  const double lmf = lf + 2 * mf;

  BesselFunctor HL(Hn, n, KL_m(), r_), HT(Hn, n, KT_m(), r_);
  BesselFunctor JL(Jn, n, KL(), r_), JT(Jn, n, KT(), r_);

  dcomp HL_dr_r = HL.dr_r(r_), HL_rr = HL._rr(r_);
  dcomp HT_dr_r = HT.dr_r(r_), HT_rr = HT._rr(r_);
  dcomp JL_dr_r = JL.dr_r(r_), JL_rr = JL._rr(r_);
  dcomp JT_dr_r = JT.dr_r(r_), JT_rr = JT._rr(r_);

  Eigen::Matrix4cd Ti;

  Ti(0, 0) = -HL.dr(r_);
  Ti(0, 1) = -ii * n * HT._r(r_);
  Ti(0, 2) = JL.dr(r_);
  Ti(0, 3) = ii * n * JT._r(r_);

  Ti(1, 0) = -ii * n * HL._r(r_);
  Ti(1, 1) = HT.dr(r_);
  Ti(1, 2) = ii * n * JL._r(r_);
  Ti(1, 3) = -JT.dr(r_);

  Ti(2, 0) = -lmm * HL.ddr(r_) - lm * (HL_dr_r - n * n * HL_rr);
  Ti(2, 1) = 2 * ii * n * mm * (HT_rr - HT_dr_r);
  Ti(2, 2) = lmf * JL.ddr(r_) + lf * (JL_dr_r - n * n * JL_rr);
  Ti(2, 3) = 2 * ii * n * mf * (-JT_rr + JT_dr_r);

  Ti(3, 0) = -2 * ii * n * mm * (-HL_rr + HL_dr_r);
  Ti(3, 1) = -mm * (HT_dr_r - HT.ddr(r_) - n * n * HT_rr);
  Ti(3, 2) = 2 * ii * n * mf * (-JL_rr + JL_dr_r);
  Ti(3, 3) = mf * (JT_dr_r - JT.ddr(r_) - n * n * JT_rr);

  return Ti.inverse();
}
template <>
auto FiberConfig<StateAP>::tm(int n) const {
  BesselFunctor H(Hn, n, KT_m(), r_), J(Jn, n, KT(), r_);
  dcomp mH = H.dr(r_) * Matrix()->Material().Mu();
  dcomp mJ = J.dr(r_) * Material().Mu();
  Eigen::Matrix2cd Ti;
  Ti << -H(r_), J(r_), -mH, mJ;
  return Ti.inverse().eval();
}
template <typename T>
void FiberConfig<T>::com_QR() {
  Q_.resize(NumCoeff(), NumBv());
  R_.resize(NumCoeff(), NumBv());

  int N = T::NumBv / 2;

#ifdef NDEBUG
#pragma omp parallel for
#endif
  for (int n = -N_; n <= N_; n++) {
    auto t             = tm(n);
    Eigen::MatrixXcd g = IntMat(T::NumBv, NumNode(), n);

    Q_.block((n + N_) * N, 0, N, NumBv()) = t.block(0, 0, N, 2 * N) * g;
    R_.block((n + N_) * N, 0, N, NumBv()) = t.block(N, 0, N, 2 * N) * g;
  }
  QR_nc = false;
}
template <>
dcomp FiberConfig<StateAP>::TT(int n) const {
  BesselFunctor Jf(Jn, n, KT(), r_), Jm(Jn, n, KT_m(), r_);
  BesselFunctor Hm(Hn, n, KT_m(), r_);
  dcomp mJf = Jf.dr(r_) * Material().Mu();
  dcomp mJm = Jm.dr(r_) * Matrix()->Material().Mu();
  dcomp mHm = Hm.dr(r_) * Matrix()->Material().Mu();

  return (mJm * Hm(r_) - mHm * Jm(r_)) / (mJm * Jf(r_) - mJf * Jm(r_));
}
template <>
void FiberConfig<StateAP>::com_CQ() {
  CQ_.resize(NumBv(), NumCoeff());

#ifdef NDEBUG
#pragma omp parallel for
#endif
  for (int n = -N_; n <= N_; n++) {
    dcomp tn = TT(n);
    EigenFunctor J(Jn, n, KT(), r_), H(Hn, n, KT_m(), r_);
    for (size_t i = 0; i < P_; i++) {
      StateAP s = ModeT<StateAP>(nullptr, node_[i], J, Material()) * tn -
                  ModeT<StateAP>(nullptr, node_[i], H, Material_m());
      CQ_.block<2, 1>(2 * i, n + N_) = s.BV();
    }
  }
  CQ_nc = false;
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
Eigen::VectorXcd FiberConfig<T>::InciVect(const InciCPtrs<T>& inc) const {
  Eigen::VectorXcd rst(NumBv());
  rst.setZero();
  for (auto& i : inc) rst += i->EffectBV(Node());
  return rst;
}
template <typename T>
Eigen::VectorXcd FiberConfig<T>::Solve(const InciCPtrs<T>& inc) const {
  switch (method_) {
    case COLLOCATION:
      return CSolve(inc);
    case DFT:
      return DSolve(inc);
    default:
      error_msg({"Unknown method."});
      exit(EXIT_FAILURE);
  }
}
template <typename T>
Eigen::VectorXcd FiberConfig<T>::CSolve(const InciCPtrs<T>& inc) const {
  return CQ_.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV)
      .solve(InciVect(inc));
}
template <typename T>
Eigen::VectorXcd FiberConfig<T>::DSolve(const InciCPtrs<T>& inc) const {
  Eigen::VectorXcd tmp = Q_ * InciVect(inc);
  return Q_ * InciVect(inc);
}

}  // namespace mss

#endif
