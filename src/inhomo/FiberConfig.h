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
  FiberConfig(const std::string& ID, size_t N_max, size_t P, double R,
              const Material& material, const Matrix* matrix)
      : ID_(ID),
        N_(N_max),
        NumCoeff_((2 * N_ + 1) * T::NumBv / 2),
        P_(P),
        r_(R),
        material_(material),
        kl_(matrix->Frequency() / material_.CL()),
        kt_(matrix->Frequency() / material_.CT()),
        matrix_(matrix) {
    add_node();
  }

  FiberConfig(const input::FiberConfig& input, const Matrix* matrix)
      : FiberConfig(input.ID, input.N_max, input.P, input.radius,
                    *input.material, matrix) {}

  virtual ~FiberConfig() { del_node(); }

  const MatrixXcd& ColloMat();   // Collocation of boundary values.
  const MatrixXcd& ColloDMat();  // Collocation of boundary displacement.
  const MatrixXcd& TransMat();
  const MatrixXcd& RefraMat();
  size_t NumNode() const { return P_; }
  size_t NumBv() const { return P_ * T::NumBv; }
  size_t NumDv() const { return P_ * T::NumDv; }
  size_t NumCoeff() const { return NumCoeff_; }
  int TopOrder() const { return N_; }

  const std::string& ID() const { return ID_; }
  double Radius() const { return r_; }
  const class Material& Material() const { return material_; }
  double KL() const { return kl_; }
  double KT() const { return kt_; }
  const class Matrix* Matrix() const { return matrix_; }
  const class Material& Material_m() const { return matrix_->Material(); }
  double KL_m() const { return matrix_->KL(); }
  double KT_m() const { return matrix_->KT(); }
  const CSCPtrs& Node() const { return node_; }

  // Functions for the factor T_n: B_n = T_n A_n.
  // tL is for longitude modes and tT is for transverse modes.
  dcomp TL(int n) const;  // TODO: T-matrix for in-plane problem.
  dcomp TT(int n) const;
  dcomp T_sc_in_T(int n) const;

  VectorXcd Solve(const VectorXcd& incBv, SolveMethod method);
  VectorXcd CSolve(const VectorXcd& incBv);
  VectorXcd DSolve(const VectorXcd& incBv);

  VectorXcd IncVec(const InciCPtrs<T>& inc) const;
  VectorXcd Solve(const InciCPtrs<T>& inc, SolveMethod method);
  VectorXcd CSolve(const InciCPtrs<T>& inc);
  VectorXcd DSolve(const InciCPtrs<T>& inc);

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
  MatrixXcd Q_;                    // Transform matrix.
  MatrixXcd R_;                    // Inner transform matrix.
  MatrixXcd CQ_;                   // Collocation matrix.
  MatrixXcd CQD_;  // Collocation matrix of boundary displacement.
  bool qr_computed_{false}, cq_computed_{false}, cqd_computed_{false};

  void add_node();
  void del_node();
  void com_QR();
  MatrixNcd<T> tm(int n) const;
};

template <typename T>
using FiberConfigPtrs = std::vector<FiberConfig<T>*>;

template <typename T>
using FiberConfigCPtrs = std::vector<const FiberConfig<T>*>;

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T>
const MatrixXcd& FiberConfig<T>::ColloMat() {
  if (cq_computed_) return CQ_;

  CQ_.resize(NumBv(), NumCoeff());

#ifdef NDEBUG
#pragma omp parallel for
#endif
  for (int n = -N_; n <= N_; n++) {
    // For antiplane only. TODO: in-plane.
    dcomp tn = TT(n);
    EigenFunctor J(Jn, n, KT(), r_), H(Hn, n, KT_m(), r_);
    for (size_t i = 0; i < P_; i++) {
      StateAP s = ModeT<AP>(nullptr, node_[i], J, Material()) * tn -
                  ModeT<AP>(nullptr, node_[i], H, Material_m());
      CQ_.block<T::NumBv, 1>(T::NumBv * i, n + N_) = s.Bv();
    }
  }

  cq_computed_ = true;

  return CQ_;
}
template <typename T>
const MatrixXcd& FiberConfig<T>::ColloDMat() {
  if (cqd_computed_) return CQD_;

  CQD_.resize(NumDv(), NumCoeff());

#ifdef NDEBUG
#pragma omp parallel for
#endif
  for (int n = -N_; n <= N_; n++) {
    // For antiplane only. TODO: in-plane.
    dcomp tn = TT(n);
    EigenFunctor J(Jn, n, KT(), r_), H(Hn, n, KT_m(), r_);
    for (size_t i = 0; i < P_; i++) {
      StateAP s = ModeT<AP>(nullptr, node_[i], J, Material()) * tn -
                  ModeT<AP>(nullptr, node_[i], H, Material_m());
      CQD_.block<T::NumDv, 1>(T::NumDv * i, n + N_) = s.Dv();
    }
  }

  cqd_computed_ = true;

  return CQD_;
}
template <typename T>
const MatrixXcd& FiberConfig<T>::TransMat() {
  if (qr_computed_) return Q_;
  com_QR();
  qr_computed_ = true;
  return Q_;
}
template <typename T>
const MatrixXcd& FiberConfig<T>::RefraMat() {
  if (qr_computed_) return R_;
  com_QR();
  qr_computed_ = true;
  return R_;
}
template <typename T>
void FiberConfig<T>::add_node() {
  node_.reserve(P_);
  for (size_t i = 0; i < P_; i++)
    node_.push_back(
        new CS(PosiVect(r_, i * pi2 / P_).Cartesian(), i * pi2 / P_));
}
template <typename T>
void FiberConfig<T>::del_node() {
  for (auto i : node_) delete i;
}
template <>
Matrix4cd FiberConfig<IP>::tm(int n) const {
  double lm = Matrix()->Material().Lambda();
  double mm = Matrix()->Material().Mu();
  double lf = Material().Lambda();
  double mf = Material().Mu();
  const double lmm = lm + 2 * mm;
  const double lmf = lf + 2 * mf;

  BesselFunctor HL(Hn, n, KL_m(), r_), HT(Hn, n, KT_m(), r_);
  BesselFunctor JL(Jn, n, KL(), r_), JT(Jn, n, KT(), r_);

  dcomp HL_dr_r = HL.dr_r(r_), HL_rr = HL._rr(r_);
  dcomp HT_dr_r = HT.dr_r(r_), HT_rr = HT._rr(r_);
  dcomp JL_dr_r = JL.dr_r(r_), JL_rr = JL._rr(r_);
  dcomp JT_dr_r = JT.dr_r(r_), JT_rr = JT._rr(r_);

  Matrix4cd Ti;

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
Matrix2cd FiberConfig<AP>::tm(int n) const {
  BesselFunctor H(Hn, n, KT_m(), r_), J(Jn, n, KT(), r_);
  dcomp mH = H.dr(r_) * Matrix()->Material().Mu();
  dcomp mJ = J.dr(r_) * Material().Mu();
  Matrix2cd Ti;
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
    auto t = tm(n);
    MatrixXcd g = DFT_m(T::NumBv, NumNode(), n);

    Q_.block((n + N_) * N, 0, N, NumBv()) = t.block(0, 0, N, 2 * N) * g;
    R_.block((n + N_) * N, 0, N, NumBv()) = t.block(N, 0, N, 2 * N) * g;
  }
}
template <>
dcomp FiberConfig<AP>::TT(int n) const {
  BesselFunctor Jf(Jn, n, KT(), r_), Jm(Jn, n, KT_m(), r_);
  BesselFunctor Hm(Hn, n, KT_m(), r_);
  dcomp mJf = Jf.dr(r_) * Material().Mu();
  dcomp mJm = Jm.dr(r_) * Matrix()->Material().Mu();
  dcomp mHm = Hm.dr(r_) * Matrix()->Material().Mu();

  return (mJm * Hm(r_) - mHm * Jm(r_)) / (mJm * Jf(r_) - mJf * Jm(r_));
}
template <>
dcomp FiberConfig<IP>::TL(int) const {
  return 0;  // TODO: tL of in-plane problem
}
template <>
dcomp FiberConfig<IP>::TT(int) const {
  return 0;  // TODO: tT of in-plane problem
}
template <>
dcomp FiberConfig<AP>::T_sc_in_T(int n) const {
  // Transformation from the nth scattering wave expansion coefficient to the
  // nth incident wave expansion coefficient.

  BesselFunctor Jf(Jn, n, KT(), r_), Jm(Jn, n, KT_m(), r_);
  BesselFunctor Hm(Hn, n, KT_m(), r_);
  dcomp mJf = Jf.dr(r_) * Material().Mu();
  dcomp mJm = Jm.dr(r_) * Matrix()->Material().Mu();
  dcomp mHm = Hm.dr(r_) * Matrix()->Material().Mu();

  return (mJf * Hm(r_) - mHm * Jf(r_)) / (mJm * Jf(r_) - mJf * Jm(r_));
}
template <>
dcomp FiberConfig<IP>::T_sc_in_T(int) const {
  return 0;  // TODO: in-plane problem
}

template <typename T>
VectorXcd FiberConfig<T>::Solve(const VectorXcd& incBv, SolveMethod method) {
  switch (method) {
    case COLLOCATION:
      return CSolve(incBv);
    case DFT:
      return DSolve(incBv);
    default:
      error_msg({"Unknown method."});
      exit(EXIT_FAILURE);
  }
}
template <typename T>
VectorXcd FiberConfig<T>::CSolve(const VectorXcd& incBv) {
  return ColloMat()
      .jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV)
      .solve(incBv);
}
template <typename T>
VectorXcd FiberConfig<T>::DSolve(const VectorXcd& incBv) {
  return TransMat() * incBv;
}
template <typename T>
VectorXcd FiberConfig<T>::IncVec(const InciCPtrs<T>& inc) const {
  VectorXcd rst(NumBv());
  rst.setZero();
  for (auto& i : inc) rst += i->EffectBv(Node());
  return rst;
}
template <typename T>
VectorXcd FiberConfig<T>::Solve(const InciCPtrs<T>& inc, SolveMethod method) {
  return Solve(inciVect(inc), method);
}
template <typename T>
VectorXcd FiberConfig<T>::CSolve(const InciCPtrs<T>& inc) {
  return CSolve(inciVect(inc));
}
template <typename T>
VectorXcd FiberConfig<T>::DSolve(const InciCPtrs<T>& inc) {
  return DSolve(inciVect(inc));
}

}  // namespace mss

#endif
