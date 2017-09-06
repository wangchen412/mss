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

#ifndef MSS_ASSEMBLYCONFIG_H
#define MSS_ASSEMBLYCONFIG_H

#include "../pre/Input.h"
#include "Boundary.h"
#include "Fiber.h"
#include "FiberConfig.h"
#include "Inhomo.h"

namespace mss {

template <typename T>
class AssemblyConfig;
template <typename T>
using AsmConfigPtrs = std::vector<AssemblyConfig<T>*>;
template <typename T>
using AsmConfigCPtrs = std::vector<const AssemblyConfig<T>*>;

template <typename T>
class AssemblyConfig {
 public:
  AssemblyConfig(const input::AssemblyConfig& input, const Matrix* matrix)
      : ID_(input.ID),
        width_(input.width),
        height_(input.height),
        matrix_(matrix),
        boundary_(input.pointDensity, height_, width_, matrix),
        input_(input) {
    add_inhomo();
  }

  virtual ~AssemblyConfig() { delete_inhomo(); }

  double CharLength() const { return height_ + width_; }
  size_t NumNode() const { return Node().size(); }
  size_t NumBv() const { return T::NumBv * NumNode(); }
  size_t NumCoeff() const { return num_coeff_; }
  size_t NumBv_in() const { return num_bv_in_; }

  double Height() const { return height_; }
  double Width() const { return width_; }
  const std::string& ID() const { return ID_; }
  const CSCPtrs& Node() const { return boundary_.Node(); }

  // TODO: in-plane
  MatrixXcd BdIntMatT() const { return boundary_.EffectMatT(inhomoC_); }

  const MatrixXcd& ColloMat();
  const MatrixXcd& DcMat();
  const MatrixXcd& TransMat();
  MatrixXcd GramMat();

  VectorXcd IncVec(const InciCPtrs<T>& incident) const;
  VectorXcd Trans_IncVec(const InciCPtrs<T>& incident) const;
  VectorXcd Trans_IncVec(const VectorXcd& incBv) const;
  MatrixXcd Trans_BiMat(const MatrixXcd& B) const;

  void Solve(const VectorXcd& incBv, SolveMethod method);
  void CSolve(const VectorXcd& incBv);
  void DSolve(const VectorXcd& incBv);

  void Solve(const InciCPtrs<T>& incident, SolveMethod method);
  void CSolve(const InciCPtrs<T>& incident);
  void DSolve(const InciCPtrs<T>& incident);

  Inhomo<T>* InWhich(const CS* objCS) const;
  T Resultant(const CS* objCS, const Inhomo<T>* inhomo,
              const InciCPtrs<T>& incident) const;
  T Resultant(const CS* objCS, const InciCPtrs<T>& incident) const;

  void PrintCoeff(std::ostream& os) const;

  const InhomoCPtrs<T>& inhomo() const { return inhomoC_; }
  const Inhomo<T>* inhomo(size_t sn) const { return inhomoC_[sn]; }

 protected:
  const std::string ID_;
  size_t num_coeff_{0}, num_bv_in_{0};
  InhomoPtrs<T> inhomo_;
  InhomoCPtrs<T> inhomoC_;
  FiberConfigPtrs<T> fiber_config_;
  AsmConfigPtrs<T> assembly_config_;

  // const size_t P_;
  const double width_, height_;

  const class Matrix* matrix_;
  Boundary<T> boundary_;
  const input::AssemblyConfig& input_;

  MatrixXcd cc_;
  MatrixXcd dc_;
  MatrixXcd Q_;
  bool cc_computed_{false}, dc_computed_{false}, Q_computed_{false};

  void add_inhomo();
  void add_fiber();
  void add_fiber_config();
  void delete_inhomo();
  void delete_fiber_config();
  MatrixXcd comb_trans_mat() const;  // Combined trans-matrix.
  void dist_solution(const VectorXcd& solution);
};

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T>
void AssemblyConfig<T>::Solve(const VectorXcd& incBv, SolveMethod method) {
  switch (method) {
    case COLLOCATION:
      CSolve(incBv);
      break;
    case DFT:
      DSolve(incBv);
      break;
    default:
      exit_error_msg({"Unknown method."});
  }
}

template <typename T>
void AssemblyConfig<T>::CSolve(const VectorXcd& incBv) {
  // Jacobi SVD:
  auto svd = ColloMat().jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
  VectorXcd solution = svd.solve(incBv);
  dist_solution(solution);
}

template <typename T>
void AssemblyConfig<T>::DSolve(const VectorXcd& incBv) {
  VectorXcd solution = DcMat().lu().solve(Trans_IncVec(incBv));
  dist_solution(solution);
}

template <typename T>
void AssemblyConfig<T>::Solve(const InciCPtrs<T>& incident,
                              SolveMethod method) {
  switch (method) {
    case COLLOCATION:
      CSolve(incident);
      break;
    case DFT:
      DSolve(incident);
      break;
    default:
      exit_error_msg({"Unknown method."});
  }
}

template <typename T>
void AssemblyConfig<T>::CSolve(const InciCPtrs<T>& incident) {
  CSolve(IncVec(incident));
}

template <typename T>
void AssemblyConfig<T>::DSolve(const InciCPtrs<T>& incident) {
  VectorXcd solution = DcMat().lu().solve(Trans_IncVec(incident));
  dist_solution(solution);
}

template <typename T>
const MatrixXcd& AssemblyConfig<T>::ColloMat() {
  if (cc_computed_) return cc_;

  cc_.resize(NumBv_in(), NumCoeff());

#ifdef NDEBUG
#pragma omp parallel for
#endif
  for (size_t v = 0; v < inhomo_.size(); v++) {
    int i = 0, j = 0, Nv = inhomo_[v]->NumCoeff();
    for (size_t k = 0; k < v; k++) j += inhomo_[k]->NumCoeff();
    for (size_t u = 0; u < inhomo_.size(); u++) {
      int Nu = inhomo_[u]->NumBv();

      cc_.block(i, j, Nu, Nv) =
          u == v ? inhomo_[v]->ColloMat() : -inhomo_[v]->ModeMat(inhomo_[u]);
      i += Nu;
    }
  }

  cc_computed_ = true;

  return cc_;
}

template <typename T>
const MatrixXcd& AssemblyConfig<T>::TransMat() {
  if (Q_computed_) return Q_;
  Q_ = DcMat().inverse() * comb_trans_mat() * BdIntMatT();  // TODO: in-plane
  Q_computed_ = true;
  return Q_;
}

template <typename T>
const MatrixXcd& AssemblyConfig<T>::DcMat() {
  if (dc_computed_) return dc_;

  dc_.resize(NumCoeff(), NumCoeff());
#ifdef NDEBUG
#pragma omp parallel for
#endif
  for (size_t v = 0; v < inhomo_.size(); v++) {
    int i = 0, j = 0, Nv = inhomo_[v]->NumCoeff();
    MatrixXcd I(Nv, Nv);
    I.setIdentity();
    for (size_t k = 0; k < v; k++) j += inhomo_[k]->NumCoeff();
    for (size_t u = 0; u < inhomo_.size(); u++) {
      int Nu = inhomo_[u]->NumCoeff();
      dc_.block(i, j, Nu, Nv) =
          u == v ? I
                 : -inhomo_[u]->TransMat() * inhomo_[v]->ModeMat(inhomo_[u]);
      i += Nu;
    }
  }
  dc_computed_ = true;

  return dc_;
}

template <typename T>
MatrixXcd AssemblyConfig<T>::GramMat() {
  return ColloMat().transpose() * ColloMat();
}

template <typename T>
MatrixXcd AssemblyConfig<T>::comb_trans_mat() const {
  MatrixXcd rst(NumCoeff(), NumBv_in());
  rst.setZero();
  size_t u = 0, v = 0;
  for (auto& i : inhomo_) {
    rst.block(u, v, i->NumCoeff(), i->NumBv()) = i->TransMat();
    u += i->NumCoeff();
    v += i->NumBv();
  }
  return rst;
}

template <typename T>
VectorXcd AssemblyConfig<T>::IncVec(const InciCPtrs<T>& incident) const {
  // The effect vector of incident wave along all the interfaces inside the
  // assembly.

  VectorXcd rst(NumBv_in());
  size_t u = 0;
  for (auto& i : inhomo_) {
    rst.segment(u, i->NumBv()) = i->IncVec(incident);
    u += i->NumBv();
  }
  return rst;
}

template <typename T>
VectorXcd AssemblyConfig<T>::Trans_IncVec(
    const InciCPtrs<T>& incident) const {
  VectorXcd rst(NumCoeff());
  size_t u = 0;
  for (auto& i : inhomo_) {
    rst.segment(u, i->NumCoeff()) = i->TransIncVec(incident);
    u += i->NumCoeff();
  }
  return rst;
}

template <typename T>
VectorXcd AssemblyConfig<T>::Trans_IncVec(const VectorXcd& incBv) const {
  VectorXcd rst(NumCoeff());
  size_t u = 0, v = 0;
  for (auto& i : inhomo_) {
    rst.segment(u, i->NumCoeff()) =
        i->TransMat() * incBv.segment(v, i->NumBv());
    u += i->NumCoeff();
    v += i->NumBv();
  }
  return rst;
}

template <typename T>
MatrixXcd AssemblyConfig<T>::Trans_BiMat(const Eigen::MatrixXcd& B) const {
  MatrixXcd rst(NumCoeff(), NumBv());
  for (size_t u = 0; u < inhomo().size(); u++) {
    size_t i = 0, j = 0;
    for (size_t k = 0; k < u; k++) {
      i += inhomo(k)->NumCoeff();
      j += inhomo(k)->NumBv();
    }
    rst.block(i, 0, inhomo(u)->NumCoeff(), NumBv()) =
        inhomo(u)->TransMat() * B.block(j, 0, inhomo(u)->NumBv(), NumBv());
  }
  return rst;
}

template <typename T>
void AssemblyConfig<T>::dist_solution(const VectorXcd& solution) {
  size_t u = 0;
  for (auto& i : inhomo_) {
    i->SetCoeff(solution.segment(u, i->NumCoeff()));
    u += i->NumCoeff();
  }
}

template <typename T>
Inhomo<T>* AssemblyConfig<T>::InWhich(const CS* objCS) const {
  // Return the pointer to the inhomogeneity in which the objCS is.
  // The local CS is considered as the global CS.

  Inhomo<T>* rst = nullptr;
  for (auto& i : inhomo_)
    if (i->Contain(objCS)) rst = i;
  return rst;
}

template <typename T>
T AssemblyConfig<T>::Resultant(const CS* objCS, const Inhomo<T>* in,
                               const InciCPtrs<T>& incident) const {
  T rst(objCS);
  if (in)
    rst = in->Inner(objCS);
  else {
    for (auto& i : incident) rst += i->Effect(objCS);
    for (auto& i : inhomo_) rst += i->Scatter(objCS);
  }
  return rst;
}

template <typename T>
T AssemblyConfig<T>::Resultant(const CS* objCS,
                               const InciCPtrs<T>& incident) const {
  return Resultant(objCS, InWhich(objCS), incident);
}

template <typename T>
void AssemblyConfig<T>::PrintCoeff(std::ostream& os) const {
  for (auto& i : inhomo_) i->PrintCoeff(os);
}

template <typename T>
void AssemblyConfig<T>::add_inhomo() {
  add_fiber();
  for (auto& i : inhomo_) {
    num_coeff_ += i->NumCoeff();
    num_bv_in_ += i->NumBv();
    inhomoC_.push_back(i);
  }
}

template <typename T>
void AssemblyConfig<T>::add_fiber() {
  add_fiber_config();
  for (auto& i : input_.fiber)
    inhomo_.push_back(
        new Fiber<T>(FindPtrID(fiber_config_, i.configID), i.position));
}

template <typename T>
void AssemblyConfig<T>::add_fiber_config() {
  for (auto& i : *input_.fiber_config)
    fiber_config_.push_back(new FiberConfig<T>(i, matrix_));
}

template <typename T>
void AssemblyConfig<T>::delete_inhomo() {
  for (auto& i : inhomo_) delete i;
  delete_fiber_config();
}

template <typename T>
void AssemblyConfig<T>::delete_fiber_config() {
  for (auto& i : fiber_config_) delete i;
}

}  // namespace mss

#endif
