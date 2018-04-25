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

#include "../incident/IncidentPlane.h"
#include "../pre/Input.h"
#include "Boundary.h"
#include "Fiber.h"
#include "FiberConfig.h"
#include "Inhomo.h"

namespace mss {

template <typename T>
class Assembly;

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

  virtual ~AssemblyConfig() { del_inhomo(); }

  size_t NumNode() const { return Node().size(); }
  size_t NumBv() const { return T::NumBv * NumNode(); }
  size_t NumDv() const { return T::NumDv * NumNode(); }
  size_t NumCoeff() const { return num_coeff_; }
  size_t NumCoeff(size_t i) const { return inhomo(i)->NumCoeff(); }
  size_t NumBv_in() const { return num_bv_in_; }

  double KL_m() const { return matrix_->KL(); }
  double KT_m() const { return matrix_->KT(); }
  double Height() const { return height_; }
  double Width() const { return width_; }
  const std::string& ID() const { return ID_; }
  const CSCPtrs& Node() const { return boundary_.Node(); }
  const CS* Node(size_t i) const { return boundary_.Node(i); }
  const std::vector<CSCPtrs>& Edge() const { return boundary_.Edge(); }
  const CSCPtrs& Edge(size_t i) const { return boundary_.Edge(i); }
  CSCPtrs EdgeNode() const;
  const CSCPtrs& Node_in() const { return node_in_; }
  Boundary<T>& Boundary() { return boundary_; }

  // TODO: in-plane
  // Boundary integral matrix. Transfers the effect on the outter boundary to
  // the effect on the inner interfaces.
  MatrixXcd BdIntMatT() const { return boundary_.EffectMatT(inhomoC_); }

  const MatrixXcd& PlaneEDMat() const;
  const MatrixXcd& ColloMat();  // The collocation matrix.
  const MatrixXcd& DcMat();     // The combined matrix for DFT method.
  const MatrixXcd& TransMat();

  // Incident wave effects on the inner interfaces.
  VectorXcd IncVec(const InciCPtrs<T>& incident) const;
  VectorXcd Trans_IncVec(const InciCPtrs<T>& incident) const;
  VectorXcd Trans_IncVec(const VectorXcd& incBv) const;
  // MatrixXcd Trans_BiMat(const MatrixXcd& B) const;

  void Solve(const VectorXcd& incBv, SolveMethod method);
  void CSolve(const VectorXcd& incBv);
  void DSolve(const VectorXcd& incBv);

  void Solve(const InciCPtrs<T>& incident, SolveMethod method);
  void CSolve(const InciCPtrs<T>& incident);
  void DSolve(const InciCPtrs<T>& incident);

  VectorXcd ScatterCoeff() const;
  MatrixXcd ScatterBvMat(const CSCPtrs& objCSs) const;
  MatrixXcd ScatterDvMat(const CSCPtrs& objCSs) const;

  MatrixXcd CylinEDMat(const CSCPtrs& objCSs);
  MatrixXcd CylinEBMat(const CSCPtrs& objCSs);
  MatrixXcd ResBvMat(const CSCPtrs& objCSs);
  MatrixXcd ResDvMat(const CSCPtrs& objCSs);

  const MatrixXcd& InToRstMat();

  const Inhomo<T>* InWhich(const CS* objCS) const;
  T Resultant(
      const CS* objCS, const Inhomo<T>* inhomo,
      const InciCPtrs<T>& incident) const;  // TODO incident maybe not needed.
  T Resultant(const CS* objCS, const InciCPtrs<T>& incident) const;

  T Resultant(const CS* objCS, const VectorXcd& coeff);

  void PrintCoeff(std::ostream& os) const;

  const InhomoCPtrs<T>& inhomo() const { return inhomoC_; }
  const Inhomo<T>* inhomo(size_t sn) const { return inhomoC_[sn]; }

 protected:
  const std::string ID_;
  size_t num_coeff_{0}, num_bv_in_{0};
  InhomoPtrs<T> inhomo_;
  InhomoCPtrs<T> inhomoC_;
  CSCPtrs node_in_;

  FiberConfigPtrs<T> fiber_config_;
  AsmConfigPtrs<T> assembly_config_;

  // const size_t P_;
  const double width_, height_;
  const class Matrix* matrix_;
  mss::Boundary<T> boundary_;
  const input::AssemblyConfig& input_;

  MatrixXcd cc_;
  MatrixXcd dc_;
  MatrixXcd Q_;
  bool cc_computed_{false}, dc_computed_{false};
  bool Q_computed_{false};

  void add_inhomo();
  void add_fiber();
  void add_assembly();
  void del_inhomo();
  void del_fiber();
  void del_assembly();
  MatrixXcd comb_trans_mat() const;  // Combined trans-matrix.
  MatrixXcd inter_identity_mat() const;
  void com_z_mat();
  const Inhomo<T>* nearest(const CS* objCS) const;
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
  dist_solution(ColloMat()
                    .jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV)
                    .solve(incBv));
}

template <typename T>
void AssemblyConfig<T>::DSolve(const VectorXcd& incBv) {
  dist_solution(DcMat().lu().solve(Trans_IncVec(incBv)));
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
VectorXcd AssemblyConfig<T>::ScatterCoeff() const {
  VectorXcd rst(NumCoeff());
  size_t u = 0;
  for (size_t i = 0; i < inhomo_.size(); u += inhomo_[i]->NumCoeff(), i++)
    rst.segment(u, inhomo_[i]->NumCoeff()) = inhomo_[i]->ScatterCoeff();
  return rst;
}

template <typename T>
MatrixXcd AssemblyConfig<T>::ScatterBvMat(const CSCPtrs& objCSs) const {
  MatrixXcd rst(objCSs.size() * T::NumBv, NumCoeff());

  size_t v = 0;
  for (size_t i = 0; i < inhomo_.size(); i++) {
    rst.block(0, v, rst.rows(), inhomo(i)->NumCoeff()) =
        inhomo(i)->ScatterBvMat(objCSs);
    v += inhomo(i)->NumCoeff();
  }

  return rst;
}

template <typename T>
MatrixXcd AssemblyConfig<T>::ScatterDvMat(const CSCPtrs& objCSs) const {
  MatrixXcd rst(objCSs.size() * T::NumDv, NumCoeff());

  size_t v = 0;
  for (size_t i = 0; i < inhomo_.size(); i++) {
    rst.block(0, v, rst.rows(), inhomo(i)->NumCoeff()) =
        inhomo(i)->ScatterDvMat(objCSs);
    v += inhomo(i)->NumCoeff();
  }

  return rst;
}

template <typename T>
const MatrixXcd& AssemblyConfig<T>::PlaneEDMat() const {
  return boundary_.PlaneEDMat(node_in_);
};

template <typename T>
MatrixXcd AssemblyConfig<T>::CylinEBMat(const CSCPtrs& objCSs) {
  // Extrapolate the boundary values at nth edge points with cylindrical
  // waves. The matrix transforms the scattering coefficients to the
  // displacement of incident wave. The boundary points are rearranged for PBC
  // computation.

  MatrixXcd E(objCSs.size() * T::NumBv, NumCoeff());
  E.setZero();

  for (size_t p = 0; p < objCSs.size(); p++) {
    const Inhomo<T>* ni = nearest(objCSs[p]);
    size_t v = 0;
    for (auto& i : inhomo_) {
      if (i == ni) {
        E.block(p * T::NumBv, v, T::NumBv, i->NumCoeff()) =
          i->PsInBvT(objCSs[p]);
        break;
      }
      v += i->NumCoeff();
    }
  }

  // for (size_t p = 0; p < objCSs.size(); p++) {
  //   size_t v = 0;
  //   for (auto& i : inhomo_) {
  //     E.block(p * T::NumBv, v, T::NumBv, i->NumCoeff()) =
  //         i->PsInBvT(objCSs[p]);
  //     v += i->NumCoeff();
  //   }
  // }

  return E * DcMat();
}

template <typename T>
MatrixXcd AssemblyConfig<T>::CylinEDMat(const CSCPtrs& objCSs) {
  // Extrapolate the boundary values at nth edge points with cylindrical
  // waves. The matrix transforms the scattering coefficients to the
  // displacement of incident wave. The boundary points are rearranged for PBC
  // computation.

  MatrixXcd E(objCSs.size() * T::NumDv, NumCoeff());
  E.setZero();

  for (size_t p = 0; p < objCSs.size(); p++) {
    const Inhomo<T>* ni = nearest(objCSs[p]);
    size_t v = 0;
    for (auto& i : inhomo_) {
      if (i == ni)
        E.block(p * T::NumDv, v, T::NumDv, i->NumCoeff()) =
            i->PsInDvT(objCSs[p]);
      v += i->NumCoeff();
    }
  }
  return E * DcMat();
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

// template <typename T>
// MatrixXcd AssemblyConfig<T>::Trans_BiMat(const Eigen::MatrixXcd& B) const {
//   MatrixXcd rst(NumCoeff(), NumBv());
//   for (size_t u = 0; u < inhomo().size(); u++) {
//     size_t i = 0, j = 0;
//     for (size_t k = 0; k < u; k++) {
//       i += inhomo(k)->NumCoeff();
//       j += inhomo(k)->NumBv();
//     }
//     rst.block(i, 0, inhomo(u)->NumCoeff(), NumBv()) =
//         inhomo(u)->TransMat() * B.block(j, 0, inhomo(u)->NumBv(), NumBv());
//   }
//   return rst;
// }

template <typename T>
void AssemblyConfig<T>::dist_solution(const VectorXcd& solution) {
  size_t u = 0;
  for (auto& i : inhomo_) {
    i->SetCoeff(solution.segment(u, i->NumCoeff()));
    u += i->NumCoeff();
  }
}

template <typename T>
const Inhomo<T>* AssemblyConfig<T>::InWhich(const CS* objCS) const {
  const Inhomo<T>* rst = nullptr;
  for (auto& i : inhomo_) {
    rst = i->Contains(objCS);
    if (rst) break;
  }
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
T AssemblyConfig<T>::Resultant(const CS* objCS, const VectorXcd& coeff) {
  dist_solution(coeff);
  T rst(objCS);
  if (InWhich(objCS) == nullptr) {
    for (auto& i : inhomo_) rst += i->Scatter(objCS);
  }
}

template <typename T>
MatrixXcd AssemblyConfig<T>::ResDvMat(const CSCPtrs& objCSs) {
  return ScatterDvMat(objCSs) + CylinEDMat(objCSs);
}

template <typename T>
MatrixXcd AssemblyConfig<T>::ResBvMat(const CSCPtrs& objCSs) {
  // Transformation from scattering coefficients to resultant boundary values.

  return ScatterBvMat(objCSs) + CylinEBMat(objCSs);
}

template <typename T>
void AssemblyConfig<T>::PrintCoeff(std::ostream& os) const {
  for (auto& i : inhomo_) i->PrintCoeff(os);
}

template <typename T>
void AssemblyConfig<T>::add_inhomo() {
  add_fiber();
  add_assembly();
  for (auto& i : inhomo_) {
    num_coeff_ += i->NumCoeff();
    num_bv_in_ += i->NumBv();
    node_in_.insert(node_in_.end(), i->Node().begin(), i->Node().end());
    inhomoC_.push_back(i);
  }
}

template <typename T>
void AssemblyConfig<T>::add_fiber() {
  for (auto& i : input_.fiber) {
    FiberConfig<T>* cp = FindPtrID(fiber_config_, i.configID);
    if (cp == nullptr) {
      cp = new FiberConfig<T>(*(i.config), matrix_);
      fiber_config_.push_back(cp);
    }
    inhomo_.push_back(new Fiber<T>(cp, i.position));
  }
}

template <typename T>
void AssemblyConfig<T>::add_assembly() {
  for (auto& i : input_.assembly) {
    AssemblyConfig<T>* cp = FindPtrID(assembly_config_, i.configID);
    if (cp == nullptr) {
      cp = new AssemblyConfig<T>(*(i.config), matrix_);
      assembly_config_.push_back(cp);
    }
    inhomo_.push_back(new Assembly<T>(cp, i.position, i.angle));
  }
}

template <typename T>
void AssemblyConfig<T>::del_inhomo() {
  for (auto& i : inhomo_) delete i;
  for (auto& i : fiber_config_) delete i;
  for (auto& i : assembly_config_) delete i;
}

template <typename T>
const Inhomo<T>* AssemblyConfig<T>::nearest(const CS* objCS) const {
  // Find the nearest inhomogeneity.

  const Inhomo<T>* rst = inhomo_[0];
  double d = inhomo_[0]->LocalCS()->Distance(objCS);
  for (auto& i : inhomo_) {
    double di = i->LocalCS()->Distance(objCS);
    if (di < d) {
      d = di;
      rst = i;
    }
  }

  return rst;
}

template <typename T>
CSCPtrs AssemblyConfig<T>::EdgeNode() const {
  CSCPtrs rst;
  for (const CSCPtrs& i : Edge()) rst.insert(rst.end(), i.begin(), i.end());
  return rst;
}

}  // namespace mss

#endif
