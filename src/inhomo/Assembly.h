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

#ifndef MSS_ASSEMBLY_H
#define MSS_ASSEMBLY_H

#include "AssemblyConfig.h"
#include "Inhomo.h"

namespace mss {

template <typename T>
class Assembly : public Inhomo<T> {
 public:
  Assembly(AssemblyConfig<T>* config, const PosiVect& position = {0, 0},
           double angle = 0, const CS* basis = nullptr)
      : Inhomo<T>(position, ASSEMBLY, angle, basis), config_(config) {
    if (!basis) add_node();
    add_inhomo();
  }
  Assembly(const Assembly* p, const CS* basis)
      : Inhomo<T>(p->Position(), ASSEMBLY, p->Angle(), basis),
        config_(p->config_) {
    add_inhomo();
  }
  ~Assembly() {
    for (auto& i : node_) delete i;
  }

  const Inhomo<T>* Contains(const CS* objCS) const override;

  T Inner(const CS* objCS) const override { return T(objCS); }
  T Scatter(const CS* objCS) const override;
  T ScatterMode(const CS* objCS, size_t sn) const override;

  // TODO
  T PsInMode(const CS*, size_t) const override { return T(); }
  MatrixXcd PsInBvT(const CS*) const override { return MatrixXcd(); }
  MatrixXcd PsInDvT(const CS*) const override { return MatrixXcd(); }
  MatrixXcd PsInBvMatT(const CSCPtrs&) const override { return MatrixXcd(); }
  MatrixXcd PsInDvMatT(const CSCPtrs&) const override { return MatrixXcd(); }

  // void SetCoeff(const VectorXcd& solution) override { cSc_ = solution; }
  // const VectorXcd& ScatterCoeff() const override { return cSc_; }

  const CSCPtrs& Node() const override { return node_; }
  const CS* Node(size_t i) const override { return node_[i]; }

  // The collocation matrix for assemblies are only available when the nodes
  // are the inner nodes, which are the nodes of the enclosed
  // inhomogeneities. So, temporarily, the collocation matrix of assembly
  // class is disabled. TODO: Add collocation matrix compatability.
  // MatrixXcd ColloMat() const override { return config_->ColloMat(); }

  const MatrixXcd& TransMat() const override { return config_->TransMat(); }
  size_t NumNode() const override { return config_->NumNode(); }
  size_t NumBv() const override { return config_->NumBv(); }
  size_t NumCoeff() const override { return config_->NumCoeff(); }
  double Width() const { return config_->Width(); }
  double Height() const { return config_->Height(); }
  const InhomoCPtrs<T>& inhomo() const { return inhomoC_; }
  const Inhomo<T>* inhomo(size_t sn) const { return inhomoC_[sn]; }

  VectorXcd DSolve(const InciCPtrs<T>& incident) const;

  void SetCoeff(const VectorXcd& solution) override;
  const VectorXcd& ScatterCoeff() const override { return cSc_; }

  using Inhomo<T>::LocalCS;
  using Inhomo<T>::IncVec;

 private:
  AssemblyConfig<T>* config_;
  CSCPtrs node_;
  InhomoPtrs<T> inhomo_;
  InhomoCPtrs<T> inhomoC_;
  VectorXcd cSc_;

  void add_node();
  void add_inhomo();
  void add_fiber(const Inhomo<T>* p);
  void add_assembly(const Inhomo<T>* p);
};

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T>
VectorXcd Assembly<T>::DSolve(const InciCPtrs<T>& incident) const {
  return config_->TransMat() * IncVec(incident);
}
template <typename T>
const Inhomo<T>* Assembly<T>::Contains(const CS* objCS) const {
  const Inhomo<T>* rst = nullptr;
  for (auto& i : inhomo()) {
    rst = i->Contains(objCS);
    if (rst) break;
  }
  return rst;
}
template <typename T>
T Assembly<T>::Scatter(const CS* objCS) const {
  T rst(objCS);
  for (auto& i : inhomo_) rst += i->Scatter(objCS);
  return rst;
}
template <typename T>
T Assembly<T>::ScatterMode(const CS* objCS, size_t sn) const {
  size_t n = 0;
  while (sn > inhomo(n)->NumCoeff()) sn -= inhomo(n++)->NumCoeff();
  return inhomo(n)->ScatterMode(objCS, sn);
}
template <typename T>
void Assembly<T>::SetCoeff(const VectorXcd& solution) {
  size_t j = 0;
  for (auto& i : inhomo_) {
    i->SetCoeff(solution.segment(j, i->NumCoeff()));
    j += i->NumCoeff();
  }
  cSc_ = solution;
}
template <typename T>
void Assembly<T>::add_node() {
  node_.reserve(config_->NumNode());
  for (auto& i : config_->Node()) node_.push_back(new CS(*i, LocalCS()));
}
template <typename T>
void Assembly<T>::add_inhomo() {
  inhomo_.reserve(config_->inhomo().size());
  for (auto& i : config_->inhomo()) {
    switch (i->Type()) {
      case FIBER:
        add_fiber(i);
        break;
      case ASSEMBLY:
        add_assembly(i);
        break;
      default:
        exit_error_msg({"Unknown type."});
    }
  }
  for (auto& i : inhomo_) inhomoC_.push_back(i);
}
template <typename T>
void Assembly<T>::add_fiber(const Inhomo<T>* p) {
  const Fiber<T>* fp = dynamic_cast<const Fiber<T>*>(p);
  inhomo_.push_back(new Fiber<T>(fp, LocalCS()));
}
template <typename T>
void Assembly<T>::add_assembly(const Inhomo<T>* p) {
  const Assembly<T>* ap = dynamic_cast<const Assembly<T>*>(p);
  inhomo_.push_back(new Assembly<T>(ap, LocalCS()));
}

}  // namespace mss

#endif
