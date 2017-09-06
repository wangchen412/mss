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
  using Inhomo<T>::LocalCS;

 public:
  Assembly(const AssemblyConfig<T>* config, const PosiVect& position = {0, 0},
           double angle = 0, const CS* basis = nullptr)
      : Inhomo<T>(position, ASSEMBLY, angle), config_(config) {
    if (!basis) add_node();
  }
  ~Assembly() {
    for (auto& i : node_) delete i;
  }

  bool Contain(const CS* objCS) const override;

  T Inner(const CS* objCS) const override;
  T Scatter(const CS* objCS) const override;
  T ScatterMode(const CS* objCS, size_t sn) const override;

  // void SetCoeff(const VectorXcd& solution) override { cSc_ = solution; }
  // const VectorXcd& ScatterCoeff() const override { return cSc_; }

  const CSCPtrs& Node() const override { return node_; }
  const CS* Node(size_t i) const override { return node_[i]; }

  // The collocation matrix for assemblies are only available when the nodes
  // are the inner nodes, which are the nodes of the enclosed
  // inhomogeneities. So, temporarily, the collocation matrix of assembly
  // class is disabled. TODO: Add collocation matrix compatability.
  // MatrixXcd ColloMat() const override { return config_->ColloMat(); }

  MatrixXcd TransMat() const override { return config_->TransMat(); }
  size_t NumNode() const override { return config_->NumNode(); }
  size_t NumBv() const override { return config_->NumBv(); }
  size_t NumCoeff() const override { return config_->NumCoeff(); }
  double Width() const override { return config_->Width(); }
  double Height() const override { return config_->Height(); }
  const InhomoCPtrs<T>& inhomo() { return config_->inhomo(); }
  const Inhomo<T>* inhomo(size_t sn) const { return config_->inhomo(sn); }

 private:
  const AssemblyConfig<T>* config_;
  CSCPtrs node_;
  InhomoPtrs<T> inhomo_;

  void add_node();
  void add_inhomo();
};

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T>
bool Assembly<T>::Contain(const CS* objCS) const {
  PosiVect r = objCS->PositionIn(LocalCS());
  return r.x > 0 && r.x < Width() && r.y > 0 && r.y < Height();
}
// template <typename T>
// T Assembly<T>::Scatter(const CS* objCS) const {

// }
template <typename T>
T Assembly<T>::ScatterMode(const CS* objCS, size_t sn) const {
  size_t n = 0;
  while (sn > inhomo(n)->NumCoeff()) sn -= inhomo(n++)->NumCoeff();
  return inhomo(n)->ScatterMode(objCS, sn);
}
template <typename T>
void Assembly<T>::add_node() {
  node_.reserve(config_->NumNode());
  for (auto& i : config_->Node()) node_.push_back(new CS(*i, LocalCS()));
}
template <typename T>
void Assembly<T>::add_inhomo() {
  inhomo_.reserve(config_->inhomo().size());
  // for (auto& i : config_->inhomo()) inhomo_.push_back(new )
}

}  // namespace mss

#endif
