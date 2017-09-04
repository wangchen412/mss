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
           double angle = 0)
      : Inhomo<T>(position, ASSEMBLY, angle), config_(config) {
    add_node();
  }

  virtual ~Assembly() { del_node(); }

  // TODO: The interactions among assemblies
  bool Contain(const CS* objCS) const override;
  T Scatter(const CS* objCS) const override;
  T Inner(const CS* objCS) const override;
  T ScatterMode(const CS* objCS, size_t sn) const override;
  T InnerMode(const CS* objCS, size_t sn) const override;
  T ScatterModeL(const CS* objCS, int n) const;
  T ScatterModeT(const CS* objCS, int n) const;
  T InnerModeL(const CS* local, int n) const;
  T InnerModeT(const CS* local, int n) const;

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

  VectorXcd Solve(const VectorXcd& incBv, SolveMethod method) const override;
  VectorXcd CSolve(const VectorXcd& incBv) const override;
  VectorXcd DSolve(const VectorXcd& incBv) const override;
  VectorXcd IncVec(const InciCPtrs<T>& inc) const override;
  VectorXcd Solve(const InciCPtrs<T>& inc, SolveMethod method) const override;
  VectorXcd CSolve(const InciCPtrs<T>& inc) const override;
  VectorXcd DSolve(const InciCPtrs<T>& inc) const override;

 private:
  const AssemblyConfig<T>* config_;
  CSCPtrs node_;

  void add_node();
  void del_node();
};

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T>
void Assembly<T>::add_node() {
  node_.reserve(config_->NumNode());
  for (auto& i : config_->Node()) node_.push_back(new CS(*i, LocalCS()));
}
template <typename T>
void Assembly<T>::del_node() {
  for (auto& i : node_) delete i;
}

}  // namespace mss

#endif
