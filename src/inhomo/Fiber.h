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

#ifndef MSS_FIBER_H
#define MSS_FIBER_H

#include "FiberConfig.h"
#include "Inhomo.h"

namespace mss {

template <typename T>
class Fiber : public Inhomo<T> {
public:
  Fiber(FiberConfig<T>* config, const PosiVect& position = {0, 0},
        const CS* basis = nullptr)
      : Inhomo<T>(position, FIBER, 0, basis),
        config_(config),
        cSc_(NumCoeff()),
        cIn_(NumCoeff()) {
    if (!basis) add_node();
  }
  Fiber(const Fiber* p, const CS* basis)
      : Inhomo<T>(p->Position(), FIBER, 0, basis),
        config_(p->config_),
        cSc_(p->NumCoeff()),
        cIn_(p->NumCoeff()) {}
  ~Fiber() {
    for (auto& i : node_) delete i;
  }

  const MatrixXcd& ColloMat() const override { return config_->ColloMat(); }
  const MatrixXcd& TransMat() const override { return config_->TransMat(); }
  size_t NumNode() const override { return config_->NumNode(); }
  size_t NumBv() const override { return config_->NumBv(); }
  size_t NumCoeff() const override { return config_->NumCoeff(); }
  double Radius() const { return config_->Radius(); }

  void SetCoeff(const VectorXcd& solution) override;

  // Check if the position of the objective CS is inside the fiber.
  bool Contain(const CS* objCS) const override;

  T Scatter(const CS* objCS) const override;
  T Inner(const CS* objCS) const override;

  // The nth Modes. The n should be the serial number, instead of the order.
  // The transformation from the serial number to the order should be done in
  // the derived class.
  // n starts at zero.
  T ScatterMode(const CS* objCS, size_t sn) const override;
  T InnerMode(const CS* objCS, size_t sn) const;

  VectorXcd Solve(const VectorXcd& incBv, SolveMethod method) const;
  VectorXcd CSolve(const VectorXcd& incBv) const;
  VectorXcd DSolve(const VectorXcd& incBv) const;

  VectorXcd Solve(const InciCPtrs<T>& incident, SolveMethod method) const;
  VectorXcd CSolve(const InciCPtrs<T>& incident) const;
  VectorXcd DSolve(const InciCPtrs<T>& incident) const;

  FiberConfig<T>* Config() const { return config_; }
  const CSCPtrs& Node() const override { return node_; }
  const CS* Node(size_t i) const override { return node_[i]; }
  const VectorXcd& ScatterCoeff() const override { return cSc_; }

  using Inhomo<T>::LocalCS;
  using Inhomo<T>::IncVec;

 private:
  // static FiberConfigPtrs<T> fiber_config_;
  FiberConfig<T>* config_;
  CSCPtrs node_;
  VectorXcd cSc_, cIn_;

  T scatterModeL(const CS* objCS, int n) const;
  T scatterModeT(const CS* objCS, int n) const;
  T innerModeL(const CS* objCS, int n) const;
  T innerModeT(const CS* objCS, int n) const;

  void add_node();
  int od(size_t sn) const { return sn - config_->TopOrder(); }
};

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T>
bool Fiber<T>::Contain(const CS* objCS) const {
  return objCS->PositionIn(LocalCS()).Length() < config_->Radius();
}
template <typename T>
T Fiber<T>::Scatter(const CS* objCS) const {
  T rst(objCS);
  for (size_t i = 0; i < NumCoeff(); i++)
    rst += ScatterMode(objCS, i) * cSc_(i);
  return rst;
}
template <typename T>
T Fiber<T>::Inner(const CS* objCS) const {
  T rst(objCS);
  for (size_t i = 0; i < NumCoeff(); i++)
    rst += InnerMode(objCS, i) * cIn_(i);
  return rst;
}
template <typename T>
void Fiber<T>::SetCoeff(const VectorXcd& solution) {
  assert(solution.size() == long(NumCoeff()));
  for (long i = 0; i < solution.size(); i++) {
    cSc_(i) = solution(i);
    cIn_(i) = cSc_(i) * config_->TT(od(i));  // TODO: in-plane problem.
  }
}
template <>
inline StateIP Fiber<IP>::ScatterMode(const CS* objCS, size_t sn) const {
  if (sn <= NumCoeff() / 2)
    return scatterModeL(objCS, od(sn));
  else
    return scatterModeT(objCS, od(sn));
}
template <>
inline StateAP Fiber<AP>::ScatterMode(const CS* objCS, size_t sn) const {
  return scatterModeT(objCS, od(sn));
}
template <>
inline StateIP Fiber<IP>::InnerMode(const CS* objCS, size_t sn) const {
  if (sn <= NumCoeff() / 2)
    return innerModeL(objCS, od(sn));
  else
    return innerModeT(objCS, od(sn));
}
template <>
inline StateAP Fiber<AP>::InnerMode(const CS* objCS, size_t sn) const {
  return innerModeT(objCS, od(sn));
}
template <typename T>
T Fiber<T>::scatterModeL(const CS* objCS, int n) const {
  return ModeL<T>(LocalCS(), objCS,
                  EigenFunctor(Hn, n, config_->KL_m(), Radius()),
                  config_->Material_m());
}
template <typename T>
T Fiber<T>::innerModeL(const CS* objCS, int n) const {
  return ModeL<T>(LocalCS(), objCS,
                  EigenFunctor(Jn, n, config_->KL(), Radius()),
                  config_->Material());
}
template <typename T>
T Fiber<T>::scatterModeT(const CS* objCS, int n) const {
  return ModeT<T>(LocalCS(), objCS,
                  EigenFunctor(Hn, n, config_->KT_m(), Radius()),
                  config_->Material_m());
}
template <typename T>
T Fiber<T>::innerModeT(const CS* objCS, int n) const {
  return ModeT<T>(LocalCS(), objCS,
                  EigenFunctor(Jn, n, config_->KT(), Radius()),
                  config_->Material());
}
template <typename T>
void Fiber<T>::add_node() {
  node_.reserve(config_->NumNode());
  for (auto& i : config_->Node()) node_.push_back(new CS(*i, LocalCS()));
}
template <typename T>
VectorXcd Fiber<T>::Solve(const VectorXcd& incBv, SolveMethod method) const {
  return config_->Solve(incBv, method);
}
template <typename T>
VectorXcd Fiber<T>::CSolve(const VectorXcd& incBv) const {
  return config_->CSolve(incBv);
}
template <typename T>
VectorXcd Fiber<T>::DSolve(const VectorXcd& incBv) const {
  return config_->DSolve(incBv);
}
template <typename T>
VectorXcd Fiber<T>::Solve(const InciCPtrs<T>& inc, SolveMethod method) const {
  return config_->Solve(IncVec(inc), method);
}
template <typename T>
VectorXcd Fiber<T>::CSolve(const InciCPtrs<T>& inc) const {
  return config_->CSolve(IncVec(inc));
}
template <typename T>
VectorXcd Fiber<T>::DSolve(const InciCPtrs<T>& inc) const {
  return config_->DSolve(IncVec(inc));
}

}  // namespace mss

#endif
