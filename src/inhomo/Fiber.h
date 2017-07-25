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

#include "ConfigFiber.h"
#include "Inhomogeneity.h"

namespace mss {

template <typename T>
class Fiber : public Inhomogeneity<T> {
 public:
  Fiber(const ConfigFiber<T>* config, const PosiVect& position = 0)
      : Inhomogeneity<T>(position),
        config_(config),
        cSc_(NoC()),
        cIn_(NoC()) {
    add_node();
  }

  virtual ~Fiber() { delete_node(); }

  size_t NoN() const override { return config_->NoN(); }
  size_t NoE() const override { return config_->NoE(); }
  size_t NoC() const override { return config_->NoC(); }
  const double& Radius() { return config_->Radius(); }

  void SetCoeff(const Eigen::VectorXcd& solution) override;
  void PrintCoeff(std::ostream& os) const override;

  // Check if the position of the objective CS is inside the fiber.
  bool Contain(const CS* objCS) const override;

  // Resultant states.
  T Scatter(const CS* objCS) const override;
  T Inner(const CS* objCS) const override;

  // The nth Modes. The n should be the serial number, instead of the order.
  // The transformation from the serial number to the order should be done in
  // the derived class.
  // n starts at zero.
  T ScatterMode(const CS* objCS, const size_t& sn) const override;
  T InnerMode(const CS* objCS, const size_t& sn) const override;

  // Return the effect of modes of other source inhomogeneity at the
  // collocation points.
  Eigen::MatrixXcd ModeMatrix(const Inhomogeneity<T>* source) const override;

  Eigen::VectorXcd InciVect(const InciPtrs<T>& incident) const override;
  Eigen::VectorXcd Solve(const InciPtrs<T>& incident) const override;

  const ConfigFiber<T>* Config() const { return config_; }
  const std::vector<CS*>& Node() const override { return node_; }
  const Eigen::VectorXcd& ScatterCoeff() const override { return cSc_; }

 private:
  const ConfigFiber<T>* config_;
  std::vector<CS*> node_;
  Eigen::VectorXcd cSc_, cIn_;

  T scatterModeL(const CS* objCS, int n) const;
  T scatterModeT(const CS* objCS, int n) const;
  T innerModeL(const CS* objCS, int n) const;
  T innerModeT(const CS* objCS, int n) const;

  void add_node();
  void delete_node();
  int od(const size_t& sn) const { return sn - config_->TopOrder(); }
};  // namespace mss

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T>
inline bool Fiber<T>::Contain(const CS* objCS) const {
  return objCS->PositionIn(this->LocalCS()).Length() < config_->Radius();
}
template <typename T>
inline T Fiber<T>::Scatter(const CS* objCS) const {
  T rst(objCS);
  for (size_t i = 0; i < NoC(); i++) rst += ScatterMode(objCS, i) * cSc_(i);
  return rst;
}
template <typename T>
inline T Fiber<T>::Inner(const CS* objCS) const {
  T rst(objCS);
  for (size_t i = 0; i < NoC(); i++) rst += InnerMode(objCS, i) * cIn_(i);
  return rst;
}
template <typename T>
inline void Fiber<T>::SetCoeff(const Eigen::VectorXcd& solution) {
  assert(solution.size() == long(NoC()));
  for (long i = 0; i < solution.size(); i++) {
    cSc_(i) = solution(i);
    cIn_(i) = cSc_(i) * config_->TT(od(i));  // TODO: in-plane problem.
  }
}
template <typename T>
inline void Fiber<T>::PrintCoeff(std::ostream& os) const {
  os << setMaxPrecision << cSc_ << std::endl;
}
template <>
inline StateIP Fiber<StateIP>::ScatterMode(const CS* objCS,
                                           const size_t& sn) const {
  if (sn <= NoC() / 2)
    return scatterModeL(objCS, od(sn));
  else
    return scatterModeT(objCS, od(sn));
}
template <>
inline StateAP Fiber<StateAP>::ScatterMode(const CS* objCS,
                                           const size_t& sn) const {
  return scatterModeT(objCS, od(sn));
}
template <>
inline StateIP Fiber<StateIP>::InnerMode(const CS* objCS,
                                         const size_t& sn) const {
  if (sn <= NoC() / 2)
    return innerModeL(objCS, od(sn));
  else
    return innerModeT(objCS, od(sn));
}
template <>
inline StateAP Fiber<StateAP>::InnerMode(const CS* objCS,
                                         const size_t& sn) const {
  return innerModeT(objCS, od(sn));
}
template <typename T>
inline T Fiber<T>::scatterModeL(const CS* objCS, int n) const {
  return config_->ModeL(this->LocalCS(), objCS,
                        EigenFunctor(Hn, n, config_->KL_m()),
                        config_->Material_m());
}
template <typename T>
inline T Fiber<T>::innerModeL(const CS* objCS, int n) const {
  return config_->ModeL(this->LocalCS(), objCS,
                        EigenFunctor(Jn, n, config_->KL()),
                        config_->Material());
}
template <typename T>
inline T Fiber<T>::scatterModeT(const CS* objCS, int n) const {
  return config_->ModeT(this->LocalCS(), objCS,
                        EigenFunctor(Hn, n, config_->KT_m()),
                        config_->Material_m());
}
template <typename T>
inline T Fiber<T>::innerModeT(const CS* objCS, int n) const {
  return config_->ModeT(this->LocalCS(), objCS,
                        EigenFunctor(Jn, n, config_->KT()),
                        config_->Material());
}
template <typename T>
inline void Fiber<T>::add_node() {
  node_.reserve(config_->NoN());
  for (auto& i : config_->Node())
    node_.push_back(new CS(*i, this->LocalCS()));
}
template <typename T>
inline void Fiber<T>::delete_node() {
  for (auto& i : node_) delete i;
}
template <typename T>
inline Eigen::MatrixXcd Fiber<T>::ModeMatrix(
    const Inhomogeneity<T>* source) const {
  if (source == this) return config_->TransMatrix();

  Eigen::MatrixXcd M(this->NoE(), source->NoC());
  for (size_t sn = 0; sn < source->NoC(); sn++)
    for (size_t i = 0; i < this->NoN(); i++)
      M.block<T::NoBV, 1>(T::NoBV * i, sn) =
          source->ScatterMode(node_[i], sn).BV();
  return M *= -1;
}
template <typename T>
inline Eigen::VectorXcd Fiber<T>::InciVect(
    const InciPtrs<T>& incident) const {
  Eigen::VectorXcd rst(NoE());
  rst.setZero();
  for (auto& i : incident) rst += i->EffectBV(Node());
  return rst;
}
template <typename T>
inline Eigen::VectorXcd Fiber<T>::Solve(const InciPtrs<T>& incident) const {
  return config_->Solve(InciVect(incident));
}

}  // namespace mss

#endif
