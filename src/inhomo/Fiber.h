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
  explicit Fiber(const ConfigFiber<T>* config, const PosiVect& position)
      : Inhomogeneity<T>(position), config_(config) {
    add_node();
  }

  // The nth Modes. The n should be the serial number, instead of the order.
  // The transformation from the serial number to the order should be done in
  // the derived class.
  // n starts at zero.
  T ScatterModeL(const CS* objCS, const size_t& n) const override;
  T ScatterModeT(const CS* objCS, const size_t& n) const override;
  T InnerModeL(const CS* objCS, const size_t& n) const override;
  T InnerModeT(const CS* objCS, const size_t& n) const override;

  Eigen::VectorXcd InVector(
      const std::vector<Incident<T>*>& incident) const override;

  // Return the effect of modes of other source inhomogeneity at the
  // collocation points.
  Eigen::MatrixXcd ModeMatrix(const Inhomogeneity<T>* source) const override;

 protected:
  const ConfigFiber<T>* config_;
  std::vector<CS> node_;

  // Transform serial number n to order.
  int od(const size_t& n) const { return n - config_->TopOrder(); }
  void add_node();
};  // namespace mss

// ---------------------------------------------------------------------------
// Inline functions:

template <>
inline StateIP Fiber<StateIP>::ScatterModeL(const CS* objCS,
                                            const size_t& n) const {
  return config_->ModeL(LocalCS(), objCS,
                        EigenFunctor(Hn, od(n), config_->KL_m()),
                        config_->Material_m());
}
template <>
inline StateIP Fiber<StateIP>::ScatterModeT(const CS* objCS,
                                            const size_t& n) const {
  return config_->ModeT(LocalCS(), objCS,
                        EigenFunctor(Hn, od(n), config_->KT_m()),
                        config_->Material_m());
}
template <>
inline StateIP Fiber<StateIP>::InnerModeL(const CS* objCS,
                                          const size_t& n) const {
  return config_->ModeL(LocalCS(), objCS,
                        EigenFunctor(Jn, od(n), config_->KL()),
                        config_->Material());
}
template <>
inline StateIP Fiber<StateIP>::InnerModeT(const CS* objCS,
                                          const size_t& n) const {
  return config_->ModeT(LocalCS(), objCS,
                        EigenFunctor(Jn, od(n), config_->KT()),
                        config_->Material());
}
template <>
inline StateAP Fiber<StateAP>::InnerModeT(const CS* objCS,
                                          const size_t& n) const {
  return config_->ModeT(LocalCS(), objCS,
                        EigenFunctor(Jn, od(n), config_->KT()),
                        config_->Material());
}
template <typename T>
inline void Fiber<T>::add_node() {
  node_.reserve(config_->NoN());
  for (auto& i : config_->Node()) node_.emplace_back(i, this->LocalCS());
}
template <typename T>
inline Eigen::MatrixXcd Fiber<T>::ModeMatrix(
    const Inhomogeneity<T>* source) const {
  size_t P = this->NoE(), N = source->NoC();
  Eigen::MatrixXcd M(P, N);
  for (size_t n = 0; n < N; n++) {
  }
}

}  // namespace mss

#endif
