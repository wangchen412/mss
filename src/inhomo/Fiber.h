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
      : Inhomogeneity<T>(config, position),
        config_(config),
        kl_f(config->KL()),
        kt_f(config->KT()) {
    add_node();
  }

  T ScatterModeL(const CS* objCS, int n) const override;
  T ScatterModeT(const CS* objCS, int n) const override;
  T InnerModeL(const CS* objCS, int n) const override;
  T InnerModeT(const CS* objCS, int n) const override;

  Eigen::VectorXcd InVector(
      const std::vector<Incident<T>*>& incident) const override;

  // Return the effect of modes of other source inhomogeneity at the
  // collocation points.
  Eigen::MatrixXcd ModeMatrix(const Inhomogeneity<T>* source) const override;

  const ConfigFiber<T>* Config() const override { return config_; }

 protected:
  const ConfigFiber<T>* config_;
  const double &kl_f, &kt_f;
  std::vector<CS> node_;

  void add_node();
};  // namespace mss

// ---------------------------------------------------------------------------
// Inline functions:

template <>
inline StateIP Fiber<StateIP>::ScatterModeL(const CS* objCS, int n) const {
  return Config()->ModeL(LocalCS(), objCS, EigenFunctor(Hn, n, kl_m),
                         Config()->Material_m());
}
template <>
inline StateIP Fiber<StateIP>::ScatterModeT(const CS* objCS, int n) const {
  return Config()->ModeT(LocalCS(), objCS, EigenFunctor(Hn, n, kt_m),
                         Config()->Material_m());
}
template <>
inline StateIP Fiber<StateIP>::InnerModeL(const CS* objCS, int n) const {
  return Config()->ModeL(LocalCS(), objCS, EigenFunctor(Jn, n, kl_f),
                         Config()->Material());
}
template <>
inline StateIP Fiber<StateIP>::InnerModeT(const CS* objCS, int n) const {
  return Config()->ModeT(LocalCS(), objCS, EigenFunctor(Jn, n, kt_f),
                         Config()->Material());
}
template <>
inline StateAP Fiber<StateAP>::InnerModeT(const CS* objCS, int n) const {
  return Config()->ModeT(LocalCS(), objCS, EigenFunctor(Jn, n, kt_f),
                         Config()->Material());
}
template <typename T>
inline void Fiber<T>::add_node() {
  node_.reserve(Config()->NoN());
  for (auto& i : Config()->Node()) node_.emplace_back(i, this->LocalCS());
}

}  // namespace mss

#endif
