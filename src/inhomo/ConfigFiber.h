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

#ifndef MSS_CONFIGFIBER_H
#define MSS_CONFIGFIBER_H

#include "Configuration.h"

namespace mss {

template <typename T>
class ConfigFiber : public Configuration<T> {
 public:
  ConfigFiber(const input::ConfigFiber& input, const Matrix* matrix)
      : Configuration<T>(matrix),
        ID_(input.ID),
        N_(input.N_max),
        P_(input.P),
        R_(input.radius),
        material_(input.material),
        kl_(matrix->Frequency() / material_.CL()),
        kt_(matrix->Frequency() / material_.CT()),
        Q_(NoE(), NoC()) {
    add_node();
    computeQ();
  }

  const Eigen::MatrixXcd& TransMatrix() const override;

  const double& CharLength() const override { return R_; }
  int NoN() const override;
  int NoE() const override;
  int NoC() const override;

  const std::string& ID() const override { return ID_; }
  const double& Radius() const { return R_; }
  const Material& Material() const { return material_; }
  const double& KL() const { return kl_; }
  const double& KT() const { return kt_; }
  const std::vector<CS>& Node() const { return node_; }

  T ModeT(const CS* localCS, const CS* objCS, const EigenFunctor& f,
          const class Material& mat) const;
  T ModeL(const CS* localCS, const CS* objCS, const EigenFunctor& f,
          const class Material& mat) const;

 protected:
  const std::string ID_;           // The ID.
  const int N_;                    // The top order of the series.
  const int P_;                    // Number of the collocation points.
  const double R_;                 // Radius of the fiber.
  const class Material material_;  // Material of the fiber.
  const double kl_, kt_;           // Wave numbers of the fiber.
  std::vector<CS> node_;           // Collocation points.
  Eigen::MatrixXcd Q_;             // Transform matrix.

  // Funtions for the factor T_n: B_n = T_n A_n.
  // tL is for longitude modes and tT is for transverse modes.
  dcomp tL(int n) const;
  dcomp tT(int n) const;

  void add_node();
  void computeQ();
};

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T>
inline int ConfigFiber<T>::NoN() const {
  return P_;
}
template <>
inline int ConfigFiber<StateIP>::NoC() const {
  return 4 * N_ + 2;
}
template <>
inline int ConfigFiber<StateAP>::NoC() const {
  return 2 * N_ + 1;
}
template <>
inline int ConfigFiber<StateIP>::NoE() const {
  return P_ * 4;
}
template <>
inline int ConfigFiber<StateAP>::NoE() const {
  return P_ * 2;
}
template <typename T>
inline const Eigen::MatrixXcd& ConfigFiber<T>::TransMatrix() const {
  // The transform matrix for the least square collocation.
  return Q_;
}
template <typename T>
inline void ConfigFiber<T>::add_node() {
  node_.reserve(P_);
  for (int i = 0; i < P_; i++)
    node_.emplace_back(PosiVect(R_, i * pi2 / P_).Cartesian(), i * pi2 / P_);
}

}  // namespace mss

#endif
