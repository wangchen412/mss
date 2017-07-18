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

#include "../core/Matrix.h"
#include "../core/State.h"
#include "../pre/Input.h"

namespace mss {

template <typename T>
class ConfigFiber {
 public:
  ConfigFiber(const input::ConfigFiber& input,
              const Matrix* matrix)
      : ID_(ID),
        N_(input.N_max),
        NoC_(2 * N_ + 1 * T::NoBV / 2),
        P_(input.P),
        R_(input.radius),
        material_(*input.material),
        kl_(matrix->Frequency() / material_.CL()),
        kt_(matrix->Frequency() / material_.CT()),
        matrix_(matrix),
        Q_(NoE(), NoC()) {
    add_node();
    compute_MatrixQ();
  }

  const Eigen::MatrixXcd& TransMatrix() const { return Q_; }
  const double& CharLength() const { return R_; }
  size_t NoN() const { return P_; }
  size_t NoE() const { return P_ * T::NoBV; }
  size_t NoC() const { return NoC_; }
  int TopOrder() const { return N_; }

  const std::string& ID() const { return ID_; }
  const double& Radius() const { return R_; }
  const class Material& Material() const { return material_; }
  const double& KL() const { return kl_; }
  const double& KT() const { return kt_; }
  const class Matrix* Matrix() const { return matrix_; }
  const class Material& Material_m() const { return matrix_->Material(); }
  const double& KL_m() const { return matrix_->KL(); }
  const double& KT_m() const { return matrix_->KT(); }
  const std::vector<CS>& Node() const { return node_; }

  T ModeT(const CS* localCS, const CS* objCS, const EigenFunctor& f,
          const class Material& mat) const;
  T ModeL(const CS* localCS, const CS* objCS, const EigenFunctor& f,
          const class Material& mat) const;

  // Functions for the factor T_n: B_n = T_n A_n.
  // tL is for longitude modes and tT is for transverse modes.
  dcomp TL(int n) const;  // TODO: T-matrix for in-plane problem.
  dcomp TT(int n) const;

 protected:
  const std::string ID_;           // The ID.
  const int N_;                    // The top order of the series.
  const int NoC_;                  // Number of the scattering coefficients.
  const size_t P_;                 // Number of the collocation points.
  const double R_;                 // Radius of the fiber.
  const class Material material_;  // Material of the fiber.
  const double kl_, kt_;           // Wave numbers of the fiber.
  const class Matrix* matrix_;     // The matrix.
  std::vector<CS> node_;           // Collocation points.
  Eigen::MatrixXcd Q_;             // Transform matrix.

  void add_node();
  void compute_MatrixQ();
};

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T>
inline void ConfigFiber<T>::add_node() {
  node_.reserve(P_);
  for (size_t i = 0; i < P_; i++)
    node_.emplace_back(PosiVect(R_, i * pi2 / P_).Cartesian(), i * pi2 / P_);
}

}  // namespace mss

#endif
