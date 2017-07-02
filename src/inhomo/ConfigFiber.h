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
      : Configuration<T>(2 * input.N_max + 1, matrix),
        R_(input.radius),
        material_(input.material),
        kl_f(matrix->Frequency() / material_.CL()),
        kt_f(matrix->Frequency() / material_.CT()) {}

  const Eigen::MatrixXcd& TransMatrix() const override;

  virtual const double& CharLength() const override { return R_; }
  const double& Radius() const { return R_; }
  const Material& Material() const { return material_; }
  const double& KL() const { return kl_f; }
  const double& KT() const { return kt_f; }

 protected:
  const double R_;
  const class Material material_;
  const double kl_f, kt_f;
};

}  // namespace mss

#endif
