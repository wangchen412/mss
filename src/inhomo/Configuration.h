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

#ifndef MSS_CONFIGURATION_H
#define MSS_CONFIGURATION_H

#include "../core/Input.h"
#include "../core/Matrix.h"
#include "../core/State.h"

namespace mss {

template <typename T>
class Configuration {
 public:
  explicit Configuration(size_t N, const Matrix* matrix)
      : N_(N), matrix_(matrix) {}

  virtual const Eigen::MatrixXcd& TransMatrix() const = 0;

  virtual const double& CharLength() const = 0;

  const Matrix* Matrix() const { return matrix_; }
  const Material& Material_m() const { return matrix_->Material(); }
  const double& KL_m() const { return matrix_->KL(); }
  const double& KT_m() const { return matrix_->KT(); }

 protected:
  const size_t N_;              // Number of the unknown coefficients.
  const class Matrix* matrix_;  // The matrix.
};

}  // namespace mss

#endif
