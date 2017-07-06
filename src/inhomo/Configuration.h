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
  explicit Configuration(const Matrix* matrix) : matrix_(matrix) {}
  virtual ~Configuration() {}

  virtual const std::string& ID() const = 0;
  virtual const Eigen::MatrixXcd& TransMatrix() const = 0;

  virtual const double& CharLength() const = 0;
  virtual int NoN() const = 0;  // Number of collocation points (node).
  virtual int NoC() const = 0;  // Number of unknown coefficients.
  virtual int NoE() const = 0;  // Number of equations.

  const Matrix* Matrix() const { return matrix_; }
  const Material& Material_m() const { return matrix_->Material(); }
  const double& KL_m() const { return matrix_->KL(); }
  const double& KT_m() const { return matrix_->KT(); }

 protected:
  const class Matrix* matrix_;  // The matrix.
};

}  // namespace mss

#endif
