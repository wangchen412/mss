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

// Matrix class. Matrix here is the binder component in heterogeneous
// materials. The matrix object records not only the matrix material
// constants, but also the incident frequency, which is the only
// frequency in the whole case. So that the wave number of the matrix
// can be provided.

#ifndef MSS_MATRIX_H
#define MSS_MATRIX_H

#include "../pre/Input.h"
#include "Material.h"

namespace mss {

class Matrix {
 public:
  Matrix(const Material& material, double frequency)
      : material_(material),
        omega_(frequency),
        kl_(omega_ / material_.CL()),
        kt_(omega_ / material_.CT()) {
    assert(frequency > 0);
  }
  explicit Matrix(const input::Matrix& input)
      : Matrix(*input.material, input.frequency) {
    assert(input.kl == kl_ && input.kt == kt_);
  }
  explicit Matrix(const input::Solution& input) : Matrix(input.matrix()) {}

  virtual ~Matrix() {}

  const class Material& Material() const { return material_; }
  double Frequency() const { return omega_; }
  double KL() const { return kl_.real(); }
  double KT() const { return kt_.real(); }

  dcomp KL_comp() const { return kl_; }
  dcomp KT_comp() const { return kt_; }

 private:
  const class Material material_;
  const double omega_;
  const dcomp kl_, kt_;
};

}  // namespace mss

#endif
