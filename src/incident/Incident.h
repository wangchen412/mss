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

// Incident base class.

#ifndef MSS_INCIDENT_H
#define MSS_INCIDENT_H

#include "../core/Matrix.h"
#include "../core/State.h"
#include "../pre/Input.h"

namespace mss {

template <typename T>
class Incident {
 public:
  Incident(const Matrix& matrix, const double& amplitude = 1,
           const double& phase = 0)
      : amp_(amplitude), phase_(phase), m(matrix.Material()) {}

  virtual ~Incident() {}

  // The effect on the position in global CS.
  virtual T Effect(const PosiVect& position) const = 0;

  // The effect on the origin of the local CS.
  virtual T Effect(const CS* localCS) const = 0;

  // The effect vector along the boundary.
  Eigen::VectorXcd EffectBV(const std::vector<CS*>& localCS) const {
    Eigen::VectorXcd rst(localCS.size() * T::NoBV);
    for (size_t i = 0; i < localCS.size(); i++)
      rst.segment(i * T::NoBV, T::NoBV) = Effect(localCS[i]).BV();
    return rst;
  }

  const double& Amplitude() const { return amp_; }
  const double& Phase() const { return phase_; }

 protected:
  double amp_, phase_;
  const Material& m;
};

template <typename T>
using InciPtrs = std::vector<Incident<T>*>;

}  // namespace mss

#endif
