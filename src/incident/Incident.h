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
  Incident(const Matrix& matrix, double amplitude = 1, double phase = 0)
      : amp_(amplitude),
        phase_(phase),
        m(matrix.Material()),
        kl_(matrix.KL()),
        kt_(matrix.KT()) {}

  virtual ~Incident() {}

  // The normalization factor.
  T Norm() const;

  // The effect on the position in global CS.
  virtual T Effect(const PosiVect& position) const = 0;

  // The effect on the origin of the local CS.
  virtual T Effect(const CS* localCS) const = 0;

  // The effect vector along the boundary.
  VectorXcd EffectBv(const CSCPtrs& localCS) const {
    VectorXcd rst(localCS.size() * T::NumBv);
    for (size_t i = 0; i < localCS.size(); i++)
      rst.segment(i * T::NumBv, T::NumBv) = Effect(localCS[i]).Bv();
    return rst;
  }

  double Amplitude() const { return amp_; }
  double Phase() const { return phase_; }

 protected:
  double amp_, phase_;
  const Material& m;

 private:
  const double kl_, kt_;
};

template <typename T>
using InciPtrs = std::vector<Incident<T>*>;

template <typename T>
using InciCPtrs = std::vector<const Incident<T>*>;

// ---------------------------------------------------------------------------
// Inline functions:

template <>
StateIP Incident<IP>::Norm() const {
  return StateIP(1, 1, m.C(kl_, kl_, kl_)) * amp_;
}
template <>
StateAP Incident<AP>::Norm() const {
  return StateAP(1, m.C(kt_, kt_)) * amp_;
}

}  // namespace mss

#endif
