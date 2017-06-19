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

// State class.
// Combining displacement, stress and the coordinate system information.

#ifndef MSS_STATE_H
#define MSS_STATE_H

#include "CS.h"

namespace mss {

template <typename T1, typename T2>
class State {
 public:
  State() : displacement_(), stress_(), basis_(nullptr) {}
  explicit State(const T1& disp, const T2& stress, const CS* basis = nullptr)
      : displacement_(disp), stress_(stress), basis_(basis) {}
  explicit State(const dcomp& w, const dcomp& szx, const dcomp& szy,
                 const CS* basis = nullptr)
      : displacement_(w), stress_(szx, szy), basis_(basis) {}
  explicit State(const dcomp& u, const dcomp& v, const dcomp& sxx,
                 const dcomp& syy, const dcomp& sxy,
                 const CS* basis = nullptr)
      : displacement_(u, v), stress_(sxx, syy, sxy), basis_(basis) {}
  State(const State& other)
      : displacement_(other.displacement_),
        stress_(other.stress_),
        basis_(other.basis_) {}
  virtual ~State() {}

  bool operator==(const State& other) const;

  const T1& Displacement() const { return displacement_; }
  const T2& Stress() const { return stress_; }
  const CS* Basis() const { return basis_; }

  State in(const CS* otherBasis) const;

 private:
  T1 displacement_;
  T2 stress_;
  const CS* basis_;
};

typedef State<DispAP, StressAP> StateAP;
typedef State<DispIP, StressIP> StateIP;

// ---------------------------------------------------------------------------
// Inline functions:

template <>
inline StateAP StateAP::in(const CS* otherBasis) const {
  double d = otherBasis->in(basis_).Angle();
  return StateAP(displacement_, stress_.Rotate(d), otherBasis);
}

template <>
inline StateIP StateIP::in(const CS* otherBasis) const {
  double d = otherBasis->in(basis_).Angle();
  return StateIP(displacement_.Rotate(d), stress_.Rotate(d), otherBasis);
}

template <typename T1, typename T2>
inline bool State<T1, T2>::operator==(const State<T1, T2>& other) const {
  return (displacement_ == other.displacement_) &&
         (stress_ == other.stress_) && (basis_ == other.basis_);
}

}  // namespace mss

#endif
