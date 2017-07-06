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
#include "Functors.h"

namespace mss {

template <typename T1, typename T2>
class State {
 public:
  State() : displacement_(), stress_(), basis_(nullptr) {}
  explicit State(const T1& disp, const T2& stress, const CS* basis = nullptr)
      : displacement_(disp), stress_(stress), basis_(basis) {}
  explicit State(const dcomp& w, const StressAP& stress,
                 const CS* basis = nullptr)
      : displacement_(w), stress_(stress), basis_(basis) {}
  explicit State(const dcomp& u, const dcomp& v, const StressIP& stress,
                 const CS* basis = nullptr)
      : displacement_(u, v), stress_(stress), basis_(basis) {}
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

  State& operator+=(const State& other) {
    displacement_ += other.displacement_;
    stress_ += other.stress_;
    return *this;
  }
  State& operator-=(const State& other) {
    displacement_ -= other.displacement_;
    stress_ -= other.stress_;
    return *this;
  }
  State& operator*=(const dcomp& n) {
    displacement_ *= n;
    stress_ *= n;
    return *this;
  }
  State& operator/=(const dcomp& n) {
    displacement_ /= n;
    stress_ /= n;
    return *this;
  }
  State operator+(const State& other) const { return State(*this) += other; }
  State operator-(const State& other) const { return State(*this) -= other; }
  State operator*(const dcomp& n) const { return State(*this) *= n; }
  State operator/(const dcomp& n) const { return State(*this) /= n; }
  bool operator==(const State& other) const {
    State tmp(in(other.basis_));
    return (tmp.displacement_ == other.displacement_) &&
           (tmp.stress_ == other.stress_);
  }
  friend std::ostream& operator<<(std::ostream& os, const State& st) {
    return os << st.AngleGLB() << "\t" << st.displacement_ << "\t"
              << st.stress_;
  }
  friend std::istream& operator>>(std::istream& is, State& st) {
    // Since the object is read from file, the only CS it can be based on is
    // the global CS. So after the displacement and stress are read, the state
    // is rotated reversely to fit it into the global CS.

    assert(st == State());  // Write in empty state only.
    double angle;
    is >> angle >> st.displacement_ >> st.stress_;
    st._rotate(-angle);
    return is;
  }

  const T1& Displacement() const { return displacement_; }
  const T2& Stress() const { return stress_; }
  const CS* Basis() const { return basis_; }

  double AngleGLB() const;
  State in(const CS* otherBasis) const;

 private:
  T1 displacement_;
  T2 stress_;
  const CS* basis_;
  State& _rotate(const double& angle);
};

typedef State<DispAP, StressAP> StateAP;
typedef State<DispIP, StressIP> StateIP;

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T1, typename T2>
inline State<T1, T2>& State<T1, T2>::_rotate(const double& angle) {
  displacement_.RotateInPlace(angle);
  stress_.RotateInPlace(angle);
  return *this;
}
template <typename T1, typename T2>
inline double State<T1, T2>::AngleGLB() const {
  if (basis_)
    return basis_->AngleGLB();
  else
    return 0;
}
template <typename T1, typename T2>
inline State<T1, T2> State<T1, T2>::in(const mss::CS* otherBasis) const {
  if (otherBasis == basis_) return *this;
  double d = 0;
  if (otherBasis) d = otherBasis->AngleGLB();
  return State<T1, T2>(displacement_, stress_, otherBasis)
      ._rotate(d - AngleGLB());
}

}  // namespace mss

#endif
