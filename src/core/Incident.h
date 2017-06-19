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

// Incident base class and its derivatives.

#ifndef MSS_INCIDENT_H
#define MSS_INCIDENT_H

#include "Matrix.h"
#include "State.h"

namespace mss {

template <typename StateType>
class Incident {
 public:
  Incident(const Matrix& matrix, const double& amplitude = 1,
           const double& phase = 0)
      : matrix_(matrix), amp_(amplitude), phase_(phase) {}

  virtual ~Incident() {}

  virtual StateType Effect(const CS* position) const = 0;

  double Amplitude() const { return amp_; }
  double Phase() const { return phase_; }

  // protected:
  const Matrix& matrix_;
  double amp_;
  double phase_;
};

template <typename StateType>
class IncidentPlane : public Incident<StateType> {
 public:
  IncidentPlane(const Matrix& matrix, const double& angle = 0,
                const double& amplitude = 1, const double& phase = 0)
      : Incident<StateType>(matrix, amplitude, phase), angle_(angle) {
    c_ = cos(angle_);
    s_ = sin(angle_);
  }

  virtual StateType Effect(const CS* localCS) const = 0;

  const double& Angle() const { return angle_; }

 protected:
  double angle_;
  double c_, s_;

  using Incident<StateType>::matrix_;
  using Incident<StateType>::phase_;

  dcomp _PhaseGLB(const CS* localCS, const double& k) const {
    PosiVect positionGLB = localCS->inGLB().Position();
    double x = positionGLB.x, y = positionGLB.y;
    return exp(ii * (k * (c_ * x + s_ * y) + phase_));
  }
  StateIP _StateGLB(const dcomp& u, const dcomp& v, const double& k) const {
    // Return the State from known displacement u and v. (The wave number is
    // needed since P-wave and SV-wave have different wave numbers.)

    const double &l = matrix_.Lambda(), &m = matrix_.Mu();
    dcomp gxx = ii * k * c_ * u;
    dcomp gyy = ii * k * s_ * v;
    dcomp gxy = ii * k * (s_ * u + c_ * v);
    return StateIP(u, v, (l + 2 * m) * gxx + l * gyy,
                   (l + 2 * m) * gyy + l * gxx, m * gxy);
  }
  StateAP _StateGLB(const dcomp& w) const {
    // Return the State from known displacement w.

    double km = matrix_.KT() * matrix_.Mu();
    return StateAP(w, ii * km * c_ * w, ii * km * s_ * w);
  }
};

class IncidentPlaneP : public IncidentPlane<StateIP> {
 public:
  IncidentPlaneP(const Matrix& matrix, const double& angle = 0,
                 const double& amplitude = 1, const double& phase = 0)
      : IncidentPlane<StateIP>(matrix, angle, amplitude, phase) {}

  StateIP Effect(const CS* localCS) const override {
    // The effect of the incident plane P-wave at the point with localCS.

    dcomp e = _PhaseGLB(localCS, matrix_.KL());
    return _StateGLB(amp_ * c_ * e, amp_ * s_ * e, matrix_.KL()).in(localCS);
  }
};

class IncidentPlaneSV : public IncidentPlane<StateIP> {
 public:
  IncidentPlaneSV(const Matrix& matrix, const double& angle = 0,
                  const double& amplitude = 1, const double& phase = 0)
      : IncidentPlane<StateIP>(matrix, angle, amplitude, phase) {}

  StateIP Effect(const CS* localCS) const override {
    // The effect of the incident plane SV-wave at the point with localCS.

    dcomp e = _PhaseGLB(localCS, matrix_.KT());
    return _StateGLB(-amp_ * s_ * e, amp_ * c_ * e, matrix_.KT()).in(localCS);
  }
};

class IncidentPlaneSH : public IncidentPlane<StateAP> {
 public:
  IncidentPlaneSH(const Matrix& matrix, const double& angle = 0,
                  const double& amplitude = 1, const double& phase = 0)
      : IncidentPlane<StateAP>(matrix, angle, amplitude, phase) {}

  StateAP Effect(const CS* localCS) const override {
    // The effect of the incident plane SH-wave at the point with localCS.

    dcomp e = _PhaseGLB(localCS, matrix_.KT());
    return _StateGLB(amp_ * e).in(localCS);
  }
};

}  // namespace mss

#endif
