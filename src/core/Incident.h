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

#include "Input.h"
#include "Matrix.h"
#include "State.h"

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
  Eigen::VectorXcd Effect(const std::vector<CS*>& localCS) {
    Eigen::VectorXcd rst(localCS.size() * T::NoBV);
    for (size_t i = 0; i < localCS.size(); i++)
      rst.segment(i * T::NoBV, T::NoBV) = Effect(localCS[i]).DispTracVect();
    return rst;
  }

  const double& Amplitude() const { return amp_; }
  const double& Phase() const { return phase_; }

 protected:
  double amp_, phase_;
  const Material& m;
};

class IncidentP : virtual public Incident<StateIP> {
 public:
  IncidentP(const Matrix& matrix, const double& amplitude = 1,
            const double& phase = 0)
      : Incident<StateIP>(matrix, amplitude, phase), k_(matrix.KL()) {}

  virtual StateIP Effect(const PosiVect& position) const = 0;
  virtual StateIP Effect(const CS* localCS) const = 0;

 protected:
  const double& k_;
};

template <typename T>
class IncidentS : virtual public Incident<T> {
 public:
  IncidentS(const Matrix& matrix, const double& amplitude = 1,
            const double& phase = 0)
      : Incident<T>(matrix, amplitude, phase), k_(matrix.KT()) {}

  virtual T Effect(const PosiVect& position) const = 0;
  virtual T Effect(const CS* localCS) const = 0;

 protected:
  const double& k_;
};

template <typename T>
class IncidentPlane : virtual public Incident<T> {
 public:
  IncidentPlane(const Matrix& matrix, const double& angle = 0,
                const double& amplitude = 1, const double& phase = 0)
      : Incident<T>(matrix, amplitude, phase), angle_(angle) {
    c_ = cos(angle_);
    s_ = sin(angle_);
  }

  virtual T Effect(const PosiVect& position) const = 0;
  virtual T Effect(const CS* localCS) const = 0;

  const double& Angle() const { return angle_; }

 protected:
  double angle_;
  double c_, s_;

  using Incident<T>::phase_;

  dcomp _phaseGLB(const PosiVect& position, const double& k) const {
    // Return the phase of the incident wave at the position in global CS.

    return exp(ii * (k * (c_ * position.x + s_ * position.y) + phase_));
  }
  StateIP _stateGLB(const dcomp& u, const dcomp& v, const double& k) const {
    // Return the State from known displacement u and v. (The wave number is
    // needed since P-wave and SV-wave have different wave numbers.)

    dcomp gxx = ii * k * c_ * u;
    dcomp gyy = ii * k * s_ * v;
    dcomp gxy = ii * k * (s_ * u + c_ * v);
    StressIP t = this->m.C(gxx, gyy, gxy);
    return StateIP(u, v, t);
  }
  StateAP _stateGLB(const dcomp& w, const double& k) const {
    // Return the State from known displacement w.

    StressAP t = this->m.C(c_, s_) * ii * k * w;
    return StateAP(w, t);
  }
};

class IncidentPlaneP : public IncidentPlane<StateIP>, public IncidentP {
 public:
  IncidentPlaneP(const Matrix& matrix, const double& angle = 0,
                 const double& amplitude = 1, const double& phase = 0)
      : Incident<StateIP>(matrix, amplitude, phase),
        IncidentPlane<StateIP>(matrix, angle, amplitude, phase),
        IncidentP(matrix, amplitude, phase) {}
  IncidentPlaneP(const Matrix& matrix, const input::IncidentPlane& input)
      : IncidentPlaneP(matrix, input.angle, input.amplitude, input.phase) {}

  StateIP Effect(const PosiVect& position) const override {
    // The effect of the incident plane P-wave at the point with the poisition
    // in the global CS.

    dcomp e = _phaseGLB(position, k_);
    return _stateGLB(amp_ * c_ * e, amp_ * s_ * e, k_);
  }
  StateIP Effect(const CS* localCS) const override {
    // The effect of the incident plane P-wave at the point with the localCS.

    return Effect(localCS->PositionGLB()).in(localCS);
  }
};

class IncidentPlaneSV : public IncidentPlane<StateIP>,
                        public IncidentS<StateIP> {
 public:
  IncidentPlaneSV(const Matrix& matrix, const double& angle = 0,
                  const double& amplitude = 1, const double& phase = 0)
      : Incident<StateIP>(matrix, amplitude, phase),
        IncidentPlane<StateIP>(matrix, angle, amplitude, phase),
        IncidentS(matrix, amplitude, phase) {}
  IncidentPlaneSV(const Matrix& matrix, const input::IncidentPlane& input)
      : IncidentPlaneSV(matrix, input.angle, input.amplitude, input.phase) {}

  StateIP Effect(const PosiVect& position) const override {
    // The effect of the incident plane S-wave at the point with the position
    // in the global CS.

    dcomp e = _phaseGLB(position, k_);
    return _stateGLB(-amp_ * s_ * e, amp_ * c_ * e, k_);
  }
  StateIP Effect(const CS* localCS) const override {
    // The effect of the incident plane SV-wave at the point with the localCS.

    return Effect(localCS->PositionGLB()).in(localCS);
  }
};

class IncidentPlaneSH : public IncidentPlane<StateAP>,
                        public IncidentS<StateAP> {
 public:
  IncidentPlaneSH(const Matrix& matrix, const double& angle = 0,
                  const double& amplitude = 1, const double& phase = 0)
      : Incident<StateAP>(matrix, amplitude, phase),
        IncidentPlane<StateAP>(matrix, angle, amplitude, phase),
        IncidentS(matrix, amplitude, phase) {}
  IncidentPlaneSH(const Matrix& matrix, const input::IncidentPlane& input)
      : IncidentPlaneSH(matrix, input.angle, input.amplitude, input.phase) {}

  StateAP Effect(const PosiVect& position) const override {
    // The effect of the incident plane S-wave at the point with the position
    // in the global CS.

    dcomp e = _phaseGLB(position, k_);
    return _stateGLB(amp_ * e, k_);
  }
  StateAP Effect(const CS* localCS) const override {
    // The effect of the incident plane SH-wave at the point with localCS.

    return Effect(localCS->PositionGLB()).in(localCS);
  }
};

// ---------------------------------------------------------------------------
// Input functor:

template <typename T>
class InputIncident {
 public:
  InputIncident(const Matrix& matrix) : matrix_(matrix) {
    funcMap["PlaneP"] = [this](const input::IncidentPlane& input) {
      return new IncidentPlaneP(matrix_, input);
    };
    funcMap["PlaneSV"] = [this](const input::IncidentPlane& input) {
      return new IncidentPlaneSV(matrix_, input);
    };
    funcMap["PlaneSH"] = [this](const input::IncidentPlane& input) {
      return new IncidentPlaneSH(matrix_, input);
    };
  }

  Incident<T>* operator()(const input::IncidentPlane& input) {
    return funcMap[input.type](input);
  }

 private:
  typedef std::function<Incident<T>*(const input::IncidentPlane&)> funcType;
  std::map<std::string, funcType, ci_comp> funcMap;
  const Matrix& matrix_;
};

}  // namespace mss

#endif
