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

#ifndef MSS_INCIDENTPLANE_H
#define MSS_INCIDENTPLANE_H

#include "Incident.h"

namespace mss {

template <typename T>
class IncidentPlane : virtual public Incident<T> {
 public:
  IncidentPlane(const Matrix& matrix, double angle = 0, double amplitude = 1,
                double phase = 0)
      : Incident<T>(matrix, amplitude, phase), angle_(angle) {
    c_ = cos(angle_);
    s_ = sin(angle_);
  }

  virtual T Effect(const PosiVect& position) const = 0;
  virtual T Effect(const CS* localCS) const = 0;

  double Angle() const { return angle_; }

 protected:
  double angle_;
  double c_, s_;

  using Incident<T>::phase_;

  dcomp _phaseGLB(const PosiVect& position, double k) const {
    // Return the phase of the incident wave at the position in global CS.

    return exp(ii * (k * (c_ * position.x + s_ * position.y) + phase_));
  }
  StateIP _stateGLB(const dcomp& u, const dcomp& v, double k) const {
    // Return the State from known displacement u and v. (The wave number is
    // needed since P-wave and SV-wave have different wave numbers.)

    dcomp gxx = ii * k * c_ * u;
    dcomp gyy = ii * k * s_ * v;
    dcomp gxy = ii * k * (s_ * u + c_ * v);
    StressIP t = this->m.C(gxx, gyy, gxy);
    return StateIP(u, v, t);
  }
  StateAP _stateGLB(const dcomp& w, double k) const {
    // Return the State from known displacement w.

    std::cout << "w: " << w << std::endl;
    std::cout << "k: " << k << std::endl;
    std::cout << "1: " << this->m.C(c_, s_) << std::endl;
    std::cout << "2: " << this->m.C(c_, s_) * ii << std::endl;
    // std::cout << "3: " << this->m.C(c_, s_) * ii * k << std::endl;
    // std::cout << "3: " << this->m.C(c_, s_) * ii * w << std::endl;

    StressAP t = this->m.C(c_, s_) * ii * k * w;
    return StateAP(w, t);
  }
};

}  // namespace mss

#endif
