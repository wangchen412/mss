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

#ifndef MSS_INCIDENTPLANEP_H
#define MSS_INCIDENTPLANEP_H

#include "IncidentP.h"
#include "IncidentPlane.h"

namespace mss {

class IncidentPlaneP : public IncidentPlane<IP>, public IncidentP {
 public:
  IncidentPlaneP(const Matrix& matrix, double angle = 0, double amplitude = 1,
                 double phase = 0)
      : Incident<IP>(matrix, amplitude, phase),
        IncidentPlane<IP>(matrix, angle, amplitude, phase),
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

}  // namespace mss

#endif
