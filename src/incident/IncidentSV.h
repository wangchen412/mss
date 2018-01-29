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

#ifndef MSS_INCIDENTSV_H
#define MSS_INCIDENTSV_H

#include "IncidentPlane.h"
#include "IncidentS.h"

namespace mss {

class IncidentPlaneSV : public IncidentPlane<IP>, public IncidentS<IP> {
 public:
  IncidentPlaneSV(const Matrix& matrix, double angle = 0,
                  double amplitude = 1, double phase = 0)
      : Incident<IP>(matrix, amplitude, phase),
        IncidentPlane<IP>(matrix, angle, amplitude, phase),
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

}  // namespace mss

#endif
