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

#ifndef MSS_INCIDENTPLANESH_H
#define MSS_INCIDENTPLANESH_H

#include "IncidentPlane.h"
#include "IncidentS.h"

namespace mss {

class IncidentPlaneSH : public IncidentPlane<AP>, public IncidentS<AP> {
 public:
  IncidentPlaneSH(const Matrix& matrix, double angle = 0,
                  double amplitude = 1, double phase = 0)
      : Incident<AP>(matrix, amplitude, phase),
        IncidentPlane<AP>(matrix, angle, amplitude, phase),
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

}  // namespace mss

#endif
