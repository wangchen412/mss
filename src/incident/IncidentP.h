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

#ifndef MSS_INCIDENTP_H
#define MSS_INCIDENTP_H

#include "Incident.h"

namespace mss {

class IncidentP : virtual public Incident<IP> {
 public:
  IncidentP(const Matrix& matrix, double amplitude = 1, double phase = 0)
      : Incident<IP>(matrix, amplitude, phase), k_(matrix.KL()) {}

  virtual StateIP Effect(const PosiVect& position) const = 0;
  virtual StateIP Effect(const CS* localCS) const = 0;

 protected:
  double k_;
};

}  // namespace mss

#endif
