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

#ifndef MSS_INCIDENTINPUT_H
#define MSS_INCIDENTINPUT_H

#include "IncidentPlaneP.h"
#include "IncidentPlaneSH.h"
#include "IncidentPlaneSV.h"

namespace mss {

template <typename T>
class IncidentInput {
 public:
  IncidentInput(const Matrix& matrix);
  Incident<T>* operator()(const input::IncidentPlane& input) {
    return funcMap[input.type](input);
  }

 private:
  typedef std::function<Incident<T>*(const input::IncidentPlane&)> funcType;
  std::map<std::string, funcType, ci_comp> funcMap;
  const Matrix& matrix_;
};

template <>
IncidentInput<IP>::IncidentInput(const Matrix& matrix) : matrix_(matrix) {
  funcMap["PlaneP"] = [this](const input::IncidentPlane& input) {
    return new IncidentPlaneP(matrix_, input);
  };
  funcMap["PlaneSV"] = [this](const input::IncidentPlane& input) {
    return new IncidentPlaneSV(matrix_, input);
  };
}

template <>
IncidentInput<AP>::IncidentInput(const Matrix& matrix) : matrix_(matrix) {
  funcMap["PlaneSH"] = [this](const input::IncidentPlane& input) {
    return new IncidentPlaneSH(matrix_, input);
  };
}

}  // namespace mss

#endif
