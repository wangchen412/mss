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

#ifndef MSS_INPUT_H
#define MSS_INPUT_H

#include <string>
#include <vector>

namespace mss {

namespace input {

struct Material {
  std::string ID;
  double rho, lambda, mu;
};
struct Matrix {
  std::string materialID;
  double frequency;
};
struct IncidentPlane {
  std::string ID;
  std::string type;
  double amplitude, phase, angle;
};
struct ConfigFiber {
  std::string ID;
  Material material;
  double radius;
  int N_max;
};
struct Fiber {
  std::string configID;
  double x, y;
};
struct Assembly {
  std::string configID;
  double x, y, angle;
};
struct ConfigAssembly {
  std::string ID;
  std::vector<Fiber> fiber;
  std::vector<Assembly> assembly;
};
struct Solution {
  Matrix matrix;
  std::vector<IncidentPlane> incidentPlane;
  Assembly assembly;
};

} // namespace input

} // namespace mss

#endif
