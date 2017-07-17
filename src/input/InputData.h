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

#ifndef MSS_INPUTDATA_H
#define MSS_INPUTDATA_H

#include <string>
#include <typeindex>
#include <vector>
#include "../core/Tensor.h"
#include "../tools/FileIO.h"

namespace mss {

namespace input {

struct Material {
  std::string ID;
  double rho, lambda, mu;
};
struct Matrix {
  double frequency, delta;
  std::string materialID;
};
struct IncidentPlane {
  std::string type;
  double amplitude, phase, angle;
};
struct ConfigFiber {
  std::string ID;
  std::string materialID;
  double radius;
  int N_max;
};
struct Fiber {
  std::string configID;
  PosiVect position;
};
// TODO
// struct Assembly {
//   std::string configID;
//   PosiVect position;
//   double angle;
// };
struct ConfigAssembly {
  std::string ID;
  std::vector<Fiber> fiber;
  // std::vector<ConfigAssembly> a;   // TODO
  // std::vector<Assembly> assembly;  // TODO
};

// ---------------------------------------------------------------------------
// Inline functions:

// I/O of the input data structs.

inline void operator>>(std::istream& is, Material& m) {
  is >> m.ID >> m.rho >> m.mu >> m.lambda;
}
inline std::ostream& operator<<(std::ostream& os, const Material& m) {
  return os << m.ID << "\t\t" << m.rho << "\t\t" << m.mu << "\t" << m.lambda;
}

inline void operator>>(std::istream& is, Matrix& m) {
  is >> m.frequency >> m.delta >> m.materialID;
}
inline std::ostream& operator<<(std::ostream& os, const Matrix& m) {
  return os << m.frequency << "\t\t" << m.delta << "\t\t" << m.materialID;
}

inline void operator>>(std::istream& is, IncidentPlane& i) {
  is >> i.type >> i.angle >> i.amplitude >> i.phase;
}
inline std::ostream& operator<<(std::ostream& os, const IncidentPlane& i) {
  return os << i.type << "\t\t" << i.angle << "\t\t" << i.amplitude << "\t\t"
            << i.phase;
}

inline void operator>>(std::istream& is, ConfigFiber& c) {
  is >> c.ID >> c.radius >> c.N_max >> c.materialID;
}
inline std::ostream& operator<<(std::ostream& os, const ConfigFiber& c) {
  return os << c.ID << "\t\t" << c.radius << "\t\t" << c.N_max << "\t\t"
            << c.materialID;
}

inline void operator>>(std::istream& is, Fiber& f) {
  is >> f.position >> f.configID;
}
inline std::ostream& operator<<(std::ostream& os, const Fiber& f) {
  return os << f.position << "\t\t" << f.configID;
}

inline std::ostream& operator<<(std::ostream& os, const ConfigAssembly& c) {
  os << separator("=", 45) << "ID" << std::endl << c.ID << std::endl;
  os << separator("-", 45) << "[Fibers]" << std::endl << separator("-", 45);
  os << "X               Y               Configuration" << std::endl;
  for (const auto& i : c.fiber) os << i << std::endl;
  return os;
}

}  // namespace input

}  // namespace mss

#endif
