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

#include "../core/Tensor.h"
#include "../tools/FileIO.h"

namespace mss {

namespace input {

struct Material {
  std::string ID;
  double rho, lambda, mu;
  double cl, ct;
};
struct Matrix {
  double frequency, delta;
  std::string materialID;
  const Material* material;
  double kl, kt;
};
struct IncidentPlane {
  std::string type;
  double amplitude, phase, angle;
};
struct FiberConfig {
  std::string ID;
  std::string materialID;
  double radius;
  int N_max;
  size_t P;
  const Material* material;
};
struct Fiber {
  std::string configID;
  PosiVect position;
  const FiberConfig* config;
};
// TODO
// struct Assembly {
//   std::string configID;
//   PosiVect position;
//   double angle;
// };
struct AssemblyConfig {
  std::string ID;
  double pointDensity;  // Point density.
  double width, height;
  std::vector<Fiber> fiber;
  const std::vector<FiberConfig>* fiber_config;
  bool nsolve;
  // std::vector<AssemblyConfig> a;   // TODO
  // std::vector<Assembly> assembly;  // TODO
};
struct Solve {
  std::string configID;
  std::string method;
};

// ---------------------------------------------------------------------------
// Inline functions:

// I/O of the input data structs.

inline void operator>>(std::istream& is, Material& m) {
  is >> m.ID >> m.rho >> m.lambda >> m.mu;
}
inline std::ostream& operator<<(std::ostream& os, const Material& m) {
  return os << m.ID << "\t\t" << m.rho << "\t\t" << m.lambda << "\t" << m.mu;
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

inline void operator>>(std::istream& is, FiberConfig& c) {
  is >> c.ID >> c.radius >> c.N_max >> c.materialID;
}
inline std::ostream& operator<<(std::ostream& os, const FiberConfig& c) {
  return os << c.ID << "\t\t" << c.radius << "\t\t" << c.N_max << "\t\t"
            << c.materialID;
}

inline void operator>>(std::istream& is, Fiber& f) {
  is >> f.position >> f.configID;
}
inline std::ostream& operator<<(std::ostream& os, const Fiber& f) {
  return os << f.position << "\t\t" << f.configID;
}

inline std::ostream& operator<<(std::ostream& os, const AssemblyConfig& c) {
  os << separator("=", 45);
  os << "ID              Width           Height" << std::endl;
  os << c.ID << "\t\t" << c.width << "\t\t" << c.height << std::endl;
  os << separator("-", 45) << "[Fibers]" << std::endl << separator("-", 45);
  os << "X               Y               Configuration" << std::endl;
  for (const auto& i : c.fiber) os << i << std::endl;
  return os;
}

inline void operator>>(std::istream& is, Solve& s) {
  is >> s.configID >> s.method;
}
inline std::ostream& operator<<(std::ostream& os, const Solve& s) {
  return os << s.configID << "\t\t" << s.method;
}

}  // namespace input

}  // namespace mss

#endif
