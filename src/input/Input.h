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
struct Assembly {
  std::string configID;
  PosiVect position;
  double angle;
};
struct ConfigAssembly {
  std::string ID;
  std::vector<ConfigFiber> configFiber;
  std::vector<Fiber> fiber;
  std::vector<ConfigAssembly> a;
  std::vector<Assembly> assembly;
};

class Solution {
 public:
  Solution(const std::string& file) : fn_(file) {
    add_keyword();
    add(material_);
    add(matrix_);
    add(incident_);
    add(configFiber_);
  }

  friend std::ostream& operator<<(std::ostream& os, const Solution& s);

 private:
  std::string fn_;
  std::vector<Material> material_;
  std::vector<Matrix> matrix_;
  std::vector<IncidentPlane> incident_;
  std::vector<ConfigFiber> configFiber_;
  std::vector<ConfigAssembly> configAssembly_;
  std::vector<Assembly> assembly_;
  std::map<std::type_index, std::string> keyword_;
  void add_keyword();

  template <typename T>
  void add(std::vector<T>& vec);
};

// ---------------------------------------------------------------------------
// Inline functions:

inline void operator>>(std::istream& is, Material& m) {
  is >> m.ID >> m.rho >> m.mu >> m.lambda;
}
inline std::ostream& operator<<(std::ostream& os, const Material& m) {
  return os << m.ID << "\t" << m.rho << "\t" << m.mu << "\t" << m.lambda;
}

inline void operator>>(std::istream& is, Matrix& m) {
  is >> m.frequency >> m.delta >> m.materialID;
}
inline std::ostream& operator<<(std::ostream& os, const Matrix& m) {
  return os << m.frequency << "\t" << m.delta << "\t" << m.materialID;
}

inline void operator>>(std::istream& is, IncidentPlane& i) {
  is >> i.type >> i.angle >> i.amplitude >> i.phase;
}
inline std::ostream& operator<<(std::ostream& os, const IncidentPlane& i) {
  return os << i.type << "\t" << i.angle << "\t" << i.amplitude << "\t"
            << i.phase;
}

inline void operator>>(std::istream& is, ConfigFiber& c) {
  is >> c.ID >> c.radius >> c.N_max >> c.materialID;
}
inline std::ostream& operator<<(std::ostream& os, const ConfigFiber& c) {
  return os << c.ID << "\t" << c.radius << "\t" << c.N_max << "\t"
            << c.materialID;
}

inline void operator>>(std::istream& is, Fiber& f) {
  is >> f.position >> f.configID;
}
inline std::ostream& operator<<(std::ostream& os, const Fiber& f) {
  return os << f.position << "\t" << f.configID;
}

inline std::ostream& operator<<(std::ostream& os, const Solution& s) {
  for (const auto& i : s.material_) os << i << std::endl;
  for (const auto& i : s.matrix_) os << i << std::endl;
  for (const auto& i : s.incident_) os << i << std::endl;
  for (const auto& i : s.configFiber_) os << i << std::endl;
  return os;
}

void Solution::add_keyword() {
  keyword_[typeid(Material)]       = "[Materials]";
  keyword_[typeid(Matrix)]         = "[Matrix]";
  keyword_[typeid(IncidentPlane)]  = "[Incident Waves]";
  keyword_[typeid(ConfigFiber)]    = "[Fiber Configurations]";
  keyword_[typeid(ConfigAssembly)] = "[Assembly Configurations]";
}

template <typename T>
inline void Solution::add(std::vector<T>& vec) {
  std::ifstream file(fn_);
  skipUntil(file, keyword_[typeid(T)], 2);
  std::string tmp;
  while (getline(file, tmp)) {
    if (isWhiteSpace(tmp)) break;
    if (tmp[0] == '#') continue;
    vec.emplace_back();
    std::stringstream(tmp) >> vec.back();
  }
  file.close();
}

}  // namespace input

}  // namespace mss

#endif
