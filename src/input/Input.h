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
#include "../core/Tensor.h"
#include "../tools/FileIO.h"

namespace mss {

namespace input {

class Material {
 public:
  Material() {}
  Material(std::stringstream data);

  std::string ID;
  double rho, lambda, mu;
};
struct Matrix {
  Matrix() {}
  Matrix(const std::string& file) { read(file); }
  void read(const std::string& file);
  Material material;
  double frequency;
};
struct IncidentPlane {
  std::string type;
  double amplitude, phase, angle;
};
struct ConfigFiber {
  std::string ID;
  Material material;
  double radius;
  int N_max;
  size_t P;
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
struct Solution {
  Solution(const std::string& file) : fn_(file) {}

  Matrix matrix;
  std::vector<IncidentPlane> incident;
  ConfigAssembly config;

 private:
  std::string fn_;
  std::vector<Material*> material_;
  std::vector<ConfigFiber*> configFiber_;
  std::vector<Fiber*> fiber_;
  std::vector<ConfigAssembly*> configAssembly_;
  std::vector<Assembly*> assembly_;

  void read(const std::string& file);
  void add_material(const std::string& file);

  template <typename T>
  void add(std::vector<T*>& vec, const std::string& key);
};

}  // namespace input

}  // namespace mss

#endif
