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

#include "Input.h"

namespace mss {

namespace input {

Solution::Solution(const std::string& file) : fn_(file) {
  add_keyword();
  add(material_, matrix_, incident_, configFiber_, configAssembly_, solve_);
  link();
}

void Solution::link() {
  for (auto& i : material_) {
    i.cl = std::sqrt((i.lambda + 2 * i.mu) / i.rho);
    i.ct = std::sqrt(i.mu / i.rho);
  }
  for (auto& i : matrix_) {
    i.material = FindID(material_, i.materialID);
    i.kl       = i.frequency / i.material->cl;
    i.kt       = i.frequency / i.material->ct;
  }
  for (auto& i : configFiber_) {
    i.material = FindID(material_, i.materialID);
    i.P = std::max(size_t(matrix().kt * i.radius * matrix().delta), P_MIN);
  }
  for (auto& i : configAssembly_) {
    i.configFiber = &configFiber_;
    for (auto& j : i.fiber) j.config = FindID(configFiber_, j.configID);
  }

}

std::ostream& Solution::Print(std::ostream& os) const {
  return print(os, material_, matrix_, incident_, configFiber_,
               configAssembly_, solve_);
}

void Solution::add_keyword() {
  keyword_[typeid(Material)]       = "[Materials]";
  keyword_[typeid(Matrix)]         = "[Matrix]";
  keyword_[typeid(IncidentPlane)]  = "[Incident Waves]";
  keyword_[typeid(ConfigFiber)]    = "[Fiber Configurations]";
  keyword_[typeid(Fiber)]          = "[Fibers]";
  keyword_[typeid(ConfigAssembly)] = "[Assembly Configurations]";
  keyword_[typeid(std::string)]    = "[Solve]";
}

}  // namespace input

}  // namespace mss
