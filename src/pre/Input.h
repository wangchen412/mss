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

#include <algorithm>
#include "InputData.h"

namespace mss {

namespace input {

const size_t P_MIN = 300;

class Solution {
 public:
  Solution(const std::string& file) : fn_(file) {
    add_keyword();
    add(material_, matrix_, incident_, fiber_config_, assembly_config_,
        solve_);
    link();
  }

  std::ostream& Print(std::ostream& os) const;
  const std::string& FN() const { return fn_; }

  const std::vector<Material>& material() const { return material_; }
  const Matrix& matrix() const { return matrix_[0]; }
  double frequency() const { return matrix().frequency; }
  const std::vector<FiberConfig>& fiber_config() const {
    return fiber_config_;
  }
  const std::vector<IncidentPlane>& incident() const { return incident_; }
  const std::vector<AssemblyConfig>& assembly_config() const {
    return assembly_config_;
  }
  const AssemblyConfig& config() const {
    return *FindID(assembly_config_, solve_[0].configID);
  }
  SolveMethod method() const {
    if (iequals(solve_[0].method, "Collocation"))
      return COLLOCATION;
    else if (iequals(solve_[0].method, "DFT"))
      return DFT;
    else
      error_msg({"Unknown method: ", solve_[0].method, "."});
    exit(EXIT_FAILURE);
  }

 private:
  std::vector<Material> material_;
  std::vector<Matrix> matrix_;
  std::vector<IncidentPlane> incident_;
  std::vector<FiberConfig> fiber_config_;
  std::vector<AssemblyConfig> assembly_config_;
  std::vector<Solve> solve_;

  std::string fn_;
  std::map<std::type_index, std::string> keyword_;
  std::map<std::type_index, std::string> header_;

  // Add parsing keywords to the dictionary.
  void add_keyword();
  bool is_keyword(const std::string& val);

  template <typename T>
  void add_header(std::ifstream& file);

  // Add entries of vec's element type to the back of vec.
  // The position of the input file stream is after the header.
  // For the AssemblyConfig class, which is partially specialized, each entry
  // includes multiple Fiber entries.
  template <typename T>
  void add_entry(std::ifstream& file, std::vector<T>& vec);

  template <typename T>
  void add(std::vector<T>& vec);
  template <typename T, typename... Ts>
  void add(std::vector<T>& vec, std::vector<Ts>&... vecs);

  template <typename T>
  std::ostream& print(std::ostream& os, const std::vector<T>& vec) const;
  template <typename T, typename... Ts>
  std::ostream& print(std::ostream& os, const std::vector<T>& vec,
                      const std::vector<Ts>&... vecs) const;
  void link();
};

// ---------------------------------------------------------------------------
// Inline functions:

bool Solution::is_keyword(const std::string& val) {
  for (auto& i : keyword_) {
    if (i.first == typeid(Fiber) || i.first == typeid(Assembly)) continue;
    if (iequals(i.second, val)) return true;
  }
  return false;
}

template <typename T>
void Solution::add_header(std::ifstream& file) {
  std::string tmp;
  skipUntil(file, keyword_[typeid(T)], 2, &tmp);
  header_[typeid(T)] = tmp;
}

template <typename T>
void Solution::add_entry(std::ifstream& file, std::vector<T>& vec) {
  std::string tmp;
  while (getline(file, tmp)) {
    if (isWhiteSpace(tmp)) break;
    if (tmp[0] == '#') continue;
    vec.emplace_back();
    std::stringstream(tmp) >> vec.back();
  }
}

// Add AssemblyConfig:
template <>
inline void Solution::add_entry(std::ifstream& file,
                                std::vector<AssemblyConfig>& vec) {
  std::string tmp;
  bool new_config = false, read_extra = false;
  while (getline(file, tmp)) {
    if (is_keyword(tmp)) break;
    if (iequals(tmp.substr(0, 2), "ID")) {
      new_config = true;
      read_extra = true;
    }
    if (new_config) {
      new_config = false;
      AssemblyConfig rst;
      if (read_extra) getline(file, tmp);
      std::stringstream(tmp) >> rst.ID >> rst.width >> rst.height;
      while (getline(file, tmp)) {
        if (is_keyword(tmp)) break;
        if (iequals(tmp.substr(0, 2), "ID")) {
          new_config = true;
          read_extra = false;
          break;
        } else if (iequals(tmp, keyword_[typeid(Fiber)])) {
          skip(file, 2, &tmp);
          header_[typeid(Fiber)] = tmp;
          add_entry(file, rst.fiber);
        } else if (iequals(tmp, keyword_[typeid(Assembly)])) {
          skip(file, 2, &tmp);
          header_[typeid(Assembly)] = tmp;
          add_entry(file, rst.assembly);
        }
      }
      vec.push_back(rst);
    }
  }
}

template <typename T>
void Solution::add(std::vector<T>& vec) {
  std::ifstream file(fn_);
  add_header<T>(file);
  add_entry(file, vec);
  file.close();
}
template <typename T, typename... Ts>
void Solution::add(std::vector<T>& vec, std::vector<Ts>&... vecs) {
  add(vec);
  add(vecs...);
}

template <typename T>
std::ostream& Solution::print(std::ostream& os,
                              const std::vector<T>& vec) const {
  os << separator("=") << keyword_.at(typeid(T)) << std::endl
     << separator("-") << header_.at(typeid(T)) << std::endl;
  for (const auto& i : vec) os << i << std::endl;
  return os << std::endl << std::endl;
}
template <typename T, typename... Ts>
std::ostream& Solution::print(std::ostream& os, const std::vector<T>& vec,
                              const std::vector<Ts>&... vecs) const {
  print(os, vec);
  return print(os, vecs...);
}

inline void Solution::link() {
  for (auto& i : material_) {
    i.cl = std::sqrt((i.lambda + 2 * i.mu) / i.rho);
    i.ct = std::sqrt(i.mu / i.rho);
  }
  for (auto& i : matrix_) {
    i.material = FindID(material_, i.materialID);
    i.kl       = i.frequency / i.material->cl;
    i.kt       = i.frequency / i.material->ct;
  }
  for (auto& i : fiber_config_) {
    i.material = FindID(material_, i.materialID);
    i.P = std::max(size_t(matrix().kt * i.radius * matrix().delta), P_MIN);
  }
  for (auto& i : assembly_config_) {
    i.fiber_config = &fiber_config_;
    i.pointDensity = matrix().kt * matrix().delta;
    for (auto& j : i.fiber) j.config = FindID(fiber_config_, j.configID);
    for (auto& j : i.assembly)
      j.config = FindID(assembly_config_, j.configID);
  }
}

inline std::ostream& Solution::Print(std::ostream& os) const {
  return print(os, material_, matrix_, incident_, fiber_config_,
               assembly_config_, solve_);
}

inline void Solution::add_keyword() {
  keyword_[typeid(Material)]       = "[Materials]";
  keyword_[typeid(Matrix)]         = "[Matrix]";
  keyword_[typeid(IncidentPlane)]  = "[Incident Waves]";
  keyword_[typeid(FiberConfig)]    = "[Fiber Configurations]";
  keyword_[typeid(Fiber)]          = "[Fibers]";
  keyword_[typeid(AssemblyConfig)] = "[Assembly Configurations]";
  keyword_[typeid(Assembly)]       = "[Assemblies]";
  keyword_[typeid(Solve)]          = "[Solve]";
}

}  // namespace input

}  // namespace mss

#endif
