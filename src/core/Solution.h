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

#ifndef MSS_SOLUTION_H
#define MSS_SOLUTION_H

#include "../incident/IncidentInput.h"
#include "../inhomo/Assembly.h"
#include "../pre/Input.h"

namespace mss {

template <typename T>
class Solution {
 public:
  Solution(const input::Solution& input)
      : matrix_(input.matrix()), config_("root", input.config(), &matrix_) {
    add_incident(input);
  }

  virtual ~Solution() { delete_incident(); }

  void Solve() { config_.Solve(incident_); }
  Inhomogeneity<T>* InWhich(const CS* objCS) const;
  T Resultant(const CS* objCS, const Inhomogeneity<T>* inhomo) const;
  T Resultant(const CS* objCS) const;

  const ConfigAssembly<T>& Config() const { return config_; }
  const InciCPtrs<T>& Incident() const { return incident_; }
  const InhomoCPtrs<T>& Inhomo() const { return Config().Inhomo(); }

 protected:
  // bool solved_;
  const Matrix matrix_;
  InciCPtrs<T> incident_;
  T inciNorm_;  // The normalization factor of the incidents.

  // Only the configuration of the "root" assembly is needed.
  // The instantiation of the "root" assembly is not necessary.
  ConfigAssembly<T> config_;

  void add_incident(const input::Solution& input);
  void delete_incident();
};

typedef Solution<StateAP> SolutionAP;
typedef Solution<StateIP> SolutionIP;

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T>
void Solution<T>::add_incident(const input::Solution& input) {
  IncidentInput<T> f(matrix_);
  for (auto& i : input.incident()) incident_.push_back(f(i));
}

template <typename T>
void Solution<T>::delete_incident() {
  for (auto& i : incident_) delete i;
}

template <typename T>
Inhomogeneity<T>* Solution<T>::InWhich(const CS* objCS) const {
  return config_.InWhich(objCS);
}

template <typename T>
T Solution<T>::Resultant(const CS* objCS, const Inhomogeneity<T>* in) const {
  return config_.Resultant(objCS, in, incident_);
}

template <typename T>
T Solution<T>::Resultant(const CS* objCS) const {
  return config_.Resultant(objCS, InWhich(objCS), incident_);
}

}  // namespace mss

#endif
