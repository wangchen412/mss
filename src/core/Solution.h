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

#include "../inhomo/Assembly.h"
#include "../pre/Input.h"
#include "Incident.h"

namespace mss {

template <typename T>
class Solution {
 public:
  Solution(const input::Solution& input)
      : matrix_(input.matrix()),
        method_(input.method()),
        input_file_(input.FN()),
        original_(true) {
    config_ = new AssemblyConfig<T>(input.config(), &matrix_);
    add_incident(input);
  }
  Solution(AssemblyConfig<T>* config, const InciCPtrs<T>& incident,
           const Matrix& matrix, const std::string& file = std::string())
      : matrix_(matrix),
        incident_(incident),
        config_(config),
        input_file_(file),
        original_(false) {}

  virtual ~Solution() {
    if (original_) {
      delete config_;
      for (auto& i : incident_) delete i;
    }
  }

  const Solution& Solve();
  const Inhomo<T>* InWhich(const CS* objCS) const;
  T Resultant(const CS* objCS, const Inhomo<T>* inhomo) const;
  T Resultant(const CS* objCS) const;

  const AssemblyConfig<T>& Config() const { return *config_; }
  const InciCPtrs<T>& Incident() const { return incident_; }
  const InhomoCPtrs<T>& inhomo() const { return Config().inhomo(); }
  const Inhomo<T>* inhomo(size_t sn) const;
  size_t NumInhomo() const { return inhomo().size(); }
  const std::string& InputFN() const { return input_file_; }
  const Matrix* Matrix() const { return &matrix_; }
  double Frequency() const { return matrix_.Frequency(); }
  void ReadCoeff(std::istream& is) {
    Eigen::VectorXcd tmp(config_->NumCoeff());
    for (long i = 0; i < tmp.size(); i++) is >> tmp(i);
    config_->dist_solution(tmp);
    solved_ = true;
  }

 protected:
  bool solved_{false};
  const class Matrix matrix_;
  SolveMethod method_{DFT};
  InciCPtrs<T> incident_;

  // Only the configuration of the "root" assembly is needed.
  // The instantiation of the "root" assembly is not necessary.
  AssemblyConfig<T>* config_;
  std::string input_file_;
  bool original_;

  void add_incident(const input::Solution& input);
};

typedef Solution<AP> SolutionAP;
typedef Solution<IP> SolutionIP;

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T>
const Solution<T>& Solution<T>::Solve() {
  if (solved_) return *this;
  config_->Solve(incident_, method_);
  solved_ = true;
  return *this;
}

template <typename T>
const Inhomo<T>* Solution<T>::inhomo(size_t sn) const {
  return Config().inhomo(sn);
}

template <typename T>
void Solution<T>::add_incident(const input::Solution& input) {
  IncidentGen<T> f(matrix_);
  for (auto& i : input.incident()) incident_.push_back(f(i));
}

template <typename T>
const Inhomo<T>* Solution<T>::InWhich(const CS* objCS) const {
  return config_->InWhich(objCS);
}

template <typename T>
T Solution<T>::Resultant(const CS* objCS, const Inhomo<T>* in) const {
  return config_->Resultant(objCS, in, incident_);
}

template <typename T>
T Solution<T>::Resultant(const CS* objCS) const {
  return config_->Resultant(objCS, InWhich(objCS), incident_);
}

}  // namespace mss

#endif
