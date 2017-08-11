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

#ifndef MSS_CONTINUITY_H
#define MSS_CONTINUITY_H

#include "../../core/Solution.h"
#include "../geometry/Circle.h"

namespace mss {

namespace post {

template <typename T>
class ContCheck {
 public:
  ContCheck(const Solution<T>* solution, const size_t& np, const double& gap,
            const double& tolerance)
      : sol_(solution),
        P_(np),
        gap_(gap / 2),
        tol_(tolerance),
        mis_(T::NumBv, np) {
    for (auto& i : solution->Incident()) norm_ += i->Norm();
  }
  virtual ~ContCheck() {}

  void ComputeMismatch();

  // If any element in the mismatch matrix exceeds the tolerance, return true
  // for "not continuous".
  bool NC() const { return (mis_.array().abs() > tol_).any(); }

  std::ostream& Print(std::ostream& os) const;

 protected:
  const Solution<T>* sol_;
  const size_t P_;    // The number of points in one set.
  const double gap_;  // The half gap between two sets of points.
  const double tol_;  // The error tolerance.
  T norm_;            // The normalization factor of incident waves.

  // Two sets of states of the points along the interface, one set in each
  // side. The points should have the same sequence.
  virtual void add_state() = 0;
  std::vector<const T*> inner_, outer_;
  Eigen::MatrixXcd mis_;
};

// Continuity Check class for Fiber.
template <typename T>
class CC_Fiber : public ContCheck<T> {
  using cc = ContCheck<T>;

 public:
  CC_Fiber(const Solution<T>* solution, const Fiber<T>* fiber,
           const size_t& np = 42, const double& gap = epsilon,
           const double& tolerance = 1e-4)
      : ContCheck<T>(solution, np, gap, tolerance), f_(fiber) {
    add_circ();
    add_state();
    cc::ComputeMismatch();
  }

  ~CC_Fiber() {
    for (auto& i : circ_) delete i;
  }

 private:
  const Fiber<T>* f_;
  std::vector<post::Circle<T>*> circ_;
  void add_circ();
  void add_state() override;
};

// Continuity Check class for Solution.
template <typename T>
class CC_Solution {
 public:
  CC_Solution(const Solution<T>* solution, const size_t& np = 42,
              const double& gap = epsilon, const double& tolerance = 1e-4) {
    add_cc(solution, np, gap, tolerance);
  }
  ~CC_Solution() {
    for (auto& i : check_) delete i;
  }

  bool isCont() const { return nc_.empty(); }
  std::string WriteAll() const;
  std::string WriteNC() const;
  void Write() const;

 private:
  std::vector<ContCheck<T>*> check_;
  std::vector<const ContCheck<T>*> nc_;
  void add_cc(const Solution<T>*, const size_t&, const double&,
              const double&);
};

// ---------------------------------------------------------------------------
// Inline functions:

// ----------------------------------------
// ContCheck methods:
template <typename T>
inline void ContCheck<T>::ComputeMismatch() {
  for (size_t i = 0; i < inner_.size(); i++)
    mis_.col(i) = ((*inner_[i] - *outer_[i]) / norm_).BV();
}
template <typename T>
inline std::ostream& ContCheck<T>::Print(std::ostream& os) const {
  return os << separator("-") << mis_.transpose() << std::endl;
}

// ----------------------------------------
// CC_Fiber methods:
template <typename T>
inline void CC_Fiber<T>::add_circ() {
  circ_.push_back(new post::Circle<T>(cc::sol_, f_->PositionGLB(),
                                      f_->Radius() - cc::gap_, cc::P_));
  circ_.push_back(new post::Circle<T>(cc::sol_, f_->PositionGLB(),
                                      f_->Radius() + cc::gap_, cc::P_));
}
template <typename T>
inline void CC_Fiber<T>::add_state() {
  for (auto& i : circ_[0]->Points()) {
    if (i->In() != f_)
      exit_error_msg({"Continuity check error, inside points out."});
    this->inner_.push_back(&i->State());
  }
  for (auto& i : circ_[1]->Points()) {
    if (i->In() != nullptr)
      exit_error_msg({"Continuity check error, outside points in."});
    this->outer_.push_back(&i->State());
  }
}

// ----------------------------------------
// CC_Solution methods:
template <typename T>
inline void CC_Solution<T>::add_cc(const Solution<T>* solution,
                                   const size_t& np, const double& gap,
                                   const double& tol) {
  for (auto& i : solution->inhomo()) {
    switch (i->Type()) {
      case FIBER:
        check_.push_back(new CC_Fiber<T>(
            solution, dynamic_cast<const Fiber<T>*>(i), np, gap, tol));
        if (check_.back()->NC()) nc_.push_back(check_.back());
        break;

      default:
        exit_error_msg({"Continuity Check input error."});
    }
  }
}
template <typename T>
inline std::string CC_Solution<T>::WriteAll() const {
  std::string fn = NewFileName("Continuity_ALL", ".dat");
  std::ofstream file(fn);
  for (auto& i : check_) i->Print(file);
  file.close();
  return fn;
}
template <typename T>
inline std::string CC_Solution<T>::WriteNC() const {
  assert(!nc_.empty());
  std::string fn = NewFileName("Continuity_NC", ".dat");
  std::ofstream file(fn);
  for (auto& i : nc_) i->Print(file);
  file.close();
  return fn;
}
template <typename T>
inline void CC_Solution<T>::Write() const {
  if (isCont())
    print_msg({"Continuity of all interfaces checked."});
  else {
    print_error_msg({"One / some of the interfaces discontinuous."});
    exit_error_msg({"The mismatch data is written in ", WriteNC()});
  }
}

}  // namespace post

}  // namespace mss

#endif
