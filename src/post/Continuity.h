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

#include "../core/Solution.h"
#include "output/Circle.h"

namespace mss {

template <typename T>
class ContinuityCheck {
 public:
  virtual ~ContinuityCheck() {}

 protected:
  const double gap_{epsilon};
  const size_t P_{823};
};

// Continuity Check class for Fiber.
template <typename T>
class CcFiber : public ContinuityCheck<T> {
 public:
 private:
  output::Circle<T> inner, outer;
};

// template <typename T>
// class ContinuityCheck {
//  public:
//   ContinuityCheck(const std::vector<Inhomogeneity<T>*>& inhomo,
//                   const std::vector<Incident<T>*>& incident,
//                   const size_t& nop = 823, const double& gap = epsilon)
//     : inhomo_(inhomo), incident_(incident), nop_(nop), gap_(gap) {}
//   virtual ~ContinuityCheck() {}

//  private:
//   const double gap_;
//   const size_t nop_;
//   const std::vector<Inhomogeneity<T>*>& inhomo_;
//   const std::vector<Incident<T>*>& incident_;
// };

}  // namespace mss

#endif
