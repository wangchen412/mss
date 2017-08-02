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

#ifndef MSS_POINT_H
#define MSS_POINT_H

#include "Geometry.h"

namespace mss {

namespace post {

template <typename T>
class Point : public Geometry<T> {
 public:
  Point(const Solution<T>* solution, const PosiVect& position,
        const double& angle = 0, const std::string& id = "1")
      : Geometry<T>(solution, "Point_" + id),
        localCS_(position, angle),
        in_(solution->InWhich(&localCS_)),
        state_(solution->Resultant(&localCS_, in_)) {}
  virtual ~Point() {}

  friend std::ostream& operator<<(std::ostream& os, const Point<T>& pt) {
    return os << setMaxPrecision << pt.localCS_ << "\t" << pt.state_
              << std::endl;
  }

  const CS* LocalCS() const { return &localCS_; }
  const T& State() const { return state_; }
  const Inhomogeneity<T>* In() const { return in_; }

  std::ostream& Print(std::ostream& os) const override;

 private:
  const CS localCS_;
  const Inhomogeneity<T>* in_;  // The pointer to the inhomogeneity in which
                                // the point is, if the point is in one.
  const T state_;
};

template <typename T>
using PtPtrs = std::vector<Point<T>*>;
template <typename T>
using PtCPtrs = std::vector<const Point<T>*>;

typedef Point<StateIP> PointIP;
typedef Point<StateAP> PointAP;

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T>
std::ostream& Point<T>::Print(std::ostream& os) const {
  return os << *this;
}

}  // namespace post

}  // namespace mss

#endif
