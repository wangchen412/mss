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

namespace output {

template <typename T>
class Point : public Geometry<T> {
 public:
  Point(const PosiVect& position, const Solution<T>* solution,
        const double& angle = 0, const std::string& id = "1")
      : localCS_(position, angle),
        state_(&localCS_),
        solution_(solution),
        in_(solution_->InWhich(localCS_)),
        id_(id) {
    state_ = solution_->Resultant(&localCS_, in_);
  }

  virtual ~Point() {}

  friend std::ostream& operator<<(std::ostream& os, const Point<T>& pt) {
    return os << pt.localCS_ << "\t" << pt.state_;
  }

  void Write() const;

 private:
  CS localCS_;
  T state_;
  const Solution<T>* solution_;  // The position in the global CS.
  const Inhomogeneity<T>* in_;   // The pointer to the inhomogeneity in which
                                 // the point is, if the point is in one.
  const std::string id_;
};

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T>
void Point<T>::Write() const {
  std::string fileName = std::string("point_") + id_ + std::string(".dat");
  std::ofstream file(fileName);
  file.precision(17);
  file << *this << std::endl;
  file.close();
}

}  // namespace output

}  // namespace mss

#endif
