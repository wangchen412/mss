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

#ifndef MSS_POINTSET_H
#define MSS_POINTSET_H

#include "Point.h"

namespace mss {

namespace post {

template <typename T>
class PointSet : public Geometry<T> {
 public:
  PointSet(const Solution<T>* solution, const std::string& id)
      : Geometry<T>(solution, id) {}
  virtual ~PointSet() {
    for (auto& i : point_) delete i;
    point_.clear();
  }

  const PtCPtrs<T>& Points() const { return point_; }
  std::ostream& Print(std::ostream& os) const override {
    for (auto& i : point_) i->Print(os);
    return os;
  }

 protected:
  PtCPtrs<T> point_;
};

}  // namespace post

}  // namespace mss

#endif