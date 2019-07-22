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
  PointSet(const std::string& id) : Geometry<T>(id) {}
  virtual ~PointSet() {
    for (auto& i : point_) delete i;
    point_.clear();
  }

  const PtCPtrs<T>& Points() const { return point_; }
  const Point<T>* Point(size_t i) const { return point_[i]; }
  std::ostream& Print(std::ostream& os) const override;
  virtual std::ostream& PrintParam(std::ostream& os) const = 0;

  virtual std::string Shape() const = 0;

  std::ostream& PrintHd(std::ostream& os) const;

  double StrainEnergyDensity() const {
    double rst = 0;
    for (auto& i : point_) rst += i->StrainEnergy();
    return rst / point_.size();
  }

  double KineticEnergyDensity() const {
    double rst = 0;
    for (auto& i : point_) rst += i->KineticEnergy();
    return rst / point_.size();
  }

 protected:
  PtCPtrs<T> point_;
};

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T>
std::ostream& PointSet<T>::Print(std::ostream& os) const {
  PrintHd(os);
  for (auto& i : point_) i->Print(os);
  return os;
}

template <typename T>
std::ostream& PointSet<T>::PrintHd(std::ostream& os) const {
  os << separator("=");
  os << "Shape:\t" << Shape() << std::endl;
  os << "Type:\t" << T::Type << std::endl;
  os << "Parameters:\t";
  PrintParam(os);
  os << separator("-");
  return os;
}

}  // namespace post

}  // namespace mss

#endif
