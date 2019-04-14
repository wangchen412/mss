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

#ifndef MSS_CIRCLE_H
#define MSS_CIRCLE_H

#include "PointSet.h"

namespace mss {

namespace post {

template <typename T>
class Circle : public PointSet<T> {
  using PointSet<T>::point_;

 public:
  template <typename S>
  Circle(const S* solution, const PosiVect& center, double R, size_t N,
         const std::string& id = "1")
      : PointSet<T>("Circle_" + id), center_(center), R_(R) {
    // Add points:
    for (size_t i = 0; i < N; i++)
      point_.push_back(new Point<T>(
          solution, center + PosiVect(R, i * pi2 / N).Cartesian(),
          i * pi2 / N));
  }

  std::string Shape() const override { return "Circle"; }
  std::ostream& PrintParam(std::ostream& os) const override {
    return os << center_ << "\t\t" << R_ << "\t\t" << point_.size()
              << std::endl;
  }

 private:
  const PosiVect center_;
  const double R_;
};

typedef Circle<IP> CircleIP;
typedef Circle<AP> CircleAP;

}  // namespace post

}  // namespace mss

#endif
