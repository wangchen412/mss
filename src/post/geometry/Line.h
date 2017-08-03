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

#ifndef MSS_LINE_H
#define MSS_LINE_H

#include "PointSet.h"

namespace mss {

namespace post {

template <typename T>
class Line : public PointSet<T> {
  using PointSet<T>::point_;

 public:
  Line(const Solution<T>* solution, const PosiVect& p1, const PosiVect& p2,
       const size_t& N, const std::string& id = "1")
      : PointSet<T>(solution, "Line_" + id) {
    // Add points:
    PosiVect d = (p2 - p1) / N;
    for (size_t i = 0; i < N; i++)
      point_.push_back(new Point<T>(solution, p1 + d * i));
  }
};

typedef Line<StateIP> LineIP;
typedef Line<StateAP> LineAP;

}  // namespace post

}  // namespace mss

#endif
