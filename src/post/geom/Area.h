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

#ifndef MSS_AREA_H
#define MSS_AREA_H

#include "PointSet.h"

namespace mss {

namespace post {

template <typename T>
class Area : public PointSet<T> {
  using PointSet<T>::point_;

 public:
  Area(const Solution<T>* solution, const PosiVect& p1, const PosiVect& p2,
       const size_t& Nx, const size_t& Ny, const std::string& id = "1")
      : PointSet<T>(solution, "Area_" + id) {
    // Add points:
    PosiVect dx((p2 - p1).x / Nx, 0);
    PosiVect dy(0, (p2 - p1).y / Ny);
    for (size_t j = 0; j < Ny; j++)
      for (size_t i = 0; i < Nx; i++)
        point_.push_back(new Point<T>(solution, p1 + dx * i + dy * j));
  }
};

typedef Area<StateIP> AreaIP;
typedef Area<StateAP> AreaAP;

}  // namespace post

}  // namespace mss

#endif
