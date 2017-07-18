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

#include "Point.h"

namespace mss {

namespace output {

template <typename T>
class Line : public Geometry<T> {
 public:
  Line(const PosiVect& p1, const PosiVect& p2, const size_t& N,
       const Solution<T>* solution, const std::string& id = "1")
      : solution_(solution), id_(id) {
    // Add points:
    PosiVect d = (p2 - p1) / N;
    for (size_t i = 0; i < N; i++)
      point_.push_back(new Point<T>(p1 + d * i, solution_));
  }
  virtual ~Line() {
    // Delete points:
    for (auto& i : point_) delete i;
  }

  void Write() const;

 private:
  std::vector<Point<T>*> point_;
  const Solution<T>* solution_;
  const std::string id_;
};

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T>
void Line<T>::Write() const {
  std::string fileName = std::string("line_") + id_ + std::string(".dat");
  std::ofstream file(fileName);
  file.precision(17);
  for (auto& i : point_) file << *i << std::endl;
  file.close();
}

}  // namespace output

}  // namespace mss

#endif
