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

#include "Point.h"

namespace mss {

namespace post {

template <typename T>
class Circle : public Geometry<T> {
 public:
  Circle(const PosiVect& center, const double& R, const size_t& N,
         const Solution<T>* solution, const std::string& id = "1")
      : solution_(solution), id_(id) {
    // Add points:
    for (size_t i = 0; i < N; i++)
      point_.push_back(
          new Point<T>(center + PosiVect(R, i * pi2 / N).Cartesian(),
                       solution, i * pi2 / N));
  }
  virtual ~Circle() {
    // Delete points:
    for (auto& i : point_) delete i;
  }

  const PtCPtrs<T>& Points() const { return point_; }

  std::ostream& Print(std::ostream& os) const;
  void Write() const;

 private:
  PtCPtrs<T> point_;
  const Solution<T>* solution_;
  const std::string id_;
};

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T>
inline std::ostream& Circle<T>::Print(std::ostream& os) const {
  for (auto& i : point_) os << *i << std::endl;
  return os;
}

template <typename T>
inline void Circle<T>::Write() const {
  std::string fileName = std::string("Circle_") + id_ + std::string(".dat");
  std::ofstream file(fileName);
  file.precision(17);
  Print(file);
  file.close();
}

}  // namespace post

}  // namespace mss

#endif
