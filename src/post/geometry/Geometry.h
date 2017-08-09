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

#ifndef MSS_GEOMETRY_H
#define MSS_GEOMETRY_H

#include "../../core/Solution.h"
#include "../../core/State.h"
#include "../../inhomo/Inhomo.h"
#include "../../tools/FileIO.h"

namespace mss {

namespace post {

template <typename T>
class Geometry {
 public:
  Geometry(const Solution<T>* solution, const std::string& id)
      : solution_(solution), id_(id) {}
  virtual ~Geometry() {}
  virtual std::ostream& Print(std::ostream& os) const = 0;
  void Write() const;

  const std::string& ID() const { return id_; }

 protected:
  const Solution<T>* solution_;
  const std::string id_;
};

template <typename T>
using GeoPtrs = std::vector<Geometry<T>*>;

template <typename T>
using GeoCPtrs = std::vector<const Geometry<T>*>;

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T>
void Geometry<T>::Write() const {
  std::ofstream file(this->id_ + ".dat");
  Print(file);
  file.close();
}

}  // namespace post

}  // namespace mss

#endif
