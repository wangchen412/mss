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
#include "../../inhomo/Inhomogeneity.h"
#include "../../tools/FileIO.h"

namespace mss {

namespace post {

template <typename T>
class Geometry {
 public:
  virtual ~Geometry() {}
  virtual void Write() const = 0;
};

template <typename T>
using GeoPtrs = std::vector<Geometry<T>*>;

template <typename T>
using GeoCPtrs = std::vector<const Geometry<T>*>;

}  // namespace post

}  // namespace mss

#endif
