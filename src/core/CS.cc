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

#include "CS.h"

namespace mss {

CS CS::in(const CS* other) const {
  assert(other != this);
  if (other == basis_) return *this;
  Vector<double> r;
  double a;
  if (other) {
    CS sub(this->in(nullptr)), sup(other->in(nullptr));
    r = (sub.position_ - sup.position_).RotateInPlace(sup.angle_);
    a = sub.angle_ - sup.angle_;
  } else {
    CS basis(basis_->in(nullptr));
    r = position_.Rotate(-basis.angle_) + basis.position_;
    a = angle_ + basis.angle_;
  }
  return CS(r, a, other);
}

}  // namespace mss
