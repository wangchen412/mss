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

// Coordinate system class.

#ifndef MSS_CS_H
#define MSS_CS_H

#include "Tensor.h"

namespace mss {

class CS {
 public:
  explicit CS(const double& x = 0, const double& y = 0,
              const double& angle = 0, const CS* basis = nullptr)
      : position_(x, y), angle_(angle), basis_(basis) {}
  explicit CS(const PosiVect& position, const double& angle = 0,
              const CS* basis = nullptr)
      : position_(position), angle_(angle), basis_(basis) {}
  CS(const CS& other)
      : position_(other.position_),
        angle_(other.angle_),
        basis_(other.basis_) {}

  // For instantiation:
  CS(const CS& other, const CS* basis)
      : position_(other.position_), angle_(other.angle_), basis_(basis) {}
  virtual ~CS() {}

  bool operator==(const CS& other) const;

  // Output the 
  friend std::ostream& operator<<(std::ostream& os, const CS& cs) {
    return os << cs.PositionGLB() << "\t" << cs.AngleGLB();
  }

  // Return a new CS object which is the map of this CS in another basis CS.
  CS in(const CS* otherBasis) const;

  // Return a new CS object which is in the global CS.
  CS inGLB() const { return in(nullptr); }

  const CS* Basis() const { return basis_; }
  const PosiVect& Position() const { return position_; }
  PosiVect PositionGLB() const { return inGLB().Position(); }
  PosiVect PositionIn(const CS* otherBasis) const;
  const double& Angle() const { return angle_; }
  double AngleGLB() const { return inGLB().Angle(); }

 private:
  PosiVect position_;
  double angle_;
  const CS* basis_;
};

// ---------------------------------------------------------------------------
// Inline functions:

inline bool CS::operator==(const CS& other) const {
  return (position_ == other.position_) && angEqu(angle_, other.angle_) &&
         (basis_ == other.basis_);
}
inline PosiVect CS::PositionIn(const CS* otherBasis) const {
  return in(otherBasis).Position();
}

}  // namespace mss

#endif
