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

// Virtual base class of inhomogeneity classes.
// Contains local CS and informations about the matrix.

#ifndef MSS_INHOMOGENEITY_H
#define MSS_INHOMOGENEITY_H

#include "../core/Incident.h"

namespace mss {

template <typename T>
class Inhomogeneity {
 public:
  explicit Inhomogeneity(const PosiVect& position, const double& angle = 0)
      : localCS_(position, angle) {}
  virtual ~Inhomogeneity() {}

  virtual T Scatter(const CS* objCS) const;
  virtual T Inner(const CS* objCS) const;

  // The nth Modes. The n should be the serial number, instead of the order.
  // The transformation from the serial number to the order should be done in
  // the derived class.
  virtual T ScatterModeL(const CS* objCS, const size_t& n) const = 0;
  virtual T ScatterModeT(const CS* objCS, const size_t& n) const = 0;
  virtual T InnerModeL(const CS* objCS, const size_t& n) const = 0;
  virtual T InnerModeT(const CS* objCS, const size_t& n) const = 0;

  virtual size_t NoN() const = 0;
  virtual size_t NoC() const = 0;
  virtual size_t NoE() const = 0;

  virtual const Eigen::MatrixXcd& TransMatrix() const = 0;
  virtual Eigen::MatrixXcd ModeMatrix(const Inhomogeneity* other) const = 0;
  virtual Eigen::VectorXcd InVector(
      const std::vector<Incident<T>*>& incident) const = 0;

  const CS* LocalCS() const { return &localCS_; }

 protected:
  const CS localCS_;
};

}  // namespace mss

#endif
