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

#include "../incident/Incident.h"

namespace mss {

template <typename T>
class Inhomogeneity {
 public:
  explicit Inhomogeneity(const PosiVect& position, const double& angle = 0)
      : localCS_(position, angle) {}
  virtual ~Inhomogeneity() {}

  // Resultant states.
  virtual T Scatter(const CS* objCS) const;
  virtual T Inner(const CS* objCS) const;

  // Check if the position of the objective CS is inside the inhomogeneity.
  virtual bool Contain(const CS* objCS) const = 0;

  // The nth Modes. The n should be the serial number, instead of the order.
  // The transformation from the serial number to the order should be done in
  // the derived class.
  virtual T ScatterMode(const CS* objCS, const size_t& sn) const = 0;
  virtual T InnerMode(const CS* objCS, const size_t& sn) const = 0;

  virtual size_t NoN() const = 0;
  virtual size_t NoC() const = 0;
  virtual size_t NoE() const = 0;

  // The matrix composed of effects of modes (correspond to its unknown
  // coefficients) of the source on this inhomogeneity's nodes.
  virtual Eigen::MatrixXcd ModeMatrix(const Inhomogeneity* source) const = 0;

  virtual void SetCoeff(const Eigen::VectorXcd&) = 0;

  // This inhomogeneity's nodes.
  virtual const std::vector<CS>& Node() const = 0;

  const CS* LocalCS() const { return &localCS_; }

 private:
  const CS localCS_;
};

}  // namespace mss

#endif
