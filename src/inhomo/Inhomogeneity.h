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
#include "Configuration.h"

namespace mss {

template <typename T>
class Inhomogeneity {
 public:
  explicit Inhomogeneity(const Configuration<T>* config,
                         const PosiVect& position, const double& angle = 0)
      : localCS_(position, angle),
        kl_m(config->KL_m()),
        kt_m(config->KT_m()) {}
  virtual ~Inhomogeneity() {}

  virtual T Scatter(const CS* objCS) const;
  virtual T Inner(const CS* objCS) const;
  virtual T ScatterModeL(const CS* objCS, int n) const = 0;
  virtual T ScatterModeT(const CS* objCS, int n) const = 0;
  virtual T InnerModeL(const CS* objCS, int n) const = 0;
  virtual T InnerModeT(const CS* objCS, int n) const = 0;

  virtual const Configuration<T>* Config() const = 0;

  virtual int NoP() const { return Config()->NoP(); }
  virtual int NoC() const { return Config()->NoC(); }
  virtual int NoE() const { return Config()->NoE(); }

  virtual const Eigen::MatrixXcd& TransMatrix() {
    return Config()->TransMatrix();
  }
  virtual Eigen::MatrixXcd ModeMatrix(const Inhomogeneity* other) const = 0;
  virtual Eigen::VectorXcd InVector(
      const std::vector<Incident<T>*>& incident) const = 0;

  const CS* LocalCS() const { return &localCS_; }

 protected:
  const CS localCS_;
  const double &kl_m, &kt_m;
};

}  // namespace mss

#endif
