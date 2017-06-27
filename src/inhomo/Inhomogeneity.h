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

#ifndef MSS_INHOMOGENEITY_H
#define MSS_INHOMOGENEITY_H

#include "../core/Incident.h"
#include "Configuration.h"

namespace mss {

template <typename T>
class Inhomogeneity {
 public:
  explicit Inhomogeneity(const CS& localCS, const Matrix& matrix)
      : localCS_(localCS),
        matrix_(matrix),
        kl_m(matrix.KL()),
        kt_m(matrix.KT()),
        l_m(matrix.Lambda()),
        mu_m(matrix.Mu()) {}

  virtual T Scatter(const CS* objCS) const = 0;
  virtual T Inner(const CS* objCS) const = 0;

  virtual T ScatterModeP(const CS* objCS, int n) const = 0;
  virtual T ScatterModeS(const CS* objCS, int n) const = 0;
  virtual T InnerModeP(const CS* objCS, int n) const = 0;
  virtual T InnerModeS(const CS* objCS, int n) const = 0;

  virtual Eigen::VectorXcd InVector(
      const std::vector<Incident<T>*>& incident) const = 0;
  virtual Eigen::MatrixXcd ModeMatrix(const Inhomogeneity* other) const = 0;

  const CS* LocalCS() const { return &localCS_; }
  const Matrix& Matrix() const { return matrix_; }
  //virtual const Configuration<T>& Config() const = 0;

 protected:
  const CS localCS_;
  const class Matrix& matrix_;
  const double &kl_m, &kt_m, &l_m, &mu_m;

};

}  // namespace mss

#endif
