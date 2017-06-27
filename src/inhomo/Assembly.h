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

#ifndef MSS_ASSEMBLY_H
#define MSS_ASSEMBLY_H

#include "Inhomogeneity.h"

namespace mss {

template <typename T>
class Assembly : public Inhomogeneity<T> {
 public:
  explicit Assembly(const CS& localCS, const Matrix& matrix)
      : Inhomogeneity<T>(localCS, matrix) {}

  T Scatter(const CS* objCS) const override;
  T Inner(const CS* objCS) const override;

  T ScatterModeP(const CS* objCS, int n) const override;
  T ScatterModeS(const CS* objCS, int n) const override;
  T InnerModeP(const CS* local, int n) const override;
  T InnerModeS(const CS* local, int n) const override;

  Eigen::VectorXcd InVector(
      const std::vector<Incident<T>*>& incident) const override;
  Eigen::MatrixXcd ModeMatrix(const Inhomogeneity<T>* other) const override;

 private:
  std::vector<Inhomogeneity<T>*> inhomo_;
  std::vector<Configuration<T>*> config_;
};

}  // namespace mss

#endif
