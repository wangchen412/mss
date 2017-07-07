// ----------------------------------------------------------------------
//
// Copyright © 2017 mss authors.
//
// mss is free software: you can redistribute it and/or modify
// it under the terms of the GNU General lPublic License as published
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

#include "ConfigAssembly.h"
#include "Inhomogeneity.h"

namespace mss {

template <typename T>
class Assembly : public Inhomogeneity<T> {
 public:
  explicit Assembly(const ConfigAssembly<T>* config, const PosiVect& position,
                    const double& angle = 0)
      : Inhomogeneity<T>(position, angle), config_(config) {}

  // T Scatter(const CS* objCS) const override;
  // T Inner(const CS* objCS) const override;

  T ScatterModeL(const CS* objCS, int n) const override;
  T ScatterModeT(const CS* objCS, int n) const override;
  T InnerModeL(const CS* local, int n) const override;
  T InnerModeT(const CS* local, int n) const override;

  Eigen::VectorXcd InVector(
      const std::vector<Incident<T>*>& incident) const override;
  Eigen::MatrixXcd ModeMatrix(const Inhomogeneity<T>* other) const override;

 private:
  const ConfigAssembly<T>* config_;
  std::vector<CS> node_;
};

}  // namespace mss

#endif
