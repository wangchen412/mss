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

#ifndef MSS_FIBER_H
#define MSS_FIBER_H

#include "ConfigFiber.h"
#include "Inhomogeneity.h"

namespace mss {

template <typename T>
class Fiber : public Inhomogeneity<T> {
 public:
  explicit Fiber(const CS& localCS, const ConfigFiber<T>* config)
      : Inhomogeneity<T>(localCS, config),
        config_(config),
        kl_f(config->KL()),
        kt_f(config->KT()) {}

  T ScatterModeL(const CS* objCS, int n) const override;
  T ScatterModeT(const CS* objCS, int n) const override;
  T InnerModeL(const CS* objCS, int n) const override;
  T InnerModeT(const CS* objCS, int n) const override;

  Eigen::VectorXcd InVector(
      const std::vector<Incident<T>*>& incident) const override;
  Eigen::MatrixXcd ModeMatrix(const Inhomogeneity<T>* other) const override;

  const ConfigFiber<T>* Config() const { return config_; }

 protected:
  const ConfigFiber<T>* config_;
  const double &kl_f, &kt_f;

  T _modeT(const CS* objCS, EigenFunctor& f, const Material& mat) const;
  T _modeL(const CS* objCS, EigenFunctor& f, const Material& mat) const;
};

}  // namespace mss

#endif
