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

// Material class, in which the material constants are recorded.
// The constitutive relation can be used as the method C.

#ifndef MSS_MATERIAL_H
#define MSS_MATERIAL_H

#include "../input/Input.h"
#include "State.h"

namespace mss {

class Material {
 public:
  Material(const double& rho, const double& lambda, const double& mu)
      : rho_(rho), lambda_(lambda), mu_(mu) {
    _computeWaveSpeed();
  }
  Material(const input::Material& input)
      : Material(input.rho, input.lambda, input.mu) {}

  StressAP C(const dcomp& gzx, const dcomp& gzy) const {
    return StressAP(mu_ * gzx, mu_ * gzy);
  }
  StressIP C(const dcomp& gxx, const dcomp& gyy, const dcomp& gxy) const {
    dcomp lkk = lambda_ * (gxx + gyy);
    return StressIP(lkk + 2 * mu_ * gxx, lkk + 2 * mu_ * gyy, mu_ * gxy);
  }
  StressAP C(const StrainAP& g) const { return C(g.x, g.y); }
  StressIP C(const StrainIP& g) const { return C(g.xx, g.yy, g.xy); }

  const double& MassDensity() const { return rho_; }
  const double& Lambda() const { return lambda_; }
  const double& Mu() const { return mu_; }
  const double& CL() const { return cl_; }
  const double& CT() const { return ct_; }

 private:
  double rho_, lambda_, mu_, cl_, ct_;

  void _computeWaveSpeed() {
    assert(rho_ > 0 && lambda_ > 0 && mu_ > 0);
    cl_ = std::sqrt((lambda_ + 2 * mu_) / rho_);
    ct_ = std::sqrt(mu_ / rho_);
  }
};

}  // namespace mss

#endif
