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

#include "../pre/Input.h"
#include "State.h"

namespace mss {

class Material {
 public:
  Material(const dcomp& rho, const dcomp& lambda, const dcomp& mu)
      : rho_(rho), lambda_(lambda), mu_(mu) {
    _computeWaveSpeed();
  }
  Material(const input::Material& input)
      : Material(input.rho, input.lambda, input.mu) {
    assert(input.cl == cl_ && input.ct == ct_);
  }
  Material(const Material& other)
      : Material(other.rho_, other.lambda_, other.mu_) {}

  StressAP C(const dcomp& gzx, const dcomp& gzy) const {
    return StressAP(mu_ * gzx, mu_ * gzy);
  }
  StressIP C(const dcomp& gxx, const dcomp& gyy, const dcomp& gxy) const {
    dcomp lkk = lambda_ * (gxx + gyy);
    return StressIP(lkk + 2 * mu_ * gxx, lkk + 2 * mu_ * gyy, mu_ * gxy);
  }
  StressAP C(const StrainAP& g) const { return C(g.x, g.y); }
  StressIP C(const StrainIP& g) const { return C(g.xx, g.yy, g.xy); }

  bool operator==(const Material& other) const {
    return (rho_ == other.rho_) && (lambda_ == other.lambda_) &&
           (mu_ == other.mu_);
  }

  double MassDensity() const { return rho_.real(); }
  double Lambda() const { return lambda_.real(); }
  double Mu() const { return mu_.real(); }

  dcomp Rho_comp() const { return rho_; }
  dcomp Lambda_comp() const { return lambda_; }
  dcomp Mu_comp() const { return mu_; }
  dcomp CL() const { return cl_; }
  dcomp CT() const { return ct_; }

  Material operator*(const Eigen::Vector4d& r) const {
    return Material({rho_.real() * r(0), rho_.imag() * r(1)}, lambda_,
                    {mu_.real() * r(2), mu_.imag() * r(3)});
  }

 private:
  dcomp rho_, lambda_, mu_;
  dcomp cl_, ct_;

  void _computeWaveSpeed() {
    // assert(rho_ > 0 && lambda_ > 0 && mu_ > 0);
    cl_ = std::sqrt((lambda_ + 2 * mu_) / rho_);
    ct_ = std::sqrt(mu_ / rho_);
  }
};

}  // namespace mss

#endif
