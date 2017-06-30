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

// Material struct, in which the material constants are recorded.

#ifndef MSS_MATERIAL_H
#define MSS_MATERIAL_H

#include "Input.h"
#include "State.h"

namespace mss {

struct Material {
  explicit Material(const double& rho, const double& lambda, const double& mu)
      : rho(rho), lambda(lambda), mu(mu) {
    _computeC();
  }

  explicit Material(const input::Material& input)
      : rho(input.rho), lambda(input.lambda), mu(input.mu) {
    _computeC();
  }

  StressAP C(const StrainAP& gamma) const { return gamma * mu; }
  StressIP C(const StrainIP& gamma) const {
    dcomp gkk = gamma.xx + gamma.yy;
    return StressIP(lambda * gkk + 2 * mu * gamma.xx,
                    lambda * gkk + 2 * mu * gamma.yy, mu * gamma.xy);
  }

  double rho, lambda, mu, cl, ct;

 private:
  void _computeC() {
    assert(rho > 0 && lambda > 0 && mu > 0);
    cl = std::sqrt((lambda + 2 * mu) / rho);
    ct = std::sqrt(mu / rho);
  }
};

}  // namespace mss

#endif
