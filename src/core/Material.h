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

#include "Math.h"

namespace mss {

struct Material {
  explicit Material(const double& rho, const double& lambda, const double& mu)
      : rho(rho), lambda(lambda), mu(mu) {
    assert(rho > 0 && lambda > 0 && mu > 0);
    cl = std::sqrt((lambda + 2 * mu) / rho);
    ct = std::sqrt(mu / rho);
  }

  double rho, lambda, mu, cl, ct;
};

}  // namespace mss

#endif
