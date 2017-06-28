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

#include "Fiber.h"

namespace mss {

template <>
StateAP Fiber<StateAP>::ScatterModeS(const mss::CS* objCS, int n) const {
  CS cs(objCS->in(LocalCS()));
  const double& r = cs.Position().Polar().x;
  const double& t = cs.Position().Polar().y;

  dcomp exp_int = exp(n * t * ii);
  dcomp Hn_kr = Hankel(n, r * kt_m);
  dcomp Hnd_kr = Hankel_dv(n, r * kt_m);

  dcomp w = Hn_kr * exp_int;
  dcomp tzr = mu_m * kt_m * Hnd_kr * exp_int;
  dcomp tzt = n * mu_m / r * ii * Hn_kr * exp_int;

  dcomp norm = Hankel(n, kt_m * config_->CharLength());
  return StateAP(w, tzr, tzt, &cs).in(objCS) / norm;
}
template <>
StateAP Fiber<StateAP>::InnerModeS(const mss::CS* objCS, int n) const {
  CS cs(objCS->in(LocalCS()));
  const double& r = cs.Position().Polar().x;
  const double& t = cs.Position().Polar().y;

  dcomp exp_int = exp(n * t * ii);
  dcomp Jn_kr = Bessel(n, r * kt_m);
  dcomp Jnd_kr = Bessel_dv(n, r * kt_m);

  dcomp w = Jn_kr * exp_int;

}

}  // namespace mss
