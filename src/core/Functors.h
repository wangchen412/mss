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

#ifndef MSS_FUNCTORS_H
#define MSS_FUNCTORS_H

#include "Tensor.h"

namespace mss {

typedef decltype((Hn)) BesselFunc;

class EigenFunctor {
 public:
  EigenFunctor(BesselFunc f, int n, const double& k)
      : f(f), n(n), ndr(0), ndt(0), k(k) {}

  dcomp operator()(const PosiVect& p) const {
    return pow(k, ndr) * pow(dcomp(0, n), ndt) * d(n, k * p.x, ndr) *
           exp(n * p.y * ii);
  }

  dcomp d(int j, const double& x, int i) const {
    return i > 0 ? (d(j - 1, x, i - 1) - d(j + 1, x, i - 1)) / 2.0 : f(j, x);
  }
  EigenFunctor dr() const {
    EigenFunctor rst(*this);
    rst.ndr++;
    return rst;
  }
  EigenFunctor dt() const {
    EigenFunctor rst(*this);
    rst.ndt++;
    return rst;
  }

 private:
  BesselFunc f;
  int n, ndr, ndt;
  double k;
};

}  // namespace mss

#endif
