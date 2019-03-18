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

#include "mismatch.h"

using namespace mss;

int main() {
  double omega = 32323.674222684484;
  Eigen::VectorXcd w, t;
  read_bv("bv.dat", w, t);
  Mismatch f(omega, w, t, {{11400, 11400}, 0, {84e9, 84e9}}, 0.8, 0.8);

  while (1) {
    double rr, ri, mr, mi;
    std::cin >> rr >> ri >> mr >> mi;
    std::cout << setMaxPrecision << f({rr, ri, mr, mi}) << std::endl;
    std::cout.flush();
  }

  return 0;
}
