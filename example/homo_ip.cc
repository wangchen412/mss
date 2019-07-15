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
  VectorXcd wl(400), wt(400), tl(400), tt(400);

  std::ifstream file("ip.dat");
  std::string tmp;
  for (int i = 0; i < 400; i++) {
    std::getline(file, tmp);
    std::stringstream ss(tmp);
    double a;
    for (int j = 0; j < 12; j++) ss >> a;
    ss >> wl(i);
    ss >> wt(i);
    ss >> tl(i);
    ss >> tt(i);
  }

  Mismatch_IP f(10000, wl, tl, wt, tt, {{11400, 0}, {116e9, 0}, {84e9, 0}},
                0.2, 0.2);

  Eigen::VectorXd x0(3);
  // x0 << 1, 0.1, 1, 0.1, 1, 0.1;
  x0 << 0.5, 0.5, 0.5;
  NelderMead(f, x0, &std::cout);

  return 0;
}
