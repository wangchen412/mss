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

int main(int argc, char** argv) {
  if (argc != 6) exit_error_msg({"Ka and initial values needed."});

  double omega = 16576.243191120248 * atof(argv[1]);
  Eigen::VectorXcd w(1200), t(1200);
  read_bv("bv.dat", w, t);

  Mismatch f(omega, w, t, {{11400, 11400}, 0, {84e9, 84e9}});
  std::cout << f({atof(argv[2]), atof(argv[3]), atof(argv[4]), atof(argv[5])})
            << std::endl;

  return 0;
}