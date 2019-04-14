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

#include "BEM_Solution.h"

using namespace mss;

int main(int argc, char** argv) {
  if (argc != 6) exit_error_msg({"Effective material properties needed."});
  double omega = atof(argv[1]);
  Material steel{7670, 116e9, 84.3e9};
  Material norm{{11400, 11400}, 1, {84e9, 84e9}};
  Eigen::Vector4d eff{atof(argv[2]), atof(argv[3]), atof(argv[4]),
                      atof(argv[5])};
  Material eff_mat{norm * eff};
  BEM_Solution<AP, 10> s(omega, eff_mat);
  post::Line<AP> l1(&s, {-6, 0}, {6, 0}, 1200);
  l1.Write();
  return 0;
}
