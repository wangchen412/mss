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
  if (argc != 5) exit_error_msg({"Initial values needed."});
  double omega = 32323.674222684484;

  Solution<AP> s{input::Solution("input.txt")};
  s.Solve();
  Boundary<AP, 4> b{500, {{-1.8, 0.3}, {-0.6, -0.3}}, s.Matrix()};
  Eigen::VectorXcd w(b.NumNode()), t(b.NumNode());
#ifdef NDEBUG
#pragma omp parallel for
#endif
  for (size_t i = 0; i < b.NumNode(); i++) {
    VectorXcd bv = s.Resultant(b.Node(i)).Bv();
    w(i) = bv(0);
    t(i) = bv(1);
  }

  Mismatch f(omega, w, t, {{11400, 11400}, 0, {84e9, 84e9}});
  std::ofstream file("iterations.dat");
  NelderMead(f, &file,
             Eigen::Vector4d(atof(argv[1]), atof(argv[2]), atof(argv[3]),
                             atof(argv[4])),
             1e-3);
  file.close();

  return 0;
}
