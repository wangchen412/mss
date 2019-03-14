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

double f(const Eigen::VectorXd& x) {
  std::cout << x.transpose() << std::endl;
  double x1 = x(0), x2 = x(1);
  return 100 * (x2 - x1 * x1) * (x2 - x1 * x1) + (1 - x1) * (1 - x1);
}

int main(int argc, char** argv) {
  // if (argc != 5) exit_error_msg({"Initial values needed."});
  double omega = 16576.243191120248 * 1.2;

  Eigen::VectorXcd w(1200), t(1200);
  // if (argc == 2) compute_bv(omega, w, t);
  // if (argc == 3) read_bv(argv[2], w, t);
  read_bv("bv.dat", w, t);

  // Mismatch f(omega, w, t, {{11400, 11400}, 0, {84e9, 84e9}});
  // std::ofstream file("iterations.dat");
  // GradientDescent(f, &file, {atof(argv[1]), atof(argv[2]), atof(argv[3]),
  // atof(argv[4])});
  // file.close();

  NelderMead(f, nullptr, Eigen::Vector2d(1.2, 1));

  return 0;
}
