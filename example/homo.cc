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
  double omega = 16576.24319;

  Solution<AP> s{input::Solution("input.txt")};
  s.Solve();

  // 4 x 4 from the second.
  Boundary<AP, 4> b{500, {{-0.3, 0.3}, {0.3, -0.3}}, s.Matrix()};
  Eigen::VectorXcd w(b.NumNode()), t(b.NumNode());

  std::vector<StateAP> v(b.Node().size());
#ifdef NDEBUG
#pragma omp parallel for
#endif
  for (size_t i = 0; i < b.NumNode(); i++) {
    v[i] = s.Resultant(b.Node(i));
    w(i) = v[i].Bv()(0);
    t(i) = v[i].Bv()(1);
  }

  std::ofstream bv_out("bv.dat");
  for (size_t i = 0; i < b.NumNode(); i++)
    bv_out << setMaxPrecision << v[i].Basis()->PositionGLB() << "\t"
           << v[i].Basis()->AngleGLB() << "\t" << w(i) << "\t" << t(i)
           << std::endl;
  bv_out.close();

  Mismatch f(omega, w, t, {{11400, 11400}, 0, {84e9, 84e9}});
  std::ofstream file("iterations.dat");
  NelderMead(f, Eigen::Vector4d::Ones(), &file);
  file.close();

  return 0;
}
