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

#include "../src/post/Output.h"
#include "../src/post/check/Continuity.h"
#include "mismatch.h"

using namespace mss;

int main(int argc, char* argv[]) {
  if (argc != 2) exit_error_msg({"Input required."});

  Solution<AP> s{input::Solution(argv[1])};
  Boundary<AP, 4> b(0, {}, s.Matrix(), INPUT);

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

  Mismatch f(s.Frequency(), w, t, {{11400, 11400}, 0, {84e9, 84e9}});
  std::ofstream file("iterations.dat");
  NelderMead(f,
             Eigen::Vector4d(atof(argv[1]), atof(argv[2]), atof(argv[3]),
                             atof(argv[4])),
             &file, 1e-3);
  file.close();

  return 0;
}
