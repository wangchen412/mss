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
using namespace Eigen;

IncidentPlaneSH* in;
Boundary<AP, 14>* b0;
Boundary<AP, 14>* b1;

StateAP rst(const CS* cs) {
  switch (b1->Contains(cs)) {
  case 1:
    return b1->Effect(cs);
  case -1:
    return in->Effect(cs) + b0->Effect(cs);
  }
  return b1->Effect(cs, true);
}

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
  VectorXd p = NelderMead(f, Vector4d::Ones(), &file);
  file.close();
  std::cout << mss_msg({"Effective properties: "}) << p.transpose()
            << std::endl;

  Material steel(7670, 116e9, 84.3e9);
  mss::Matrix m{steel, s.Frequency()}, ff{f.material(p), s.Frequency()};
  in = new IncidentPlaneSH(m);
  b0 = new Boundary<AP, 14>(50 * m.KT(), {{-2, 2}, {2, -2}}, &m, RECTANGULAR,
                            true);
  b1 = new Boundary<AP, 14>(50 * m.KT(), {{-2, 2}, {2, -2}}, &ff);

  Eigen::VectorXcd bv =
    b1->DispToEffect() * (b0->MatrixH() + b0->MatrixG() * b1->DtN())
    .lu()
    .solve(in->EffectDv(b0->Node()));
  b1->SetBv(bv);
  for (long i = 0; i < bv.size(); i++)
    if (i % 2) bv(i) *= -1;
  b0->SetBv(bv);

  post::Area<AP> a2(rst, {-6, 6}, {6, -6}, 1200, 1200, "bem");
  a2.Write();

  delete in;
  delete b0;
  delete b1;
  return 0;
}
