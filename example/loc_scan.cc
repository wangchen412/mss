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

int main(int argc, char** argv) {
  if (argc != 2) exit_error_msg({"X of the center is needed."});
  VectorXcd w, t;
  double omega = 32323.674222684484;
  double x1 = atof(argv[1]) * 0.2 - 2.2;
  compute_bv(omega, w, t, x1, 0.4, x1 + 0.8, -0.4);
  Mismatch f(omega, w, t, {{11400, 11400}, 0, {84e9, 84e9}}, 0.8, 0.8);
  std::ofstream file("iterations.dat");
  VectorXd p = NelderMead(f, Vector4d::Ones(), &file);
  file.close();
  std::cout << mss_msg({"Effective properties: "}) << p.transpose()
            << std::endl;

  Material steel(7670, 116e9, 84.3e9);
  mss::Matrix m{steel, omega}, ff{f.material(p), omega};
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
