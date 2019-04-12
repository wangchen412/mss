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

#include <fstream>
#include <sstream>
#include "../src/core/Solution.h"
#include "../src/post/Output.h"

using namespace mss;

Material steel{7670, 116e9, 84.3e9};
Material norm{{11400, 11400}, 1, {84e9, 84e9}};
Material eff_mat{norm *
                 Eigen::Vector4d(-0.129907, 0.927331, 0.0357903, 0.519424)};
Matrix m{steel, 32323.674222684484}, ff{eff_mat, 32323.674222684484};

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

int main() {
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

  post::Line<AP>(rst, {-6, 0}, {6, 0}, 1200, "h_0").Write();
  post::Line<AP>(rst, {-6, 4}, {6, 4}, 1200, "h_4").Write();
  post::Line<AP>(rst, {0, 6}, {0, -6}, 1200, "v_0").Write();
  post::Line<AP>(rst, {4, 6}, {4, -6}, 1200, "v_4").Write();

  delete in;
  delete b0;
  delete b1;
  return 0;
}
