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

Material rubber{1300, 1.41908e9, 0.832e9}, lead{11400, 36.32496e9, 8.43e9};
Matrix m{rubber, 1e6}, ff{lead, 1e6};

IncidentPlaneSH* in;
FiberConfig<AP>* fc;
Fiber<AP>* f;
Boundary<AP, 14>* b0;
Boundary<AP, 14>* b1;

StateAP rst1(const CS* cs) {
  if (f->Contains(cs)) return f->Inner(cs);
  return in->Effect(cs) + f->Scatter(cs);
}

StateAP rst2(const CS* cs) {
  if (f->Contains(cs)) return b1->Effect(cs);
  return in->Effect(cs) + b0->Effect(cs);
}

int main() {
  in = new IncidentPlaneSH(m);
  fc = new FiberConfig<AP>("1", 20, 1413, 3e-3, lead, &m);
  f = new Fiber<AP>(fc);
  f->SetCoeff(f->DSolve(in->EffectBv(f->Node())));

  b0 = new Boundary<AP, 14>(50 * m.KT(), {{0, 0}, {0, 3e-3}}, &m, CIRCULAR,
                            true);
  b1 = new Boundary<AP, 14>(50 * m.KT(), {{0, 0}, {0, 3e-3}}, &ff, CIRCULAR);
  Eigen::VectorXcd bv =
      b1->DispToEffect() * (b0->MatrixH() + b0->MatrixG() * b1->DtN())
                               .lu()
                               .solve(in->EffectDv(b0->Node()));
  b1->SetBv(bv);
  for (long i = 0; i < bv.size(); i++)
    if (i % 2) bv(i) *= -1;
  b0->SetBv(bv);

  post::Area<AP> a1(rst1, {-3e-2, 3e-2}, {3e-2, -3e-2}, 300, 300, "1");
  post::Area<AP> a2(rst2, {-3e-2, 3e-2}, {3e-2, -3e-2}, 300, 300, "2");
  a1.Write();
  a2.Write();

  delete in;
  delete fc;
  delete f;
  delete b0;
  delete b1;
  return 0;
}
