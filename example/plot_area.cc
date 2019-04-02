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

const Material steel(7670, 116e9, 84.3e9), lead(11400, 36e9, 8.43e9);
const double a(0.2);
Solution<AP>* solution;

StateAP rst(const CS* cs) {
  return solution->Resultant(cs);
}

double omega(double ka) {
  return ka / a * steel.CT().real();
}
void freq_scan(double ka) {
  Matrix m(steel, omega(ka));
  IncidentPlaneSH in(m, 0, 1e-6);

  auto fc = new FiberConfig<AP>("1", 40, 400, 0.06, lead, &m);
  InhomoPtrs<AP> inhomo;
  int n = 20;
  double c = (n - 1) * a / 2;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      inhomo.push_back(new Fiber<AP>(fc, {i * a - c, j * a - c}));
  auto ac = new AssemblyConfig<AP>("1", inhomo, &m);
  solution = new Solution<AP>(ac, {&in}, m);
  solution->Solve();

  post::Area<AP> a(rst, {-4, 4}, {4, -4}, 400, 400,
                   std::to_string(omega(ka)));
  a.Write();

  delete solution;
  delete ac;
}

int main() {
  for (double ka = 1; ka < 10; ka += 0.2) {
    freq_scan(ka);
    std::cout << ka << std::endl;
  }
  return 0;
}
