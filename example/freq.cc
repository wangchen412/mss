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

// nc: number of layers of center inhomos.
// ns: number of layers of surrounding inhomos.
// Fox example, nc = 1, ns = 2 will be 5 x 5 with considering the center one
// as the RVE; nc = 2, ns = 1 will be 4 x 4 with considering the center 2 x 2
// as the RVE.
// np: number of nodes along RVE boundary.
Eigen::MatrixXd homo(double ka, int nc, int ns, size_t np) {
  double omega = ka * 16576.24319112025;

  const Material steel(7670, 116e9, 84.3e9), lead(11400, 36e9, 8.43e9);
  Matrix m(steel, omega);
  IncidentPlaneSH in(m);
  auto fc = new FiberConfig<AP>("1", 20, 200, 0.06, lead, &m);
  InhomoPtrs<AP> fibers;

  int nn = nc + ns * 2;
  double d = 0.2;
  PosiVect lbp(-(nn - 1) / 2.0 * d, -(nn - 1) / 2.0 * d);

  for (int i = 0; i < nn; i++)
    for (int j = 0; j < nn; j++)
      fibers.push_back(new Fiber<AP>(fc, lbp + PosiVect(i * d, j * d)));
  auto ac = new AssemblyConfig<AP>("1", fibers, &m);
  auto sol = new Solution<AP>(ac, {&in}, m);
  sol->Solve();

  double hw = nc / 2.0 * d;  // half width of the RVE
  double node_dens = np / hw / 8;
  Boundary<AP, 4> b{node_dens, {{-hw, hw}, {hw, -hw}}, &m};
  Eigen::VectorXcd w(b.NumNode()), t(b.NumNode());

  std::cout << "Number of nodes: " << b.NumNode() << ":\t";

  std::vector<StateAP> v(b.Node().size());
#ifdef NDEBUG
#pragma omp parallel for
#endif
  for (size_t i = 0; i < b.NumNode(); i++) {
    v[i] = sol->Resultant(b.Node(i));
    w(i) = v[i].Bv()(0);
    t(i) = v[i].Bv()(1);
  }

  Mismatch f(omega, w, t, {{11400, 11400}, 0, {84e9, 84e9}}, hw * 2, hw * 2,
             node_dens);
  std::ofstream file("iterations_" + std::to_string(ka) + ".dat");
  // Eigen::VectorXd rst = NelderMead(f, Eigen::Vector4d::Ones(), &file);
  Eigen::VectorXd rst =
      BasinHopping(10, 40, f, Eigen::Vector4d::Ones(), &file);
  file.close();

  delete sol;
  delete ac;
  return rst.transpose();
}

int main() {
  for (int i = 0; i < 9; i++)
    std::cout << homo(1.2, 2, 5, 400 + i * 100) << std::endl;
  return 0;
}
