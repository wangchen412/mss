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

// nl: number of layers around the unit cell.

Eigen::VectorXd homo(double ka, int nl, const Eigen::VectorXd& x0,
                     int n_iter = 2, int n_hop = 2, double d = 0.2) {
  Eigen::VectorXd rst(4);
  double omega = ka * 16576.24319112025;
  const Material steel(7670, 116e9, 84.3e9), lead(11400, 36e9, 8.43e9);
  const Material norm_mat({11400, 11400}, 0, {84e9, 84e9});

  Matrix m(steel, omega);
  IncidentPlaneSH in(m);

  auto fc = new FiberConfig<AP>("1", 16, 200, 0.06, lead, &m);
  InhomoPtrs<AP> fibers;

  PosiVect llp(-nl * d, -nl * d);

  for (int i = 0; i < 2 * nl + 1; i++)
    for (int j = 0; j < 2 * nl + 1; j++)
      fibers.push_back(new Fiber<AP>(fc, llp + PosiVect(i * d, j * d)));
  auto ac = new AssemblyConfig<AP>("1", fibers, &m);
  auto sol = new Solution<AP>(ac, {&in}, m);
  sol->Solve();

  int nn = 800;
  double node_dens = nn / 4.0 / d;  // 800 nodes.

  Eigen::VectorXcd w(nn), t(nn);
  Boundary<AP, 4> b{node_dens, {{-d / 2, d / 2}, {d / 2, -d / 2}}, &m};
  std::vector<StateAP> v(b.Node().size());
#ifdef NDEBUG
#pragma omp parallel for
#endif
  for (size_t k = 0; k < b.NumNode(); k++) {
    v[k] = sol->Resultant(b.Node(k));
    w(k) = v[k].Bv()(0);
    t(k) = v[k].Bv()(1);
  }
  Mismatch f(omega, w, t, norm_mat, d, d, node_dens);
  rst = BasinHopping(n_iter, n_hop, f, x0);
  delete sol;
  delete ac;
  return rst;
}

int main() {
  Eigen::VectorXd x0(4);
  x0.setOnes();

  std::ofstream file("layers.txt");
  for (int i = 0; i < 18; i++) {
    file << i << "\t";
    for (int j = 0; j < 10; j++)
      file << homo(1 + 0.05 * i, j, x0).transpose() << "\t";
    file << std::endl;
    std::cout << i << "th finished." << std::endl;
  }
  file.close();

  return 0;
}
