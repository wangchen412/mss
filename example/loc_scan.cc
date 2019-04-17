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

// nx: the number of fibers in the array in x direction.
// ny: the number of fibers in the array in y direction.
// cx: the number of fibers in the RVE in x direction.
// cy: the number of fibers in the RVE in y direction.
// xx: the lower left fiber in RVE is the nxth fiber from left.
// yy: the lower left fiber in RVE is the nyth fiber from bottom.

Eigen::VectorXd homo(double ka, int nx, int ny, int cx, int cy, int xx,
                     int yy, int ax, int ay, const Eigen::VectorXd& x0,
                     int n_iter = 1, int n_hop = 1, double d = 0.2) {
  Eigen::VectorXd rst(4);
  double omega = ka * 16576.24319112025;
  const Material steel(7670, 116e9, 84.3e9), lead(11400, 36e9, 8.43e9);
  const Material norm_mat({11400, 11400}, 0, {84e9, 84e9});

  Matrix m(steel, omega);
  IncidentPlaneSH in(m);

  auto fc = new FiberConfig<AP>("1", 20, 200, 0.06, lead, &m);
  InhomoPtrs<AP> fibers;

  PosiVect llp(-(nx - 1) / 2.0 * d, -(ny - 1) / 2.0 * d);

  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      fibers.push_back(new Fiber<AP>(fc, llp + PosiVect(i * d, j * d)));
  auto ac = new AssemblyConfig<AP>("1", fibers, &m);
  auto sol = new Solution<AP>(ac, {&in}, m);
  sol->Solve();

  double w_rve = cx * d;  // width of the RVE
  double h_rve = cy * d;  // height of the RVE

  int nn = 800;
  double node_dens = nn / 2.0 / (w_rve + h_rve);  // 800 nodes.

  for (int i = 0; i < ax; i++)
    for (int j = 0; j < ay; j++) {
      // Eigen::MatrixXcd w(nn, ax * ay), t(nn, ax * ay);
      Eigen::VectorXcd w(nn), t(nn);
      Boundary<AP, 4> b{
          node_dens,
          {{-d / 2 * nx + d * (xx + i), -d / 2 * ny + d * (yy + j) + h_rve},
           {-d / 2 * nx + d * (xx + i) + w_rve, -d / 2 * ny + d * (yy + j)}},
          &m};
      std::vector<StateAP> v(b.Node().size());
#ifdef NDEBUG
#pragma omp parallel for
#endif
      for (size_t k = 0; k < b.NumNode(); k++) {
        v[k] = sol->Resultant(b.Node(k));
        w(k) = v[k].Bv()(0);
        t(k) = v[k].Bv()(1);
      }
      Mismatch f(omega, w, t, norm_mat, w_rve, h_rve, node_dens);
      std::ofstream file;
      file.open(std::to_string(i) + "x" + std::to_string(j) + ".dat",
                std::ofstream::out | std::ofstream::app);
      rst = BasinHopping(n_iter, n_hop, f, x0);
      file << ka << ": " << rst.transpose() << std::endl;
      file.close();
    }
  delete sol;
  delete ac;
  return rst;
}

int main() {
  int N = 18;
  Eigen::VectorXd x0(4);
  x0.setOnes();

  for (int i = 0; i < N; i++) {
    double ka = 1 + 0.05 * i;
    x0 = homo(ka, 18, 18, 1, 1, 0, 0, 9, 9, x0);
  }

  return 0;
}
