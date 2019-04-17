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

Eigen::VectorXd multi_homo(double ka, int nx, int ny, int cx, int cy, int xx,
                           int yy, int ax, int ay, const Eigen::VectorXd& x0,
                           int n_iter = 1, int n_hop = 1, double d = 0.2) {
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
  Eigen::MatrixXcd w(nn, ax * ay), t(nn, ax * ay);
  for (int i = 0; i < ax; i++)
    for (int j = 0; j < ay; j++) {
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
        w(k, i * ay + j) = v[k].Bv()(0);
        t(k, i * ay + j) = v[k].Bv()(1);
      }
    }

  Mismatch f(omega, w, t, norm_mat, w_rve, h_rve, node_dens);
  std::ofstream file("iterations_" + std::to_string(ka) + ".dat");
  Eigen::VectorXd rst = BasinHopping(n_iter, n_hop, f, x0, &file);
  file.close();

  delete sol;
  delete ac;
  return rst;
}

int main(int argc, char** argv) {
  if (argc != 9) exit_error_msg({"Fiber numbers needed."});

  // nx: the number of fibers in the array in x direction.
  // ny: the number of fibers in the array in y direction.
  // cx: the number of fibers in the RVE in x direction.
  // cy: the number of fibers in the RVE in y direction.
  // ax: the number of RVEs for averaging in x direction.
  // ay: the number of RVEs for averaging in y direction.
  int ax(atoi(argv[1])), ay(atoi(argv[2])), cx(atoi(argv[3])),
      cy(atoi(argv[4])), nx(atoi(argv[5])), ny(atoi(argv[6])),
      xx(atoi(argv[7])), yy(atoi(argv[8]));

  std::string fn("");
  fn += std::to_string(ax) + "x" + std::to_string(ay) + "_of_";
  fn += std::to_string(cx) + "x" + std::to_string(cy) + "_in_";
  fn += std::to_string(nx) + "x" + std::to_string(ny) + "_start_from_";
  fn += std::to_string(xx) + "x" + std::to_string(yy) + ".txt";

  std::ofstream file(fn);

  double d = 0.2;

  int N = 18;
  Eigen::VectorXd x0(4);
  x0.setOnes();

  for (int i = 0; i < N; i++) {
    double ka = 1 + 0.05 * i;
    x0 = multi_homo(ka, nx, ny, cx, cy, xx, yy, ax, ay, x0);
    file << ka << ": " << x0.transpose() << std::endl;
  }
  file.close();

  return 0;
}
