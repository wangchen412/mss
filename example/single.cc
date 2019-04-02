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

Eigen::MatrixXd single_homo(double ka) {
  double omega = ka * 16576.24319112025;

  const Material steel(7670, 116e9, 84.3e9), lead(11400, 36e9, 8.43e9);
  Matrix m(steel, omega);
  IncidentPlaneSH in(m);
  FiberConfig<AP> fc{"1", 20, 1000, 0.06, lead, &m};
  Fiber<AP> f1{&fc};
  f1.SetCoeff(f1.DSolve({&in}));

  Boundary<AP, 4> b{500, {{-0.1, 0.1}, {0.1, -0.1}}, &m};
  Eigen::VectorXcd w(b.NumNode()), t(b.NumNode());

  std::vector<StateAP> v(b.Node().size());
#ifdef NDEBUG
#pragma omp parallel for
#endif
  for (size_t i = 0; i < b.NumNode(); i++) {
    v[i] = f1.Scatter(b.Node(i)) + in.Effect(b.Node(i));
    w(i) = v[i].Bv()(0);
    t(i) = v[i].Bv()(1);
  }

  Mismatch f(omega, w, t, {{11400, 11400}, 0, {84e9, 84e9}}, 0.2, 0.2);
  std::ofstream file("iterations_" + std::to_string(ka) + ".dat");
  Eigen::VectorXd rst = NelderMead(f, Eigen::Vector4d::Ones(), &file);
  file.close();
  return rst.transpose();
}

int main() {
  std::ofstream file("frequency.txt");
  int N = 20;
  double d = 1.0 / N;
  for (int i = 0; i < N; i++) {
    std::cout << i << std::endl;
    file << single_homo(1 + d * i) << std::endl;
  }
  file.close();
  return 0;
}
