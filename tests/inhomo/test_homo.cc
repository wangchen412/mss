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

#include "../../src/inhomo/Boundary.h"
#include "../test.h"

namespace mss {

namespace test {

class HomoTest : public Test {
 protected:
  HomoTest() : Test(__FILE__, "homo") {}

  input::Solution input{path("input.txt")};
  Matrix matrix{input};
  AssemblyConfig<AP> ac{input.config(), &matrix};
  Boundary<AP, 10> b{200, {{-0.3, 0.3}, {0.3, -0.3}}, &matrix};
};

TEST_F(HomoTest, HnG) {
  Material rubber{1300, 1.41908e9, 0.832e9};
  Matrix m{rubber, 1e5};
  Boundary<AP, 10> b{20 * m.KT(), {{-6e-2, 6e-2}, {6e-2, -6e-2}}, &m};
  MatrixXcd e(b.NumNode(), 5);
  for (size_t j = 0; j < 5; j++) {
    IncidentPlaneSH in{m, pi / 10 * j, 1e-6};
    VectorXcd w(b.NumNode()), t(b.NumNode());
    for (size_t i = 0; i < b.NumNode(); i++) {
      Vector2cd tmp = in.Effect(b.Node(i)).Bv();
      w(i) = tmp(0);
      t(i) = tmp(1);
    }
    e.col(j) = b.MatrixH() * w - b.MatrixG() * t;
  }

  std::cout << e.norm() << std::endl;
}
TEST_F(HomoTest, HnG2) {
  MatrixXcd e(b.NumNode(), 5);
  VectorXcd w(b.NumNode()), t(b.NumNode());
  for (size_t j = 1; j < 2; j++) {
    IncidentPlaneSH in{matrix, pi / 10 * j, 1e-6};
    // VectorXcd w(b.NumNode()), t(b.NumNode());
    for (size_t i = 0; i < b.NumNode(); i++) {
      Vector2cd tmp = in.Effect(b.Node(i)).Bv();
      w(i) = tmp(0);
      t(i) = tmp(1);
    }
    e.col(j) = b.MatrixH() * w - b.MatrixG() * t;
  }
  std::cout << e.norm() << std::endl;
  for (size_t i = 0; i < b.NumNode(); i++)
    std::cout << w(i) << "\t" << t(i) << std::endl;
}
TEST_F(HomoTest, AssemblySolve) {
  IncidentPlaneSH in{matrix, pi / 10, 1e-6};
  ac.DSolve({&in});

  VectorXcd w(b.NumNode()), t(b.NumNode());
  for (size_t i = 0; i < b.NumNode(); i++) {
    Vector2cd tmp = ac.Resultant(b.Node(i), {&in}).Bv();
    w(i) = tmp(0);
    t(i) = tmp(1);
  }

  for (size_t i = 0; i < b.NumNode(); i++)
    std::cout << w(i) << "\t" << t(i) << std::endl;

  // Eigen::MatrixXd err(30, 30);
  // for (size_t i = 0; i < 30; i++) {
  //   for (size_t j = 0; j < 30; j++) {
  //     double rho = 7600.0 + i * 126.67, mu = 84.3e9 - j * 25.29e8;
  //     Matrix mm({rho, 80e9, mu}, 55508.0206);
  //     Boundary<AP, 10> b{200, {{-0.3, 0.3}, {0.3, -0.3}}, &mm};
  //     err(i, j) = (b.MatrixH() * w - b.MatrixG() * t).norm();
  //   }
  // }

  // std::cout << err << std::endl;
}

}  // namespace test

}  // namespace mss
