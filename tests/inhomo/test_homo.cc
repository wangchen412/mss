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
  Material rubber{1300, 1.41908e9, 0.832e9};
  Matrix m{rubber, 1e6};

  IncidentPlaneSH in{m, pi / 3, 1e-6, 2};
  Boundary<AP, 10> b{50 * m.KT(), {{-6e-3, 6e-3}, {6e-3, -6e-3}}, &m};
};

TEST_F(HomoTest, HnG) {
  VectorXcd w(b.NumNode()), t(b.NumNode());
  for (size_t i = 0; i < b.NumNode(); i++) {
    Vector2cd tmp = in.Effect(b.Node(i)).Bv();
    w(i) = tmp(0);
    t(i) = tmp(1);
  }
  VectorXcd err = b.MatrixH() * w - b.MatrixG() * t;
  std::cout << err.norm() << std::endl;
}

}  // namespace test

}  // namespace mss
