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
#include "../../src/post/Output.h"
#include "../../src/post/check/Continuity.h"
#include "../test.h"

namespace mss {

namespace test {

class HomoTest : public Test {
 protected:
  HomoTest() : Test(__FILE__, "homo") {}

  input::Solution input{path("input.txt")};
  Matrix matrix{input};
  Boundary<AP, 12> b{500, {{-0.25, 0.25}, {0.25, -0.25}}, &matrix};
};

TEST_F(HomoTest, H_G_matrices) {
  MatrixXcd e(b.NumNode(), 5);
  for (size_t j = 0; j < 5; j++) {
    IncidentPlaneSH in{matrix, pi / 10 * j};
    VectorXcd w(b.NumNode()), t(b.NumNode());
    for (size_t i = 0; i < b.NumNode(); i++) {
      Vector2cd tmp = in.Effect(b.Node(i)).Bv();
      w(i) = tmp(0);
      t(i) = tmp(1);
    }
    e.col(j) = b.MatrixH() * w - b.MatrixG() * t;
  }
  EXPECT_TRUE(e.norm() < 1e-2);
}

}  // namespace test

}  // namespace mss
