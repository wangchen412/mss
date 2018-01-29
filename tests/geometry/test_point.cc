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

#include "../../src/post/geometry/Point.h"
#include "../test.h"

namespace mss {

namespace test {

class PointTest : public Test {
 protected:
  PointTest() : Test(__FILE__) {}
  SolutionAP s{path("input.txt")};
};

TEST_F(PointTest, Ctor_Computation) {
  s.Solve();
  post::PointAP p1(&s, {}), p2(&s, {1, 2}, pi / 3, "p2");

  EXPECT_EQ(p1.LocalCS()->PositionGLB(), PosiVect(0, 0));
  EXPECT_EQ(p1.LocalCS()->AngleGLB(), 0);
  EXPECT_EQ(p1.LocalCS()->Basis(), nullptr);

  EXPECT_EQ(p2.LocalCS()->PositionGLB(), PosiVect(1, 2));
  EXPECT_EQ(p2.LocalCS()->AngleGLB(), pi / 3);
  EXPECT_EQ(p2.LocalCS()->Basis(), nullptr);

  post::PointAP p3(&s, {1, 1.03}, 1.5707963267948966);
  post::PointAP p4(&s, {0, 0});
  post::PointAP p5(&s, {18e-3, 18e-3});
  post::PointAP p6(&s, {-18e-3, -18e-3});
  post::PointAP p7(&s, {-18e-3, 18e-3});
  post::PointAP p8(&s, {18e-3, -18e-3});

  dcomp w = {1.6587125982302423e-06, -1.7477625866189296e-06};
  dcomp tzx = {1775443.2580089821, 1989770.6209655141};
  dcomp tzy = {2696226.105811324, -846321.05973796837};
  StateAP st = {w, tzx, tzy, p3.LocalCS()};

  EXPECT_TRUE(st.isApprox(p3.State(), 1e-6));
}

}  // namespace test

}  // namespace mss
