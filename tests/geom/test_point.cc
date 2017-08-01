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

#include "../test.h"

namespace mss {

namespace test {

class PointTest : public Test {
protected:
  PointTest() : Test(__FILE__) {}

  SolutionAP s{path("input.txt")};
  post::PointAP p1{{1.2, 4.2}, &s}, p2{{5.6, 3.9}, &s, pi/3, "120"};
};

TEST_F(PointTest, Constructor) {
  EXPECT_EQ(p1.LocalCS()->PositionGLB(), PosiVect(1.2, 4.2));
  EXPECT_EQ(p1.LocalCS()->AngleGLB(), 0);
  EXPECT_EQ(p1.LocalCS()->Basis(), nullptr);

  EXPECT_EQ(p2.LocalCS()->PositionGLB(), PosiVect(5.6, 3.9));
  EXPECT_EQ(p2.LocalCS()->AngleGLB(), pi/3);
  EXPECT_EQ(p2.LocalCS()->Basis(), nullptr);
}

}

}  // namespace mss
