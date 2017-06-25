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

// Test Coordinate System Class

#include <gtest/gtest.h>
#include "../../src/core/CS.h"

namespace mss {

namespace test {

class CSTest : public testing::Test {
 protected:
  CSTest() : a(), b(3, 4, t), c(-1, -2, -t, &b), d(4, -5, -pi / 2, &a) {}

  const double at34 = atan(0.75);
  const double at54 = atan(1.25);
  const double t = pi - at34;

  CS a, b, c, d;
};

TEST_F(CSTest, Constructors) {
  EXPECT_EQ(a.Position(), PosiVect(0, 0));
  EXPECT_EQ(a.Angle(), 0);
  EXPECT_EQ(a.Basis(), nullptr);
  EXPECT_EQ(b.Position(), PosiVect(3, 4));
  EXPECT_EQ(b.Angle(), t);
  EXPECT_EQ(b.Basis(), nullptr);
  EXPECT_EQ(c.Position(), PosiVect(-1, -2));
  EXPECT_EQ(c.Angle(), -t);
  EXPECT_EQ(c.Basis(), &b);
  EXPECT_EQ(d.Position(), PosiVect(4, -5));
  EXPECT_EQ(d.Angle(), -pi / 2);
  EXPECT_EQ(d.Basis(), &a);
  EXPECT_EQ(CS(c), c);
}
TEST_F(CSTest, in) {
  // a in b:
  EXPECT_EQ(a.in(&b).Position(), PosiVect(0, 5));
  EXPECT_PRED2(angEqu, a.in(&b).Angle(), -t);
  // a in c:
  EXPECT_EQ(a.in(&c).Position(), PosiVect(-5, -5));
  EXPECT_PRED2(angEqu, a.in(&c).Angle(), 0);
  // a in d:
  EXPECT_EQ(a.in(&d).Position(), PosiVect(-5, -4));
  EXPECT_PRED2(angEqu, a.in(&d).Angle(), pi / 2);
  // a in GLB:
  EXPECT_EQ(a.inGLB().Position(), PosiVect(0, 0));
  EXPECT_PRED2(angEqu, a.inGLB().Angle(), 0);

  // b in a:
  EXPECT_EQ(b.in(&a).Position(), PosiVect(3, 4));
  EXPECT_PRED2(angEqu, b.in(&a).Angle(), t);
  // b in c:
  EXPECT_EQ(b.in(&c).Position(), PosiVect(-2, -1));
  EXPECT_PRED2(angEqu, b.in(&c).Angle(), t);
  // b in d:
  EXPECT_EQ(b.in(&d).Position(), PosiVect(-9, -1));
  EXPECT_PRED2(angEqu, b.in(&d).Angle(), t + pi / 2);
  // b in GLB:
  EXPECT_EQ(b.inGLB().Position(), PosiVect(3, 4));
  EXPECT_PRED2(angEqu, b.inGLB().Angle(), t);

  // c in a:
  EXPECT_EQ(c.in(&a).Position(), PosiVect(5, 5));
  EXPECT_PRED2(angEqu, c.in(&a).Angle(), 0);
  // c in b:
  EXPECT_EQ(c.in(&b).Position(), PosiVect(-1, -2));
  EXPECT_PRED2(angEqu, c.in(&b).Angle(), -t);
  // c in d:
  EXPECT_EQ(c.in(&d).Position(), PosiVect(-10, 1));
  EXPECT_PRED2(angEqu, c.in(&d).Angle(), pi / 2);
  // c in GLB:
  EXPECT_EQ(c.in(nullptr).Position(), PosiVect(5, 5));
  EXPECT_PRED2(angEqu, c.in(nullptr).Angle(), 0);

  // d in a:
  EXPECT_EQ(d.in(&a).Position(), PosiVect(4, -5));
  EXPECT_PRED2(angEqu, d.in(&a).Angle(), -pi / 2);
  // d in b:
  double tmp = std::sqrt(66 - 10 * std::sqrt(41) * cos(at54 - at34)) / 5;
  EXPECT_EQ(d.in(&b).Position(), PosiVect(-5 - tmp * 3, 5 + tmp * 4));
  EXPECT_PRED2(angEqu, d.in(&b).Angle(), pi / 2 + at34);
  // d in c:
  EXPECT_EQ(d.in(&c).Position(), PosiVect(-1, -10));
  EXPECT_PRED2(angEqu, d.in(&c).Angle(), -pi / 2);
  // d in GLB:
  EXPECT_EQ(d.inGLB().Position(), PosiVect(4, -5));
  EXPECT_PRED2(angEqu, d.inGLB().Angle(), -pi / 2);
}

}  // namespace test

}  // namespace mss
