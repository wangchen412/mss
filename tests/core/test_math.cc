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

// Test basic math functions and classes.

#include <gtest/gtest.h>
#include <limits.h>
#include <iomanip>
#include <iostream>
#include "../../src/core/Tensor.h"

namespace mss {

class VectorTest : public testing::Test {
 protected:
  VectorTest()
      : a(),
        b(1.2),
        c(0, 3.4),
        d(1.2, 3.4),
        aa(),
        bb(1.2 + 3.4 * ii),
        cc(0, 5.6 - 7.8 * ii),
        dd(1.2 + 3.4 * ii, 5.6 - 7.8 * ii) {}

  Vector<double> a, b, c, d;
  Vector<dcomp> aa, bb, cc, dd;
};
class TensorTest : public testing::Test {
 protected:
  TensorTest()
      : a(),
        b(1.2 + 3.4 * ii),
        c(1.2 + 3.4 * ii, 5.6 - 7.8 * ii),
        d(1.2 + 3.4 * ii, 5.6 - 7.8 * ii, 9.0 + 1.2 * ii),
        e(0, 5.6 - 7.8 * ii, 9.0 + 1.2 * ii),
        f(0, 0, 9.0 + 1.2 * ii) {}

  Tensor<dcomp> a, b, c, d, e, f;
};

TEST(MathTest, ConstantsTest) {
  EXPECT_LT(std::abs(sin(0)), epsilon);
  EXPECT_LT(std::abs(sin(pi)), epsilon);
  EXPECT_LT(std::abs(sin(pi / 6) - 0.5), epsilon);
  EXPECT_LT(std::abs(sin(pi / 3) - std::sqrt(3) / 2), epsilon);
  EXPECT_EQ(ii * ii, -1.0);
  EXPECT_EQ(1.0 / ii, -ii);
}
TEST(MathTest, angEquTest) {
  EXPECT_TRUE(angEqu(0, 0));
  EXPECT_TRUE(angEqu(pi, pi));
  EXPECT_TRUE(angEqu(2 * pi, 2 * pi));
  EXPECT_TRUE(angEqu(0, 2 * pi));
  EXPECT_TRUE(angEqu(pi, -pi));
  EXPECT_TRUE(angEqu(ee * pi, (ee + 12) * pi));
  EXPECT_FALSE(angEqu(ee * pi, (ee - 9) * pi));
}
TEST(MathTest, LpTest) {
  EXPECT_EQ(Lp<1>({3, 4}), 7);
  EXPECT_EQ(Lp<2>({3, 4}), 5);
  EXPECT_EQ(Lp<2>({1, 1}), std::sqrt(2));
  EXPECT_EQ(Lp<1>({1, 2, 3}), 6);
  EXPECT_EQ(Lp<2>({12, 15, 16}), 25);
}

TEST_F(VectorTest, ConstructorsTest) {
  EXPECT_EQ(a.x, 0.0);
  EXPECT_EQ(a.y, 0.0);
  EXPECT_EQ(b.x, 1.2);
  EXPECT_EQ(b.y, 0.0);
  EXPECT_EQ(c.x, 0.0);
  EXPECT_EQ(c.y, 3.4);
  EXPECT_EQ(d.x, 1.2);
  EXPECT_EQ(d.y, 3.4);

  EXPECT_EQ(aa.x, 0.0 + 0.0 * ii);
  EXPECT_EQ(aa.y, 0.0 + 0.0 * ii);
  EXPECT_EQ(bb.x, 1.2 + 3.4 * ii);
  EXPECT_EQ(bb.y, 0.0 + 0.0 * ii);
  EXPECT_EQ(cc.x, 0.0 + 0.0 * ii);
  EXPECT_EQ(cc.y, 5.6 - 7.8 * ii);
  EXPECT_EQ(dd.x, 1.2 + 3.4 * ii);
  EXPECT_EQ(dd.y, 5.6 - 7.8 * ii);

  // Copy Constructor:
  Vector<double> m(d);
  Vector<dcomp> mm(dd);
  EXPECT_EQ(d.x, m.x);
  EXPECT_EQ(d.y, m.y);
  EXPECT_EQ(dd.x, mm.x);
  EXPECT_EQ(dd.y, mm.y);
}
TEST_F(VectorTest, ComparisonTest) {
  EXPECT_TRUE(a == a);
  EXPECT_TRUE(d == d);
  EXPECT_TRUE(aa == aa);
  EXPECT_TRUE(dd == dd);
  EXPECT_FALSE(c == b);
  EXPECT_FALSE(cc == bb);

  Vector<double> m(d);
  Vector<dcomp> mm(dd);

  m.x += epsilon / 10;
  mm.x += epsilon / 10 * ii;
  EXPECT_TRUE(m == d);
  EXPECT_TRUE(mm == dd);

  m.y += 1.5 * epsilon;
  mm.y += 1.5 * epsilon;
  EXPECT_FALSE(m == d);
  EXPECT_FALSE(mm == dd);
}
TEST_F(VectorTest, LengthTest) {
  EXPECT_EQ(a.Length(), 0);
  EXPECT_EQ(b.Length(), 1.2);
  EXPECT_EQ(c.Length(), 3.4);
  EXPECT_EQ(d.Length(), std::sqrt(1.2 * 1.2 + 3.4 * 3.4));
}
TEST_F(VectorTest, AngleTest) {
  EXPECT_EQ(a.Angle(), 0);
  EXPECT_EQ(b.Angle(), 0);
  EXPECT_EQ(c.Angle(), pi / 2);
  EXPECT_EQ(d.Angle(), atan(d.y / d.x));

  Vector<double> nb(-b.x, -b.y);
  Vector<double> nc(-c.x, -c.y);
  EXPECT_TRUE(angEqu(nb.Angle(), pi));
  EXPECT_TRUE(angEqu(nc.Angle(), -pi / 2));

  Vector<double> d2(-d.x, d.y);
  Vector<double> d3(-d.x, -d.y);
  Vector<double> d4(d.x, -d.y);
  EXPECT_TRUE(angEqu(d2.Angle(), pi - d.Angle()));
  EXPECT_TRUE(angEqu(d3.Angle(), pi + d.Angle()));
  EXPECT_TRUE(angEqu(d4.Angle(), -d.Angle()));
}
TEST_F(VectorTest, RotateTest) {
  EXPECT_EQ(a.Rotate(0), a);
  EXPECT_EQ(a.Rotate(pi / 7), a);
  EXPECT_EQ(b.Rotate(2 * pi), b);
  EXPECT_EQ(b.Rotate(pi / 3).Angle(), -pi / 3);
  EXPECT_EQ(b.Rotate(5 * pi / 3).Angle(), pi / 3);
  EXPECT_EQ(c.Rotate(3 * pi / 7).Rotate(-3 * pi / 7), c);
  EXPECT_EQ(d.Rotate(pi / 9), d.Rotate(19 * pi / 9));

  EXPECT_EQ(aa.Rotate(0), aa);
  EXPECT_EQ(aa.Rotate(pi / 7), aa);
  EXPECT_EQ(bb.Rotate(2 * pi), bb);
  EXPECT_EQ(cc.Rotate(3 * pi / 7).Rotate(-3 * pi / 7), cc);
  EXPECT_EQ(dd.Rotate(pi / 9), dd.Rotate(19 * pi / 9));
}

TEST_F(TensorTest, ConstructorsTest) {
  EXPECT_EQ(a.xx, 0.0);
  EXPECT_EQ(a.yy, 0.0);
  EXPECT_EQ(a.xy, 0.0);
  EXPECT_EQ(b.xx, 1.2 + 3.4 * ii);
  EXPECT_EQ(b.yy, 0.0);
  EXPECT_EQ(b.xy, 0.0);
  EXPECT_EQ(c.xx, 1.2 + 3.4 * ii);
  EXPECT_EQ(c.yy, 5.6 - 7.8 * ii);
  EXPECT_EQ(c.xy, 0.0);
  EXPECT_EQ(d.xx, 1.2 + 3.4 * ii);
  EXPECT_EQ(d.yy, 5.6 - 7.8 * ii);
  EXPECT_EQ(d.xy, 9.0 + 1.2 * ii);

  // Copy Constructor:
  Tensor<dcomp> m(d);
  EXPECT_EQ(d.xx, m.xx);
  EXPECT_EQ(d.yy, m.yy);
  EXPECT_EQ(d.xy, m.xy);
}
TEST_F(TensorTest, ComparisonTest) {
  EXPECT_TRUE(a == a);
  EXPECT_TRUE(d == d);
  EXPECT_FALSE(c == d);

  Tensor<dcomp> m(d);

  m.xx += epsilon / 10 * ii;
  m.yy += epsilon / 10;
  EXPECT_TRUE(m == d);

  m.xy += 1.5 * epsilon;
  EXPECT_FALSE(m == d);
}
TEST_F(TensorTest, RotateTest) {
  EXPECT_EQ(a.Rotate(0), a);
  EXPECT_EQ(a.Rotate(pi / 7), a);
  EXPECT_EQ(b.Rotate(2 * pi), b);
  EXPECT_EQ(c.Rotate(3 * pi / 7).Rotate(-3 * pi / 7), c);
  EXPECT_EQ(d.Rotate(pi / 9), d.Rotate(19 * pi / 9));
}

}  // namespace mss
