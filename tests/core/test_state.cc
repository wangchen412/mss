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

#include <gtest/gtest.h>
#include "../../src/core/State.h"

namespace mss {

class StateTest : public testing::Test {
 protected:
  StateTest()
      : cs1(3, 4, t),
        cs2(-1, -2, -t, &cs1),
        cs3(-1, -10, -pi / 2, &cs2),
        a(),
        b(DispAP(1.9), StressAP(2.0 + 8.0 * ii, 3.7), &cs1),
        c(4.0 + 6.0 * ii, 5.5, 6.0 + 4.0 * ii, &cs3),
        d(DispAP(1.9), StressAP(2.0 + 8.0 * ii, 3.7).Rotate(pi / 2), &cs2),
        aa(),
        bb(DispIP(1.0 + 9.0 * ii, 2.8),
           StressIP(3.0 + 7.0 * ii, 4.6, 5.0 + 5.0 * ii), &cs1),
        cc(6.4, 7.0 + 3.0 * ii, 8.2, 9.0 + 1.0 * ii, 1.9, &cs3),
        dd(DispIP(1.0 + 9.0 * ii, 2.8).Rotate(pi / 2),
           StressIP(3.0 + 7.0 * ii, 4.6, 5.0 + 5.0 * ii).Rotate(pi / 2),
           &cs2) {}

  const double at34 = atan(0.75);
  const double at54 = atan(1.25);
  const double t = pi - at34;

  CS cs1, cs2, cs3, cs4;
  StateAP a, b, c, d;
  StateIP aa, bb, cc, dd;
};

TEST_F(StateTest, Constructors) {
  EXPECT_EQ(a.Displacement(), 0.0);
  EXPECT_EQ(a.Stress(), StressAP(0.0, 0.0));
  EXPECT_EQ(a.Basis(), nullptr);
  EXPECT_EQ(b.Displacement(), 1.9);
  EXPECT_EQ(b.Stress(), StressAP(2.0 + 8.0 * ii, 3.7));
  EXPECT_EQ(b.Basis(), &cs1);
  EXPECT_EQ(c.Displacement(), 4.0 + 6.0 * ii);
  EXPECT_EQ(c.Stress(), StressAP(5.5, 6.0 + 4.0 * ii));
  EXPECT_EQ(c.Basis(), &cs3);
  EXPECT_EQ(aa.Displacement(), DispIP(0.0, 0.0));
  EXPECT_EQ(aa.Stress(), StressIP(0.0, 0.0, 0.0));
  EXPECT_EQ(aa.Basis(), nullptr);
  EXPECT_EQ(bb.Displacement(), DispIP(1.0 + 9.0 * ii, 2.8));
  EXPECT_EQ(bb.Stress(), StressIP(3.0 + 7.0 * ii, 4.6, 5.0 + 5.0 * ii));
  EXPECT_EQ(bb.Basis(), &cs1);
  EXPECT_EQ(cc.Displacement(), DispIP(6.4, 7.0 + 3.0 * ii));
  EXPECT_EQ(cc.Stress(), StressIP(8.2, 9.0 + 1.0 * ii, 1.9));
  EXPECT_EQ(cc.Basis(), &cs3);
  EXPECT_EQ(StateAP(d), d);
  EXPECT_EQ(StateIP(dd), dd);
}
TEST_F(StateTest, in) {
  // a:
  EXPECT_EQ(a.in(&cs1), StateAP(0, 0, 0, &cs1));
  EXPECT_EQ(a.in(&cs2), StateAP(0, 0, 0, &cs2));
  EXPECT_EQ(a.in(&cs3), StateAP(0, 0, 0, &cs3));
  // b:
  EXPECT_EQ(b.in(&cs1), b);
  EXPECT_EQ(b.in(&cs2),
            StateAP(1.9, StressAP(2.0 + 8.0 * ii, 3.7).Rotate(-t), &cs2));
  EXPECT_EQ(
      b.in(&cs3),
      StateAP(1.9, StressAP(2.0 + 8.0 * ii, 3.7).Rotate(-pi / 2 - t), &cs3));
  // d:
  EXPECT_EQ(d.in(&cs1), b.in(&cs1));
  EXPECT_EQ(d.in(&cs2), b.in(&cs2));
  EXPECT_EQ(d.in(&cs3), b.in(&cs3));

  // aa:
  EXPECT_EQ(aa.in(&cs1), StateIP(0, 0, 0, 0, 0, &cs1));
  EXPECT_EQ(aa.in(&cs2), StateIP(0, 0, 0, 0, 0, &cs2));
  EXPECT_EQ(aa.in(&cs3), StateIP(0, 0, 0, 0, 0, &cs3));
  // bb:
  EXPECT_EQ(bb.in(&cs1), bb);
  EXPECT_EQ(bb.in(&cs2),
            StateIP(DispIP(1.0 + 9.0 * ii, 2.8).Rotate(-t),
                    StressIP(3.0 + 7.0 * ii, 4.6, 5.0 + 5.0 * ii).Rotate(-t),
                    &cs2));
  EXPECT_EQ(
      bb.in(&cs3),
      StateIP(
          DispIP(1.0 + 9.0 * ii, 2.8).Rotate(-pi / 2 - t),
          StressIP(3.0 + 7.0 * ii, 4.6, 5.0 + 5.0 * ii).Rotate(-pi / 2 - t),
          &cs3));
  // dd:
  EXPECT_EQ(dd.in(&cs1), bb.in(&cs1));
  EXPECT_EQ(dd.in(&cs2), bb.in(&cs2));
  EXPECT_EQ(dd.in(&cs3), bb.in(&cs3));
}

}  // namespace mss
