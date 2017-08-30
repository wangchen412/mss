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

class BoundaryTest : public Test {
 protected:
  Material rubber{1300, 1.41908e9, 0.832e9}, lead{11400, 36.32496e9, 8.43e9};
  Matrix m{rubber, 1.25664e6};
  IncidentPlaneSH in1{m, pi / 3, 1, 2};
  FiberConfig<StateAP> fc1{"1", 20, 213, 3e-3, lead, &m};
  Fiber<StateAP> f1{&fc1, {10e-3, 6e-3}};
  Fiber<StateAP> f2{&fc1, {10e-3, 6e-3}};
  Fiber<StateAP> f3{&fc1, {50e-3, 50e-3}};

  Boundary<StateAP, 2> b1{100 * m.KT(), {{0, 12e-3}, {23e-3, 0}}, &m};
};

TEST_F(BoundaryTest, Constructor) {
  EXPECT_EQ(b1.Node().size(), 10992);
}

TEST_F(BoundaryTest, EffectMat) {
  // The points with normal vector perpendicular to the incident plane wave
  // should have zero traction, which will cause large relative error.
  VectorXcd bv_in = in1.EffectBv(f1.Node());
  VectorXcd bv_bd = b1.EffectMatT(f1.Node()) * in1.EffectBv(b1.Node());
  EXPECT_TRUE(ApproxVectRv(bv_in, bv_bd, 1e-4));
}

TEST_F(BoundaryTest, Solve) {
  EXPECT_TRUE(ApproxVectRv(
      f1.CSolve({&in1}),
      f2.CSolve(b1.EffectBvT(&f2, in1.EffectBv(b1.Node()))), 1e-5));
  EXPECT_TRUE(ApproxVectRv(
      f1.DSolve({&in1}),
      f2.DSolve(b1.EffectBvT(&f2, in1.EffectBv(b1.Node()))), 1e-4));
}

}  // namespace test

}  // namespace mss
