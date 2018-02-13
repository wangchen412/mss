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

class BEMTest : public Test {
 protected:
  Material rubber{1300, 1.41908e9, 0.832e9}, lead{11400, 36.32496e9, 8.43e9};
  Matrix m{rubber, 1e6};
  Matrix ff{lead, 1e6};
  IncidentPlaneSH in{m, pi / 3, 1e-6, 2};
  FiberConfig<AP> fc1{"1", 20, 1413, 3e-3, lead, &m};
  Fiber<AP> f1{&fc1};

  Boundary<AP, 4> b1{60 * m.KT(), {{0, 0}, {0, 3e-3}}, &ff, CIRCULAR};
  Boundary<AP, 10> b2{20 * m.KT(), {{-6e-3, 6e-3}, {6e-3, -6e-3}}, &m};
};

TEST_F(BEMTest, Constructor) {
  EXPECT_EQ(f1.NumNode(), b1.NumNode());

  VectorXcd xf(f1.NumNode()), yf(f1.NumNode());
  VectorXcd xb(b1.NumNode()), yb(b1.NumNode());
  for (size_t i = 0; i < f1.NumNode(); i++) {
    xf(i) = f1.Node(i)->PositionGLB().x;
    yf(i) = f1.Node(i)->PositionGLB().y;
    xb(i) = b1.Node(i)->PositionGLB().x;
    yb(i) = b1.Node(i)->PositionGLB().y;
  }

  EXPECT_TRUE(ApproxVectRv(xf, xb, 1e-12));
  EXPECT_TRUE(ApproxVectRv(yf, yb, 1e-12));
}
TEST_F(BEMTest, InfluenceMatrices) {
  f1.SetCoeff(f1.DSolve(in.EffectBv(f1.Node())));

  VectorXcd w(f1.NumNode()), t(f1.NumNode());
  for (size_t i = 0; i < f1.NumNode(); i++) {
    Vector2cd tmp = (f1.Scatter(f1.Node(i)) + in.Effect(f1.Node(i))).Bv();
    w(i) = tmp(0);
    t(i) = tmp(1);
  }

  VectorXcd com = b1.MatrixH().lu().solve(b1.MatrixG() * t);
  EXPECT_TRUE(ApproxVectRv(w, com, 2e-3, 0, true));
}
TEST_F(BEMTest, DtN) {
  f1.SetCoeff(f1.DSolve(in.EffectBv(f1.Node())));

  VectorXcd ref1 = f1.ScatterBv(f1.Node()) + in.EffectBv(f1.Node());
  VectorXcd w(f1.NumNode());
  for (size_t i = 0; i < f1.NumNode(); i++) w(i) = ref1(i * 2);

  VectorXcd com1 = b1.DispToEffect() * w;
  EXPECT_TRUE(ApproxVectRv(ref1, com1, 2e-2, 0, true));
}
TEST_F(BEMTest, DispMatT) {
  VectorXcd ref = in.EffectDv(f1.Node());
  VectorXcd com = b2.DispMatT(f1.Node()) * in.EffectDv(b2.Node());
  EXPECT_TRUE(ApproxVectRv(ref, com, 2e-3, 0, true));
}

}  // namespace test

}  // namespace mss
