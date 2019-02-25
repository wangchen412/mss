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
  VectorXcd inv1 = in.EffectDv(f1.Node());
  VectorXcd inv2 = in.EffectDv(b2.Node());
  MatrixXcd B = b2.DispMatT(f1.Node());
  VectorXcd inv1_com = B * inv2;
  EXPECT_TRUE(ApproxVectRv(inv1, inv1_com, 2e-3, 0, true));

  // Pseudo inverse of the boundary integral, as the BEM-based NAH
  // No extra regularization except for the filtering of small singular
  // values.

  // MatrixXcd Bi = PseudoInverse(B);
  // VectorXcd inv2_com = Bi * inv1_com;
  // EXPECT_TRUE(ApproxVectRv(inv2, inv2_com, 2e-3, 0, true));
}
TEST_F(BEMTest, HoleSolution) {
  // Traction free BC.

  Boundary<AP, 10> b(20 * m.KT(), {{0, 0}, {0, 3e-3}}, &m, CIRCULAR_EXTERN);
  VectorXcd ws_com = b.MatrixH().lu().solve(in.EffectDv(b.Node()));

  FiberConfig<AP> fc_hole{"hole", 50, 10000, 3e-3, {1e-20, 0, 1e-30}, &m};
  Fiber<AP> hole{&fc_hole};
  hole.SetCoeff(hole.DSolve(in.EffectBv(hole.Node())));
  VectorXcd ws_ref = hole.ScatterDv(b.Node()) + in.EffectDv(b.Node());

  EXPECT_FALSE(ApproxVectRv(ws_ref, ws_com, 2e-4));
  EXPECT_TRUE(ApproxVectRv(ws_ref, ws_com, 3e-4, 0, true));
}
TEST_F(BEMTest, FixedSolution) {
  // Fix BC.
  IncidentPlaneSH in2{m, pi / 6 * 5, 1e-6};

  Boundary<AP, 10> b(50 * m.KT(), {{0, 0}, {0, 3e-3}}, &m, CIRCULAR_EXTERN);
  VectorXcd ts_com =
      -b.MatrixG().lu().solve(in.EffectDv(b.Node()) + in2.EffectDv(b.Node()));

  FiberConfig<AP> fc_rigid{"rigid", 50, 10000, 3e-3, {1e20, 0, 1e20}, &m};
  Fiber<AP> rigid{&fc_rigid};
  rigid.SetCoeff(
      rigid.DSolve(in.EffectBv(rigid.Node()) + in2.EffectBv(rigid.Node())));
  VectorXcd ts_ref(b.Node().size());
  VectorXcd bv = rigid.ScatterBv(b.Node()) + in.EffectBv(b.Node()) +
                 in2.EffectBv(b.Node());
  for (long i = 0; i < ts_ref.size(); i++) ts_ref(i) = bv(i * 2 + 1);

  EXPECT_FALSE(ApproxVectRv(ts_ref, ts_com, 0.03));
  EXPECT_TRUE(ApproxVectRv(ts_ref, ts_com, 0.04, 0, true));
}
TEST_F(BEMTest, Scattering) {
  Boundary<AP, 14> b0(50 * m.KT(), {{0, 0}, {0, 3e-3}}, &m, CIRCULAR_EXTERN);
  Boundary<AP, 14> b1(50 * m.KT(), {{0, 0}, {0, 3e-3}}, &ff, CIRCULAR);

  // Scattering problem solution.
  VectorXcd ws_com = (b0.MatrixH() + b0.MatrixG() * b1.DtN())
                         .lu()
                         .solve(in.EffectDv(b0.Node()));
  f1.SetCoeff(f1.DSolve(in.EffectBv(f1.Node())));
  VectorXcd ws_ref = f1.ScatterDv(b0.Node()) + in.EffectDv(b0.Node());
  EXPECT_TRUE(ApproxVectRv(ws_ref, ws_com, 3e-3, 0, true));

  // Traction determination.
  VectorXcd bv_ref = f1.ScatterBv(b1.Node()) + in.EffectBv(b1.Node());
  VectorXcd bv_com = b1.DispToEffect() * ws_com;
  EXPECT_TRUE(ApproxVectRv(bv_ref, bv_com, 1e-2, 0, true));

  // Inner field computation.
  Boundary<AP, 4> bb(m.KT(), {{1.5e-3, 0}, {2.5e-3, 0}}, &ff, CIRCULAR);
  CSCPtrs sp = bb.Node();

  VectorXcd inner_ref(sp.size() * 3), inner_com(sp.size() * 3);
  for (size_t i = 0; i < sp.size(); i++) {
    inner_ref.segment<3>(i * 3) = f1.Inner(sp[i]).V();
    inner_com.segment<3>(i * 3) = b1.EffectStateMatT(sp[i]) * bv_com;
  }
  EXPECT_TRUE(ApproxVectRv(inner_ref, inner_com, 3e-3, 0, true));

  // Exterior field computation.
  Boundary<AP, 4> bb2(m.KT(), {{6e-3, 6e-3}, {5e-3, 5e-3}}, &ff, CIRCULAR);
  CSCPtrs sp2 = bb2.Node();

  VectorXcd outer_ref(sp.size() * 3), outer_com(sp.size() * 3);
  for (size_t i = 0; i < sp.size(); i++) {

    StateAP tmp0 = in.Effect(sp[i]);
    Eigen::Vector3cd inci{tmp0.Displacement().x, tmp0.Stress().x,
                          tmp0.Stress().y};

    outer_ref.segment<3>(i * 3) = inci + f1.Scatter(sp[i]).V();
    outer_com.segment<3>(i * 3) = inci - b1.EffectStateMatT(sp[i]) * bv_com;
  }
  EXPECT_TRUE(ApproxVectRv(inner_ref, inner_com, 3e-3, 0, true));
}

}  // namespace test

}  // namespace mss
