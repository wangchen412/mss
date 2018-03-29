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
  Matrix m2{rubber, 1e6};
  IncidentPlaneSH in1{m, pi / 3, 1, 2};
  IncidentPlaneSH in2{m2, pi / 3, 1, 2};
  FiberConfig<AP> fc1{"1", 20, 213, 3e-3, lead, &m};
  FiberConfig<AP> fc2{"1", 25, 400, 3e-3, lead, &m2};
  Fiber<AP> f1{&fc1, {10e-3, 6e-3}};
  Fiber<AP> f2{&fc1, {10e-3, 6e-3}};
  Fiber<AP> f3{&fc1, {50e-3, 50e-3}};
  Fiber<AP> f4{&fc2, {0, 0}};

  Boundary<AP> b1{100 * m.KT(), {{0, 12e-3}, {23e-3, 0}}, &m};
  // Boundary<StateAP, 2> b2{
  //     100 * m.KT(), {{10e-3, 6e-3}, {0, 0}}, &m, CIRCULAR};
  double a{6e-3};
  Boundary<AP> b3{10 * m2.KT(), {{-a, a}, {a, -a}}, &m2};
};

TEST_F(BoundaryTest, Constructor) {
  EXPECT_EQ(b1.Node().size(), 10992);
  EXPECT_EQ(b3.Edge().size(), 4);
  EXPECT_EQ(b3.Edge(0).size(), b3.NumNode() / 4);
  EXPECT_EQ(b3.Edge(1).size(), b3.NumNode() / 4);
  EXPECT_EQ(b3.Edge(2).size(), b3.NumNode() / 4);
  EXPECT_EQ(b3.Edge(3).size(), b3.NumNode() / 4);

  b1.ReverseEdge();
  for (size_t i = 0; i < b1.Edge(0).size(); i++) {
    EXPECT_TRUE(ApproxRv(b1.Edge(0)[i]->Position().x + 23e-3,
                         b1.Edge(2)[i]->Position().x, 1e-12));
    EXPECT_TRUE(ApproxRv(b1.Edge(0)[i]->Position().y,
                         b1.Edge(2)[i]->Position().y, 1e-12));
  }
  for (size_t i = 0; i < b1.Edge(1).size(); i++) {
    EXPECT_TRUE(ApproxRv(b1.Edge(1)[i]->Position().x,
                         b1.Edge(3)[i]->Position().x, 1e-12));
    EXPECT_TRUE(ApproxRv(b1.Edge(1)[i]->Position().y + 12e-3,
                         b1.Edge(3)[i]->Position().y, 1e-12));
  }

  b3.ReverseEdge();
  for (size_t i = 0; i < b3.Edge(0).size(); i++) {
    EXPECT_TRUE(ApproxRv(b3.Edge(0)[i]->Position().x + 2 * a,
                         b3.Edge(2)[i]->Position().x, 1e-12));
    EXPECT_TRUE(ApproxRv(b3.Edge(0)[i]->Position().y,
                         b3.Edge(2)[i]->Position().y, 1e-12));
  }
  for (size_t i = 0; i < b3.Edge(1).size(); i++) {
    EXPECT_TRUE(ApproxRv(b3.Edge(1)[i]->Position().x,
                         b3.Edge(3)[i]->Position().x, 1e-12));
    EXPECT_TRUE(ApproxRv(b3.Edge(1)[i]->Position().y + 2 * a,
                         b3.Edge(3)[i]->Position().y, 1e-12));
  }
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

TEST_F(BoundaryTest, DISABLED_ColloMat_Rectangular) {
  MatrixXcd c = b1.ColloMatT();
  // MatrixXcd d = PseudoInverse(c);
  MatrixXcd m = b1.ModeMatT(f1.Node());

  VectorXcd v = c.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV)
                    .solve(in1.EffectBv(b1.Node()));

  VectorXcd bv_in = in1.EffectBv(f1.Node());
  VectorXcd bv_bd = m * v;
  // VectorXcd bv_bd = m * d * in1.EffectBv(b2.Node());
  EXPECT_TRUE(ApproxVectRv(bv_in, bv_bd, 1e-4, 0, true));
}
// TEST_F(BoundaryTest, ColloMat_Circular) {
//   MatrixXcd c = b2.ColloMatT();
//   MatrixXcd d = PseudoInverse(c);
//   MatrixXcd m = b2.ModeMatT(f1.Node());

//   VectorXcd bv_in = in1.EffectBv(f1.Node());
//   VectorXcd bv_bd = m * d * in1.EffectBv(b2.Node());
//   EXPECT_TRUE(ApproxVectRv(bv_in, bv_bd, 1e-4, 0, true));
// }

class AssemBoundaryTest : public Test {
 protected:
  AssemBoundaryTest() : Test(__FILE__, "boundary") {}

  input::Solution s{path("input.txt")};
  Matrix matrix{s};
  AssemblyConfig<AP> c1{s.config(), &matrix};
  AssemblyConfig<AP> c2{s.config(), &matrix};
  AssemblyConfig<AP> c3{s.assembly_config()[1], &matrix};
  IncidentPlaneSH inSH1{matrix, s.incident()[0]};
};
TEST_F(AssemBoundaryTest, CSolve) {
  c1.CSolve({&inSH1});
  c2.CSolve(c2.BdIntMatT() * inSH1.EffectBv(c2.Node()));
  for (int i = 0; i < 3; i++)
    EXPECT_TRUE(ApproxVectRv(c1.inhomo(i)->ScatterCoeff(),
                             c2.inhomo(i)->ScatterCoeff(), 1e-4, 15));
}
TEST_F(AssemBoundaryTest, DSolve) {
  c1.DSolve({&inSH1});
  c2.DSolve(c2.BdIntMatT() * inSH1.EffectBv(c2.Node()));
  for (int i = 0; i < 3; i++)
    EXPECT_TRUE(ApproxVectRv(c1.inhomo(i)->ScatterCoeff(),
                             c2.inhomo(i)->ScatterCoeff(), 1e-3, 15));
}
TEST_F(AssemBoundaryTest, CSolve_InvMat) {
  c1.CSolve({&inSH1});
  VectorXcd solution = c2.GramMat().inverse() * c2.ColloMat().transpose() *
                       c2.IncVec({&inSH1});
  for (int i = 0; i < 3; i++) {
    VectorXcd rr = c1.inhomo(i)->ScatterCoeff();
    VectorXcd cc = solution.segment(61 * i, 61);
    EXPECT_TRUE(ApproxVectRv(rr, cc, 1e-3, 15));
  }
}
TEST_F(AssemBoundaryTest, DSolve_InvMat) {
  c1.DSolve({&inSH1});
  VectorXcd solution = c2.DcMat().inverse() * c2.Trans_IncVec({&inSH1});
  for (int i = 0; i < 3; i++) {
    VectorXcd rr = c1.inhomo(i)->ScatterCoeff();
    VectorXcd cc = solution.segment(61 * i, 61);
    EXPECT_TRUE(ApproxVectRv(rr, cc));
  }
}
TEST_F(AssemBoundaryTest, CSolve_InvMatBi) {
  c1.CSolve({&inSH1});
  VectorXcd solution = c2.GramMat().inverse() * c2.ColloMat().transpose() *
                       c2.BdIntMatT() * inSH1.EffectBv(c2.Node());
  for (int i = 0; i < 3; i++) {
    VectorXcd rr = c1.inhomo(i)->ScatterCoeff();
    VectorXcd cc = solution.segment(61 * i, 61);
    EXPECT_TRUE(ApproxVectRv(rr, cc, 1e-3, 15));
  }
}
TEST_F(AssemBoundaryTest, DSolve_InvMatBi) {
  c1.DSolve({&inSH1});
  VectorXcd solution =
      c2.DcMat().inverse() *
      c2.Trans_IncVec(c2.BdIntMatT() * inSH1.EffectBv(c2.Node()));
  for (int i = 0; i < 3; i++) {
    VectorXcd rr = c1.inhomo(i)->ScatterCoeff();
    VectorXcd cc = solution.segment(61 * i, 61);
    EXPECT_TRUE(ApproxVectRv(rr, cc, 1e-3, 15));
  }
}
TEST_F(AssemBoundaryTest, PlaneEDMat_Single) {
  FiberConfig<AP> fc(s.fiber_config()[2], &matrix);
  Fiber<AP> f(&fc, {80e-3, 10e-3});  // This fiber is acting as a scatterer.
  f.SetCoeff(f.CSolve({&inSH1}));

  VectorXcd ref = inSH1.EffectDv(c3.Node()) + f.ScatterDv(c3.Node());
  VectorXcd in = inSH1.EffectDv(c3.Node_in()) + f.ScatterDv(c3.Node_in());

  MatrixXcd ext_m = c3.Boundary().PlaneEDMat(c3.Node_in());
  VectorXcd com = ext_m * in;

  EXPECT_TRUE(ApproxVectRv(ref, com, 6e-3));
}
TEST_F(AssemBoundaryTest, PlaneEDMat_Multiple) {
  FiberConfig<AP> fc(s.fiber_config()[2], &matrix);
  Fiber<AP> f(&fc, {140e-3, 50e-3});  // This fiber is acting as a scatterer.
  f.SetCoeff(f.CSolve({&inSH1}));

  AssemblyConfig<AP> c4{s.assembly_config()[2], &matrix};
  VectorXcd ref = inSH1.EffectDv(c4.Node()) + f.ScatterDv(c4.Node());
  VectorXcd in = inSH1.EffectDv(c4.Node_in()) + f.ScatterDv(c4.Node_in());

  MatrixXcd ext_m = c4.Boundary().PlaneEDMat(c4.Node_in());
  VectorXcd com = ext_m * in;

  EXPECT_TRUE(ApproxVectRv(ref, com, 6e-3));
}
TEST_F(AssemBoundaryTest, CylinEDMat) {
  FiberConfig<AP> fc(s.fiber_config()[2], &matrix);
  Fiber<AP> f(&fc, {140e-3, 50e-3});  // This fiber is acting as a scatterer.
  f.SetCoeff(f.CSolve({&inSH1}));

  AssemblyConfig<AP> c4{s.assembly_config()[2], &matrix};
  c4.DSolve(inSH1.EffectBv(c4.Node_in()) + f.ScatterBv(c4.Node_in()));

  for (size_t i = 0; i < 4; i++) {
    VectorXcd ref = inSH1.EffectDv(c4.Edge(i)) + f.ScatterDv(c4.Edge(i));
    VectorXcd com = c4.CylinEDMat(c4.Edge(i)) * c4.ScatterCoeff();
    EXPECT_TRUE(ApproxVectRv(ref, com, 5e-4));
  }
}
TEST_F(AssemBoundaryTest, CylinEBMat) {
  FiberConfig<AP> fc(s.fiber_config()[2], &matrix);
  Fiber<AP> f(&fc, {140e-3, 50e-3});  // This fiber is acting as a scatterer.
  f.SetCoeff(f.CSolve({&inSH1}));

  AssemblyConfig<AP> c4{s.assembly_config()[2], &matrix};
  c4.DSolve(inSH1.EffectBv(c4.Node_in()) + f.ScatterBv(c4.Node_in()));

  for (size_t i = 0; i < 4; i++) {
    VectorXcd ref = inSH1.EffectBv(c4.Edge(i)) + f.ScatterBv(c4.Edge(i));
    VectorXcd com = c4.CylinEBMat(c4.Edge(i)) * c4.ScatterCoeff();
    EXPECT_TRUE(ApproxVectRv(ref, com, 1e-3));
  }
}

}  // namespace test

}  // namespace mss
