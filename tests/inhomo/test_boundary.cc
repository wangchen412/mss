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
  FiberConfig<AP> fc1{"1", 20, 213, 3e-3, lead, &m};
  FiberConfig<AP> fc2{"1", 30, 800, 3e-3, lead, &m};
  Fiber<AP> f1{&fc1, {10e-3, 6e-3}};
  Fiber<AP> f2{&fc1, {10e-3, 6e-3}};
  Fiber<AP> f3{&fc1, {50e-3, 50e-3}};
  Fiber<AP> f4{&fc2, {0, 0}};

  Boundary<StateAP, 2> b1{100 * m.KT(), {{0, 12e-3}, {23e-3, 0}}, &m};
  Boundary<StateAP, 2> b2{
      100 * m.KT(), {{10e-3, 6e-3}, {0, 0}}, &m, CIRCULAR};
  double a{4e-3};
  Boundary<AP, 2> b3{10*m.KT(), {{-a, a}, {a, -a}}, &m};
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
TEST_F(BoundaryTest, DISABLED_ColloMat_Circular) {
  MatrixXcd c = b2.ColloMatT();
  MatrixXcd d = PseudoInverse(c);
  MatrixXcd m = b2.ModeMatT(f1.Node());

  VectorXcd bv_in = in1.EffectBv(f1.Node());
  VectorXcd bv_bd = m * d * in1.EffectBv(b2.Node());
  EXPECT_TRUE(ApproxVectRv(bv_in, bv_bd, 1e-4, 0, true));
}
TEST_F(BoundaryTest, Resultant_Interface) {
  // Along the interface
  f1.SetCoeff(f1.CSolve({&in1}));

  VectorXcd ref(f1.NumNode() * 2), cal(f1.NumNode() * 2);
  ref = in1.EffectBv(f1.Node());
  for (size_t i = 0; i < f1.NumNode(); i++)
    cal.segment<2>(2 * i) = f1.Pseudo(f1.Node(i)).Bv();

  EXPECT_TRUE(ApproxVectRv(ref, cal, 1e-9));
}
TEST_F(BoundaryTest, Resultant_Boundary) {
  // Along the outer boundary
  f4.SetCoeff(f4.CSolve({&in1}));

  VectorXcd ref(b3.NumNode()*2), cal(b3.NumNode()*2);
  ref = in1.EffectBv(b3.Node());
  for (size_t i = 0; i < b3.NumNode(); i++)
    cal.segment<2>(2 * i) = f4.Pseudo(b3.Node(i)).Bv();

  EXPECT_TRUE(ApproxVectRv(ref, cal, 1e-5, 0, true));
  EXPECT_TRUE(ApproxVectRv(ref, cal, 1e-6, 0, true));
  EXPECT_TRUE(ApproxVectRv(ref, cal, 1e-7, 0, true));
  EXPECT_TRUE(ApproxVectRv(ref, cal, 1e-8, 0, true)); // Fail

  //std::cout << f4.Radius() * m.KT() << std::endl;
}

class AssemBoundaryTest : public Test {
 protected:
  AssemBoundaryTest() : Test(__FILE__, "boundary") {}

  input::Solution s{path("input.txt")};
  Matrix matrix{s};
  AssemblyConfig<AP> c1{s.config(), &matrix};
  AssemblyConfig<AP> c2{s.config(), &matrix};
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

}  // namespace test

}  // namespace mss
