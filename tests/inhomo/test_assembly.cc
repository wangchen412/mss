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

class AssemblyTest : public Test {
 protected:
  AssemblyTest() : Test(__FILE__, "assembly") {}

  input::Solution s{path("input2.txt")};
  Matrix matrix{s};
  AssemblyConfig<AP> c1{s.config(), &matrix};
  AssemblyConfig<AP> c2{s.assembly_config()[1], &matrix};
  AssemblyConfig<AP> c3{s.assembly_config()[2], &matrix};
  AssemblyConfig<AP> c4{s.assembly_config()[3], &matrix};
  AssemblyConfig<AP> c5{s.assembly_config()[4], &matrix};
  IncidentPlaneSH inSH1{matrix, s.incident()[0]};

  Assembly<AP> a1{&c1};
  Assembly<AP> a2{&c1, {40e-3, 30e-3}, pi / 6};

  std::vector<Eigen::Vector2d> f{
      {20e-3, 15e-3}, {60e-3, 15e-3}, {20e-3, 45e-3}, {60e-3, 45e-3}};
};

TEST_F(AssemblyTest, Constructor) {
  EXPECT_EQ(c1.Node().size(), 8794);
  EXPECT_EQ(c1.Node().size(), a1.Node().size());
  EXPECT_EQ(c1.Node().size(), a2.Node().size());

  double a = pi / 6;
  Eigen::Matrix2d rot;
  rot << cos(a), -sin(a), sin(a), cos(a);
  Eigen::Vector2d r0(40e-3, 30e-3);

  for (int i = 0; i < 4; i++)
    EXPECT_EQ(a2.inhomo(i)->PositionGLB(), PosiVect(r0 + rot * f[i]));
}
TEST_F(AssemblyTest, Contain) {
  PosiVect cp(0.0792820323027551, 0.1219615242270663);
  CS cs1(cp + PosiVect(0, -1e-10));
  CS cs2(cp + PosiVect(0, 1e-10));
  EXPECT_TRUE(a2.Contain(&cs1));
  EXPECT_FALSE(a2.Contain(&cs2));
}
TEST_F(AssemblyTest, Solve) {
  c1.DSolve({&inSH1});
  VectorXcd sol1 = a1.DSolve({&inSH1});
  for (int i = 0; i < 4; i++) {
    VectorXcd rr = c1.inhomo(i)->ScatterCoeff();
    VectorXcd cc = sol1.segment(61 * i, 61);
    EXPECT_TRUE(ApproxVectRv(rr, cc, 1e-4, 15));
  }
}
TEST_F(AssemblyTest, Scatter) {
  // To create a fiber and use its nodes as sample points. R = 5e-3.
  FiberConfig<AP> fc1{s.fiber_config()[1], &matrix};
  Fiber<AP> f1{&fc1, {-10e-3, -10e-3}};

  c2.DSolve({&inSH1});
  a2.SetCoeff(a2.DSolve({&inSH1}));

  for (auto& i : f1.Node())
    EXPECT_TRUE(c2.Resultant(i, {&inSH1})
                    .isApprox(a2.Scatter(i) + inSH1.Effect(i), 1e-6));
}
TEST_F(AssemblyTest, AssemblyInAssembly) {
  auto aa1 = dynamic_cast<const Assembly<AP>*>(c3.inhomo(0));
  auto ff1 = dynamic_cast<const Fiber<AP>*>(aa1->inhomo(0));
  EXPECT_EQ(ff1->Position(), PosiVect(20e-3, 15e-3));
  EXPECT_EQ(ff1->Config()->ID(), "b-s");
  EXPECT_EQ(ff1->Config()->Material().Mu(), 8.43e9);
  EXPECT_EQ(ff1->Config()->Material_m().Mu(), 0.832e9);
  EXPECT_EQ(ff1->Config()->Radius(), 10e-3);

  // To create a fiber and use its nodes as sample points.R = 5e-3.
  FiberConfig<AP> fc1{s.fiber_config()[1], &matrix};
  Fiber<AP> f1{&fc1, {-10e-3, -10e-3}};

  c2.DSolve({&inSH1});
  c3.DSolve({&inSH1});

  for (auto& i : f1.Node())
    EXPECT_TRUE(
        c2.Resultant(i, {&inSH1}).isApprox(c3.Resultant(i, {&inSH1}), 1e-6));
}
TEST_F(AssemblyTest, Mixed_ctor) {
  EXPECT_EQ(c4.inhomo().size(), 4);
  auto aa1 = dynamic_cast<const Assembly<AP>*>(c4.inhomo(2));
  auto aa2 = dynamic_cast<const Assembly<AP>*>(c4.inhomo(3));

  double ang1 = 0.5, ang2 = -1;
  Eigen::Matrix2d rot1, rot2;
  rot1 << cos(ang1), -sin(ang1), sin(ang1), cos(ang1);
  rot2 << cos(ang2), -sin(ang2), sin(ang2), cos(ang2);
  Eigen::Vector2d r1(40e-3, 30e-3);
  Eigen::Vector2d r2(20e-3, -40e-3);

  for (int i = 0; i < 4; i++) {
    EXPECT_EQ(aa1->inhomo(i)->PositionGLB(), PosiVect(r1 + rot1 * f[i]));
    EXPECT_EQ(aa1->inhomo(i)->PositionGLB(), c5.inhomo(i)->PositionGLB());
  }
  for (int i = 0; i < 4; i++) {
    EXPECT_EQ(aa2->inhomo(i)->PositionGLB(), PosiVect(r2 + rot2 * f[i]));
    EXPECT_EQ(aa2->inhomo(i)->PositionGLB(), c5.inhomo(i + 4)->PositionGLB());
  }
}
TEST_F(AssemblyTest, Mixed_Scatter) {
  // To create a fiber and use its nodes as sample points. R = 5e-3.
  FiberConfig<AP> fc1{s.fiber_config()[1], &matrix};
  Fiber<AP> f1{&fc1, {-50e-3, -50e-3}};

  c4.DSolve({&inSH1});
  c5.DSolve({&inSH1});

  for (auto& i : f1.Node())
    EXPECT_TRUE(c4.Resultant(i, {&inSH1})
                    .isApprox(c5.Resultant(i, {&inSH1}), 1e-6, true));
}

}  // namespace test

}  // namespace mss
