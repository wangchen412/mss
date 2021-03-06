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

class AssemConfigTest : public Test {
 protected:
  AssemConfigTest() : Test(__FILE__, "assembly") {
    for (auto& i : c.inhomo())
      fiberPtrs.push_back(dynamic_cast<const Fiber<AP>*>(i));
  }

  input::Solution s{path("input.txt")};
  Matrix matrix{s};
  AssemblyConfig<AP> c{s.config(), &matrix};
  IncidentPlaneSH inSH1{matrix, s.incident()[0]};
  IncidentPlaneSH inSH2{matrix, s.incident()[1]};
  InciCPtrs<AP> incident{&inSH1, &inSH2};
  std::vector<const Fiber<AP>*> fiberPtrs;
};

TEST_F(AssemConfigTest, Constructor) {
  EXPECT_EQ(c.ID(), "Assembly_1");
  EXPECT_EQ(c.inhomo(0)->Position(), PosiVect(0, 0));
  EXPECT_EQ(c.inhomo(1)->Position(), PosiVect(18e-3, 18e-3));
  EXPECT_EQ(c.inhomo(2)->Position(), PosiVect(-18e-3, -18e-3));
  EXPECT_EQ(c.inhomo(3)->Position(), PosiVect(-18e-3, 18e-3));
  EXPECT_EQ(c.inhomo(4)->Position(), PosiVect(18e-3, -18e-3));

  EXPECT_EQ(fiberPtrs[0]->Config()->ID(), "b-s");
  EXPECT_EQ(fiberPtrs[1]->Config()->ID(), "s-h");
  EXPECT_EQ(fiberPtrs[2]->Config()->ID(), "m-l");
  EXPECT_EQ(fiberPtrs[3]->Config()->ID(), "m-l");
  EXPECT_EQ(fiberPtrs[4]->Config()->ID(), "b-s");

  EXPECT_EQ(fiberPtrs[0]->Config()->Material().Mu(), 8.43e9);
  EXPECT_EQ(fiberPtrs[1]->Config()->Material().Mu(), 80.0698e9);
  EXPECT_EQ(fiberPtrs[2]->Config()->Material().Mu(), 28.65845e9);
  EXPECT_EQ(fiberPtrs[3]->Config()->Material().Mu(), 28.65845e9);
  EXPECT_EQ(fiberPtrs[4]->Config()->Material().Mu(), 8.43e9);

  EXPECT_EQ(fiberPtrs[0]->Config()->Material_m().Mu(), 0.832e9);
  EXPECT_EQ(fiberPtrs[1]->Config()->Material_m().Mu(), 0.832e9);
  EXPECT_EQ(fiberPtrs[2]->Config()->Material_m().Mu(), 0.832e9);
  EXPECT_EQ(fiberPtrs[3]->Config()->Material_m().Mu(), 0.832e9);
  EXPECT_EQ(fiberPtrs[4]->Config()->Material_m().Mu(), 0.832e9);

  EXPECT_EQ(fiberPtrs[0]->Config()->Radius(), 10e-3);
  EXPECT_EQ(fiberPtrs[1]->Config()->Radius(), 5e-3);
  EXPECT_EQ(fiberPtrs[2]->Config()->Radius(), 8e-3);
  EXPECT_EQ(fiberPtrs[3]->Config()->Radius(), 8e-3);
  EXPECT_EQ(fiberPtrs[4]->Config()->Radius(), 10e-3);
}
TEST_F(AssemConfigTest, InWhich) {
  CS p1(0, 0), p2(12e-3, 0), p3(17e-3, 17e-3), p4(-11e-3, -11e-3),
      p5(15e-3, -15e-3), p6(50e-3, 50e-3),
      p7(21e-3 - epsilon, 22e-3 - epsilon),
      p8(21e-3 + epsilon, 22e-3 + epsilon);

  EXPECT_EQ(c.InWhich(&p1), c.inhomo(0));
  EXPECT_EQ(c.InWhich(&p2), nullptr);
  EXPECT_EQ(c.InWhich(&p3), c.inhomo(1));
  EXPECT_EQ(c.InWhich(&p4), nullptr);
  EXPECT_EQ(c.InWhich(&p5), c.inhomo(4));
  EXPECT_EQ(c.InWhich(&p6), nullptr);
  EXPECT_EQ(c.InWhich(&p7), c.inhomo(1));
  EXPECT_EQ(c.InWhich(&p8), nullptr);
}
TEST_F(AssemConfigTest, Solve) {
  VectorXcd ref(305);
  ReadCoeff("Coeff_SH.dat", ref);

  c.Solve(incident, DFT);
  for (int i = 0; i < 5; i++) {
    VectorXcd rr = ref.segment(61 * i, 61), cc = c.inhomo(i)->ScatterCoeff();
    EXPECT_TRUE(ApproxVectRv(rr, cc, 1e-3, 10));
  }
  c.Solve(incident, COLLOCATION);
  for (int i = 0; i < 5; i++) {
    VectorXcd rr = ref.segment(61 * i, 61), cc = c.inhomo(i)->ScatterCoeff();
    EXPECT_TRUE(ApproxVectRv(rr, cc, 1e-5, 10));
  }
}
TEST_F(AssemConfigTest, Scatter) {
  std::vector<AP> ref, com1, com2;
  ReadSample("line_1.dat", ref);
  EXPECT_EQ(ref.size(), 100);

  c.Solve(incident, DFT);
  for (auto& i : SamplePts()) com1.emplace_back(c.Resultant(i, incident));
  EXPECT_EQ(com1.size(), 100);
  for (size_t i = 0; i < 100; i++)
    EXPECT_TRUE(ref[i].isApprox(com1[i], 1e-3));

  c.Solve(incident, COLLOCATION);
  for (auto& i : SamplePts()) com2.emplace_back(c.Resultant(i, incident));
  EXPECT_EQ(com2.size(), 100);
  for (size_t i = 0; i < 100; i++)
    EXPECT_TRUE(ref[i].isApprox(com2[i], 1e-5));
}

// class AssemblyTest : public AssemConfigTest {
//  protected:
//   AssemblyConfig<AP> c1{s.assembly_config()[1], &matrix};
//   Assembly<AP> a1{&c1};
//   Assembly<AP> a2{&c1, {40e-3, 30e-3}, pi / 6};
// };

}  // namespace test

}  // namespace mss
