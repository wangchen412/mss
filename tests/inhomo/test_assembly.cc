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
  AssemblyTest() : Test(__FILE__, "assembly") {
    for (auto& i : c.Inhomo())
      fiberPtrs.push_back(dynamic_cast<Fiber<StateAP>*>(i));
  }

  input::Solution s{path("input.txt")};
  Matrix matrix{s};
  ConfigAssembly<StateAP> c{"Test", s.config(), &matrix};
  IncidentPlaneSH inSH1{matrix, s.incident()[0]};
  IncidentPlaneSH inSH2{matrix, s.incident()[1]};
  InciCPtrs<StateAP> incident{&inSH1, &inSH2};
  std::vector<const Fiber<StateAP>*> fiberPtrs;
};

TEST_F(AssemblyTest, Constructor) {
  EXPECT_EQ(c.ID(), "Test");
  EXPECT_EQ(c.Inhomo()[0]->Position(), PosiVect(0, 0));
  EXPECT_EQ(c.Inhomo()[1]->Position(), PosiVect(18e-3, 18e-3));
  EXPECT_EQ(c.Inhomo()[2]->Position(), PosiVect(-18e-3, -18e-3));
  EXPECT_EQ(c.Inhomo()[3]->Position(), PosiVect(-18e-3, 18e-3));
  EXPECT_EQ(c.Inhomo()[4]->Position(), PosiVect(18e-3, -18e-3));

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
TEST_F(AssemblyTest, InWhich) {
  CS p1(0, 0), p2(12e-3, 0), p3(17e-3, 17e-3), p4(-11e-3, -11e-3),
      p5(15e-3, -15e-3), p6(50e-3, 50e-3),
      p7(21e-3 - epsilon, 22e-3 - epsilon),
      p8(21e-3 + epsilon, 22e-3 + epsilon);

  EXPECT_EQ(c.InWhich(&p1), c.Inhomo()[0]);
  EXPECT_EQ(c.InWhich(&p2), nullptr);
  EXPECT_EQ(c.InWhich(&p3), c.Inhomo()[1]);
  EXPECT_EQ(c.InWhich(&p4), nullptr);
  EXPECT_EQ(c.InWhich(&p5), c.Inhomo()[4]);
  EXPECT_EQ(c.InWhich(&p6), nullptr);
  EXPECT_EQ(c.InWhich(&p7), c.Inhomo()[1]);
  EXPECT_EQ(c.InWhich(&p8), nullptr);
}
TEST_F(AssemblyTest, Solve) {
  c.Solve(incident);
  Eigen::VectorXcd ref(305);
  ReadCoeff("Coeff_SH.dat", ref);

  for (int i = 0; i < 5; i++) {
    Eigen::VectorXcd rr = ref.segment(61 * i, 61),
                     cc = c.Inhomo()[i]->ScatterCoeff();
    EXPECT_TRUE(ApproxVectRV(rr, cc, 1e-5, 10));
  }
}
TEST_F(AssemblyTest, Scatter) {
  std::vector<StateAP> ref, com;
  ReadSample("line_1.dat", ref);
  EXPECT_EQ(ref.size(), 100);

  c.Solve(incident);
  for (auto& i : SamplePts(0)) com.emplace_back(c.Resultant(i, incident));
  EXPECT_EQ(com.size(), 100);

  for (size_t i = 0; i < 100; i++) EXPECT_TRUE(ref[i].isApprox(com[i], 1e-5));
}

}  // namespace test

}  // namespace mss
