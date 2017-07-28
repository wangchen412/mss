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

class AssemblyTest : public testing::Test {
 protected:
  input::Solution s1{testDataPath(__FILE__) + "assembly/Multiple.txt"};
  Matrix matrix1{s1};
  ConfigAssembly<StateAP> c1{"c1", s1.config(), &matrix1};
  IncidentPlaneSH inSH1{matrix1, s1.incident()[0]};
};

class AssemblySingleTest : public testing::Test {
 protected:
  input::Solution s2{testDataPath(__FILE__) + "assembly/Single.txt"};
  Matrix matrix2{s2};
  ConfigAssembly<StateAP> c2{"c2", s2.config(), &matrix2};
  IncidentPlaneSH inSH2{matrix2, s2.incident()[0]};
};

void AssemblyTest_ReadCoeff(const std::string& fn, Eigen::VectorXcd& sc) {
  std::ifstream file(testDataPath(__FILE__) + fn);
  for (int i = 0; i < sc.size(); i++) file >> sc(i);
  file.close();
}

TEST_F(AssemblyTest, Constructor) {
  EXPECT_EQ(c1.ID(), "c1");
  EXPECT_EQ(c1.Inhomo()[0]->Position(), PosiVect(0, 0));
  EXPECT_EQ(c1.Inhomo()[1]->Position(), PosiVect(18e-3, 18e-3));
  EXPECT_EQ(c1.Inhomo()[2]->Position(), PosiVect(-18e-3, -18e-3));
  EXPECT_EQ(c1.Inhomo()[3]->Position(), PosiVect(-18e-3, 18e-3));
  EXPECT_EQ(c1.Inhomo()[4]->Position(), PosiVect(18e-3, -18e-3));
}
TEST_F(AssemblyTest, InWhich) {
  CS p1(0, 0), p2(12e-3, 0), p3(17e-3, 17e-3), p4(-11e-3, -11e-3),
      p5(15e-3, -15e-3), p6(50e-3, 50e-3),
      p7(21e-3 - epsilon, 22e-3 - epsilon),
      p8(21e-3 + epsilon, 22e-3 + epsilon);

  EXPECT_EQ(c1.InWhich(&p1), c1.Inhomo()[0]);
  EXPECT_EQ(c1.InWhich(&p2), nullptr);
  EXPECT_EQ(c1.InWhich(&p3), c1.Inhomo()[1]);
  EXPECT_EQ(c1.InWhich(&p4), nullptr);
  EXPECT_EQ(c1.InWhich(&p5), c1.Inhomo()[4]);
  EXPECT_EQ(c1.InWhich(&p6), nullptr);
  EXPECT_EQ(c1.InWhich(&p7), c1.Inhomo()[1]);
  EXPECT_EQ(c1.InWhich(&p8), nullptr);
}
TEST_F(AssemblySingleTest, SingleScattering) {
  c2.Solve({&inSH2});
  Eigen::VectorXcd ref(61), wref(61);
  AssemblyTest_ReadCoeff("assembly/Single_SH2.dat", ref);
  AssemblyTest_ReadCoeff("assembly/Single_SH1.dat", wref);

  EXPECT_TRUE(ApproxVectRV(ref, c2.Inhomo()[0]->ScatterCoeff(), 1e-4));
  EXPECT_FALSE(ApproxVectRV(wref, c2.Inhomo()[0]->ScatterCoeff(), 1e-4));
}
TEST_F(AssemblyTest, MultipleScattering) {
  c1.Solve({&inSH1});
  Eigen::VectorXcd ref(305);
  AssemblyTest_ReadCoeff("assembly/Multiple_SH1.dat", ref);

  for (int i = 0; i < 5; i++) {
    Eigen::VectorXcd rr = ref.segment(61 * i, 61),
                     cc = c1.Inhomo()[i]->ScatterCoeff();
    EXPECT_TRUE(ApproxVectRV(rr, cc, 1e-5, 10));
  }
}

}  // namespace test

}  // namespace mss
