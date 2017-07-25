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

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include "../../src/incident/IncidentInput.h"
#include "../../src/inhomo/ConfigAssembly.h"
#include "../../src/pre/Input.h"

namespace mss {

namespace test {

class AssemblyTest : public testing::Test {
 protected:
  input::Solution s1{testDataPath(__FILE__) + "Multiple.txt"};
  input::Solution s2{testDataPath(__FILE__) + "Single.txt"};
  Matrix matrix1{s1};
  Matrix matrix2{s2};
  ConfigAssembly<StateAP> c1{"c1", s1.config(), &matrix1};
  ConfigAssembly<StateAP> c2{"c2", s2.config(), &matrix2};
};

void AssemblyTest_ReadFile(const std::string& fn, Eigen::VectorXcd& sc) {
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
TEST_F(AssemblyTest, SingleScattering) {
  IncidentPlaneSH inSH(matrix2, s2.incident()[0]);
  c2.Solve({&inSH});
  EXPECT_EQ(c2.Inhomo()[0]->Solve({&inSH}), c2.Inhomo()[0]->ScatterCoeff());

  Eigen::MatrixXcd a;
}

}  // namespace test

}  // namespace mss
