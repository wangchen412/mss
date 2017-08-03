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

class SolutionTest : public Test {
 protected:
  SolutionTest() : Test(__FILE__) {}

  input::Solution inS{path("Single.txt")};
  Solution<StateAP> s{inS};
};

TEST_F(SolutionTest, Constructor) {
  EXPECT_EQ(s.Incident().size(), 2);
  EXPECT_EQ(s.Incident()[0]->Amplitude(), 2.3e-6);
  EXPECT_EQ(s.Incident()[1]->Amplitude(), 3.2e-6);
  EXPECT_EQ(s.Incident()[0]->Phase(), 3.4);
  EXPECT_EQ(s.Incident()[1]->Phase(), 4.3);
}
TEST_F(SolutionTest, Coefficient) {
  Eigen::VectorXcd ref(61);
  ReadCoeff("Single_SH1.dat", ref);

  s.Solve();
  EXPECT_TRUE(
      ApproxVectRV(ref, s.Config().Inhomo()[0]->ScatterCoeff(), 1e-4));
}
TEST_F(SolutionTest, SampleLine) {
  std::vector<StateAP> ref, com;
  ReadSample("SampleLine_SH1.dat", ref);
  EXPECT_EQ(ref.size(), 100);

  s.Solve();
  for (auto& i : SamplePts(0)) com.emplace_back(s.Resultant(i));
  EXPECT_EQ(com.size(), 100);

  const double re = 1e-4;
  for (size_t i = 0; i < 100; i++) EXPECT_TRUE(ref[i].isApprox(com[i], re));
}
TEST_F(SolutionTest, MsSampleLine) {
  std::vector<StateAP> ref, com;
  ReadSample("line_1.dat", ref);
  EXPECT_EQ(ref.size(), 100);

  Solution<StateAP> ss{input::Solution(path("Multiple.txt"))};
  ss.Solve();
  for (auto& i : SamplePts(0)) com.emplace_back(ss.Resultant(i));
  EXPECT_EQ(com.size(), 100);

  for (size_t i = 0; i < 100; i++) EXPECT_TRUE(ref[i].isApprox(com[i], 1e-4));
}

}  // namespace test

}  // namespace mss
