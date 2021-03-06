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
  Solution<AP> s{inS};
};

TEST_F(SolutionTest, Constructor) {
  EXPECT_EQ(s.Incident().size(), 2);
  EXPECT_EQ(s.Incident()[0]->Amplitude(), 2.3e-6);
  EXPECT_EQ(s.Incident()[1]->Amplitude(), 3.2e-6);
  EXPECT_EQ(s.Incident()[0]->Phase(), 3.4);
  EXPECT_EQ(s.Incident()[1]->Phase(), 4.3);
}
TEST_F(SolutionTest, Coefficient) {
  VectorXcd ref(61);
  ReadCoeff("Single_SH1.dat", ref);

  s.Solve();
  EXPECT_TRUE(ApproxVectRv(ref, s.Config().inhomo(0)->ScatterCoeff(), 1e-4));
}
TEST_F(SolutionTest, SampleLine) {
  std::vector<AP> ref, com;
  ReadSample("SampleLine_SH1.dat", ref);
  EXPECT_EQ(ref.size(), 100);

  s.Solve();
  for (auto& i : SamplePts()) com.emplace_back(s.Resultant(i));
  EXPECT_EQ(com.size(), 100);

  const double re = 1e-4;
  for (size_t i = 0; i < 100; i++) EXPECT_TRUE(ref[i].isApprox(com[i], re));
}
TEST_F(SolutionTest, MsSampleLine) {
  std::vector<AP> ref, com1, com2;
  ReadSample("SampleLine_SH2.dat", ref);
  EXPECT_EQ(ref.size(), 100);

  Solution<AP> sc{input::Solution(path("Multiple.txt"))};
  for (auto& i : SamplePts()) com1.emplace_back(sc.Solve().Resultant(i));
  EXPECT_EQ(com1.size(), 100);
  for (size_t i = 0; i < 100; i++)
    EXPECT_TRUE(ref[i].isApprox(com1[i], 1e-4));

  Solution<AP> sd{input::Solution(path("Multiple_DFT.txt"))};
  for (auto& i : SamplePts()) com2.emplace_back(sd.Solve().Resultant(i));
  EXPECT_EQ(com2.size(), 100);
  for (size_t i = 0; i < 100; i++)
    EXPECT_TRUE(ref[i].isApprox(com2[i], 1e-3));
}

}  // namespace test

}  // namespace mss
