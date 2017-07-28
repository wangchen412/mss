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

class SolutionSingleTest : public testing::Test {
 protected:
  ~SolutionSingleTest() {
    for (auto& i : sampleCS_) delete i;
  }
  input::Solution in{testDataPath(__FILE__) + "/solution/Single.txt"};
  Solution<StateAP> s{in};
  CSCPtrs sampleCS_;
};

template <typename T>
void SolutionTest_ReadSample(const Solution<T>& s, const std::string& fn,
                             std::vector<T>& ref, std::vector<T>& com,
                             CSCPtrs& sampleCS) {
  std::ifstream file(testDataPath(__FILE__) + fn);
  std::string ts;
  while (std::getline(file, ts)) {
    std::stringstream tss(ts);
    PosiVect r;
    T st;
    tss >> r >> st;
    sampleCS.push_back(new CS(r));
    ref.emplace_back(st);
    com.emplace_back(s.Resultant(sampleCS.back()));
  }
  file.close();
}

void SolutionTest_ReadCoeff(const std::string& fn, Eigen::VectorXcd& sc) {
  std::ifstream file(testDataPath(__FILE__) + fn);
  for (int i = 0; i < sc.size(); i++) file >> sc(i);
  file.close();
}

TEST_F(SolutionSingleTest, Constructor) {
  EXPECT_EQ(s.Incident().size(), 2);
  EXPECT_EQ(s.Incident()[0]->Amplitude(), 2.3e-6);
  EXPECT_EQ(s.Incident()[1]->Amplitude(), 3.2e-6);
  EXPECT_EQ(s.Incident()[0]->Phase(), 3.4);
  EXPECT_EQ(s.Incident()[1]->Phase(), 4.3);
}
TEST_F(SolutionSingleTest, Coefficient) {
  s.Solve();
  Eigen::VectorXcd ref(61);
  SolutionTest_ReadCoeff("/solution/Single_SH1.dat", ref);
  EXPECT_TRUE(
      ApproxVectRV(ref, s.Config().Inhomo()[0]->ScatterCoeff(), 1e-4));
}
TEST_F(SolutionSingleTest, SampleLine) {
  s.Solve();
  std::vector<StateAP> ref, com;
  SolutionTest_ReadSample(s, "/solution/SampleLine_SH1.dat", ref, com,
                          sampleCS_);
  EXPECT_EQ(ref.size(), 1000);
  const double re = 1e-4;
  for (size_t i = 0; i < 1000; i++) EXPECT_TRUE(ref[i].isApprox(com[i], re));
}

}  // namespace test

}  // namespace mss
