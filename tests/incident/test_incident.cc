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

// Test plane incident wave classes.

#include "../test.h"

namespace mss {

namespace test {

class IncidentTest : public testing::Test {
 protected:
  IncidentTest()
      : m(Material(1300, 1.41908e9, 0.832e9), 1.25664e6),
        p1(m),
        p2(m, 0, 1, 2),
        p3(m, pi / 3),
        v1(m),
        v2(m, 0, 1, 2),
        v3(m, pi / 3),
        h1(m),
        h2(m, 0, 1, 2),
        h3(m, pi / 3) {}

  Matrix m;
  IncidentPlaneP p1, p2, p3;
  IncidentPlaneSV v1, v2, v3;
  IncidentPlaneSH h1, h2, h3;
};

template <typename T>
void IncidentTest_ReadFile(const Incident<T>* inc,
                           const std::string& fileName, std::vector<T>& ref,
                           std::vector<T>& com) {
  std::ifstream file(testDataPath(__FILE__) + fileName);

  std::string ts;
  while (std::getline(file, ts)) {
    std::stringstream tss(ts);
    PosiVect r;
    T s;
    tss >> r >> s;
    ref.emplace_back(s);
    com.emplace_back(inc->Effect(r));
  }
  file.close();
}

TEST_F(IncidentTest, Constructors) {
  EXPECT_EQ(p1.Angle(), 0);
  EXPECT_EQ(p1.Amplitude(), 1);
  EXPECT_EQ(p1.Phase(), 0);
  EXPECT_EQ(p2.Angle(), 0);
  EXPECT_EQ(p2.Amplitude(), 1);
  EXPECT_EQ(p2.Phase(), 2);
  EXPECT_EQ(p3.Angle(), pi / 3);
  EXPECT_EQ(p3.Amplitude(), 1);
  EXPECT_EQ(p3.Phase(), 0);
  EXPECT_EQ(v1.Angle(), 0);
  EXPECT_EQ(v1.Amplitude(), 1);
  EXPECT_EQ(v1.Phase(), 0);
  EXPECT_EQ(v2.Angle(), 0);
  EXPECT_EQ(v2.Amplitude(), 1);
  EXPECT_EQ(v2.Phase(), 2);
  EXPECT_EQ(v3.Angle(), pi / 3);
  EXPECT_EQ(v3.Amplitude(), 1);
  EXPECT_EQ(v3.Phase(), 0);
  EXPECT_EQ(h1.Angle(), 0);
  EXPECT_EQ(h1.Amplitude(), 1);
  EXPECT_EQ(h1.Phase(), 0);
  EXPECT_EQ(h2.Angle(), 0);
  EXPECT_EQ(h2.Amplitude(), 1);
  EXPECT_EQ(h2.Phase(), 2);
  EXPECT_EQ(h3.Angle(), pi / 3);
  EXPECT_EQ(h3.Amplitude(), 1);
  EXPECT_EQ(h3.Phase(), 0);
}
TEST_F(IncidentTest, EffectSH) {
  std::vector<StateAP> ref1, ref2, ref3;
  std::vector<StateAP> com1, com2, com3;

  IncidentTest_ReadFile(&h1, "SH1.dat", ref1, com1);
  IncidentTest_ReadFile(&h2, "SH2.dat", ref2, com2);
  IncidentTest_ReadFile(&h3, "SH3.dat", ref3, com3);

  EXPECT_EQ(com1.size(), 401);
  EXPECT_EQ(com2.size(), 401);
  EXPECT_EQ(com3.size(), 401);

  EXPECT_THAT(com1, testing::ContainerEq(ref1));
  EXPECT_THAT(com2, testing::ContainerEq(ref2));
  EXPECT_THAT(com3, testing::ContainerEq(ref3));
}
TEST_F(IncidentTest, EffectP) {
  std::vector<StateIP> ref1, ref2, ref3;
  std::vector<StateIP> com1, com2, com3;

  IncidentTest_ReadFile(&p1, "P1.dat", ref1, com1);
  IncidentTest_ReadFile(&p2, "P2.dat", ref2, com2);
  IncidentTest_ReadFile(&p3, "P3.dat", ref3, com3);

  EXPECT_EQ(com1.size(), 401);
  EXPECT_EQ(com2.size(), 401);
  EXPECT_EQ(com3.size(), 401);

  EXPECT_THAT(com1, testing::ContainerEq(ref1));
  EXPECT_THAT(com2, testing::ContainerEq(ref2));
  EXPECT_THAT(com3, testing::ContainerEq(ref3));
}

TEST_F(IncidentTest, EffectSV) {
  std::vector<StateIP> ref1, ref2, ref3;
  std::vector<StateIP> com1, com2, com3;

  IncidentTest_ReadFile(&v1, "SV1.dat", ref1, com1);
  IncidentTest_ReadFile(&v2, "SV2.dat", ref2, com2);
  IncidentTest_ReadFile(&v3, "SV3.dat", ref3, com3);

  EXPECT_EQ(com1.size(), 401);
  EXPECT_EQ(com2.size(), 401);
  EXPECT_EQ(com3.size(), 401);

  EXPECT_THAT(com1, testing::ContainerEq(ref1));
  EXPECT_THAT(com2, testing::ContainerEq(ref2));
  EXPECT_THAT(com3, testing::ContainerEq(ref3));
}

}  // namespace test

}  // namespace mss
