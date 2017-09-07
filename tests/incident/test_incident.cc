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

class IncidentTest : public Test {
 protected:
  IncidentTest() : Test(__FILE__) {}
  Matrix m{Material(1300, 1.41908e9, 0.832e9), 1.25664e6};
  IncidentPlaneP p1{m}, p2{m, 0, 1, 2}, p3{m, pi / 3};
  IncidentPlaneSV v1{m}, v2{m, 0, 1, 2}, v3{m, pi / 3};
  IncidentPlaneSH h1{m}, h2{m, 0, 1, 2}, h3{m, pi / 3};

  // Read and compute sample points.
  template <typename T>
  void RnC(const Incident<T>* inc, const std::string& fn, std::vector<T>& ref,
           std::vector<T>& com) {
    ReadSample(fn, ref);
    for (auto& i : SamplePtsBack()) com.emplace_back(inc->Effect(i));
  }
};

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
  std::vector<AP> ref1, ref2, ref3;
  std::vector<AP> com1, com2, com3;

  RnC(&h1, "SH1.dat", ref1, com1);
  RnC(&h2, "SH2.dat", ref2, com2);
  RnC(&h3, "SH3.dat", ref3, com3);

  EXPECT_EQ(com1.size(), 401);
  EXPECT_EQ(com2.size(), 401);
  EXPECT_EQ(com3.size(), 401);

  EXPECT_THAT(com1, testing::ContainerEq(ref1));
  EXPECT_THAT(com2, testing::ContainerEq(ref2));
  EXPECT_THAT(com3, testing::ContainerEq(ref3));
}
TEST_F(IncidentTest, EffectP) {
  std::vector<IP> ref1, ref2, ref3;
  std::vector<IP> com1, com2, com3;

  RnC(&p1, "P1.dat", ref1, com1);
  RnC(&p2, "P2.dat", ref2, com2);
  RnC(&p3, "P3.dat", ref3, com3);

  EXPECT_EQ(com1.size(), 401);
  EXPECT_EQ(com2.size(), 401);
  EXPECT_EQ(com3.size(), 401);

  EXPECT_THAT(com1, testing::ContainerEq(ref1));
  EXPECT_THAT(com2, testing::ContainerEq(ref2));
  EXPECT_THAT(com3, testing::ContainerEq(ref3));
}
TEST_F(IncidentTest, EffectSV) {
  std::vector<IP> ref1, ref2, ref3;
  std::vector<IP> com1, com2, com3;

  RnC(&v1, "SV1.dat", ref1, com1);
  RnC(&v2, "SV2.dat", ref2, com2);
  RnC(&v3, "SV3.dat", ref3, com3);

  EXPECT_EQ(com1.size(), 401);
  EXPECT_EQ(com2.size(), 401);
  EXPECT_EQ(com3.size(), 401);

  EXPECT_THAT(com1, testing::ContainerEq(ref1));
  EXPECT_THAT(com2, testing::ContainerEq(ref2));
  EXPECT_THAT(com3, testing::ContainerEq(ref3));
}

}  // namespace test

}  // namespace mss
