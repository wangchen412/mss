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

#include "../../src/post/Output.h"
#include "../test.h"

namespace mss {

namespace test {

class OutputTest : public Test {
 protected:
  ~OutputTest() {
    std::remove("Point_1.dat");
    std::remove("Line_1.dat");
    std::remove("Circle_1.dat");
    std::remove("Area_1.dat");
  }
  SolutionAP s{path("input.txt")};

  template <typename T>
  void extract_state(const post::Geometry<T>* p, std::vector<AP>& c) const {
    const post::PointSet<T>* sp = dynamic_cast<const post::PointSet<T>*>(p);
    for (auto& i : sp->Points()) c.push_back(i->State());
  }
};

TEST_F(OutputTest, Write) {
  s.Solve();

  post::OutputAP o1{&s};
  EXPECT_EQ(o1.Geo().size(), 4);
  EXPECT_EQ(o1.Geo(0)->ID(), "Point_1");
  EXPECT_EQ(o1.Geo(1)->ID(), "Line_1");
  EXPECT_EQ(o1.Geo(2)->ID(), "Circle_1");
  EXPECT_EQ(o1.Geo(3)->ID(), "Area_1");

  o1.Write();
  std::ifstream f1("Line_1.dat"), f2("Circle_1.dat"), f3("Area_1.dat");
  std::vector<AP> r1, r2, r3;
  ReadSample(f1, r1, 5);
  ReadSample(f2, r2, 5);
  ReadSample(f3, r3, 5);

  ASSERT_EQ(r1.size(), 100);
  ASSERT_EQ(r2.size(), 100);
  ASSERT_EQ(r3.size(), 100);

  std::vector<AP> c1, c2, c3;
  extract_state(o1.Geo(1), c1);
  extract_state(o1.Geo(2), c2);
  extract_state(o1.Geo(3), c3);

  ASSERT_EQ(c1.size(), 100);
  ASSERT_EQ(c2.size(), 100);
  ASSERT_EQ(c3.size(), 100);

  EXPECT_THAT(r1, testing::Eq(c1));
  EXPECT_THAT(r2, testing::Eq(c2));
  EXPECT_THAT(r3, testing::Eq(c3));
}

}  // namespace test

}  // namespace mss
