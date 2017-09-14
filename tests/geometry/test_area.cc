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

#include "../../src/post/geometry/Area.h"
#include "../test.h"

namespace mss {

namespace test {

class AreaTest : public Test {
 protected:
  AreaTest() : Test(__FILE__) {}
  SolutionAP s{path("input.txt")};
};

TEST_F(AreaTest, Computation) {
  std::vector<AP> ref, com;
  ReadSample("Area_r1.dat", ref);
  EXPECT_EQ(ref.size(), 100);

  s.Solve();

  post::AreaAP a1{&s, {1, 2}, {3, 4}, 10, 10};
  EXPECT_EQ(a1.Points().size(), 100);

  post::AreaAP t1(&s, {-50e-3, 50e-3}, {50e-3, -50e-3}, 10, 10);
  std::string fn("Area_t1.dat");
  std::ofstream file(path(fn));
  t1.Print(file);
  file.close();
  ReadSample(fn, com, 5);
  EXPECT_EQ(com.size(), 100);
  std::remove(path(fn).c_str());

  // The tolerance is set as 1e-3 since the previous program cannot compute
  // the data for the center points accurately.
  for (size_t i = 0; i < 100; i++) EXPECT_TRUE(ref[i].isApprox(com[i], 1e-3));
}

}  // namespace test

}  // namespace mss
