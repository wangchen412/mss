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

class ArrayTest : public Test {
 protected:
  ArrayTest() : Test(__FILE__, "array") {}

  input::Solution ins{path("input.txt")};
  Solution<AP> s{ins};
};

TEST_F(ArrayTest, Constructor) {
  EXPECT_EQ(s.inhomo().size(), 15);
}
TEST_F(ArrayTest, DISABLED_Solve) {
  s.Solve();
  post::OutputAP o{&s};
  o.Write();
}

}  // namespace test

}  // namespace mss
