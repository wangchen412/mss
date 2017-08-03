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
  OutputTest() { s.Solve(); }

  SolutionAP s{path("input.txt")};
};

TEST_F(OutputTest, Constructor) {
  post::OutputAP ou(&s);
  EXPECT_EQ(ou.Geo().size(), 4);
}

}  // namespace test

}  // namespace mss
