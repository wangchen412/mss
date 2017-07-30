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

// Unit test of Matrix class.

#include "../test.h"

namespace mss {

namespace test {

class MatrixTest : public Test {
 protected:
  const Material rubber{1300, 1.41908e9, 0.832e9},
      lead{11400, 36.32496e9, 8.43e9};
  Matrix a{rubber, pi / 2}, b{lead, pi};
};

TEST_F(MatrixTest, Constructors) {
  EXPECT_DOUBLE_EQ(a.Frequency(), pi / 2);
  EXPECT_DOUBLE_EQ(a.KL(), 0.001019997614801881);
  EXPECT_DOUBLE_EQ(a.KT(), 0.0019634954084936209);
  EXPECT_DOUBLE_EQ(b.Frequency(), pi);
  EXPECT_DOUBLE_EQ(b.KL(), 0.001454480422257349);
  EXPECT_DOUBLE_EQ(b.KT(), 0.0036533267014104112);
}

// TODO: Add tests for constitutive relation methods.

}  // namespace test

}  // namespace mss
