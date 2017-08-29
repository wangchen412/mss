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

#include "../../src/inhomo/Boundary.h"
#include "../test.h"

namespace mss {

namespace test {

class BoundaryTest : public Test {
 protected:
  Matrix m{Material(1300, 1.41908e9, 0.832e9), 1.25664e6};
  IncidentPlaneSH in1{m, pi / 3, 1, 2};

  Boundary<StateAP> b1{10 * m.KT(), {{0, 12e-3}, {23e-3, 0}}, &m};
};

TEST_F(BoundaryTest, DISABLED_Constructor) {
  EXPECT_EQ(b1.Node().size(), 1098);
}

TEST_F(BoundaryTest, InfMat) {
  CS cs(10e-3, 10e-3, 1);
  Eigen::Vector2cd bv_in = in1.Effect(&cs).BV();
  Eigen::Vector2cd bv_bd = b1.InfMatT(&cs) * in1.EffectBv(b1.Node());

  // std::cout << b1.InfMatT(&cs) << std::endl;
  // std::cout << in1.EffectBv(b1.Node()) << std::endl;

  
  std::cout << bv_in << std::endl;
  std::cout << bv_bd << std::endl;

  // EXPECT_EQ(bv_in, bv_bd);

  // CS cs1(10e-3, 10e-3), cs2(30e-3, 10e-3);
  // std::cout << GreenT<StateAP>(&cs2, &cs1, &m) << std::endl;
  // std::cout << m.KT() << std::endl;

  // std::cout << b1.InfMatT(&cs1);
}

}  // namespace test

}  // namespace mss
