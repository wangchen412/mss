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

// Test input classes.

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include "../../src/input/Input.h"

namespace mss {

namespace test {

std::string f1 = testDataPath(__FILE__) + std::string("input_material.txt");

class InputTest : public testing::Test {
 protected:
  InputTest() : s(f1) {}

  input::Solution s;
};

TEST_F(InputTest, Constructors) {
  EXPECT_EQ(1, 1);

  for (const auto& i : s.material_) {
    std::cout << i << std::endl;
  }

  std::cout << s.matrix_.size() << std::endl;
  for (auto& i : s.matrix_) {
    std::cout << i.materialID << std::endl;
  }

  for (auto& i : s.incident_) {
    std::cout << i.type << " " << i.angle << " " << i.amplitude << " "
              << i.phase << std::endl;
  }

  for (auto& i : s.configFiber_) {
    std::cout << i.ID << " " << i.radius << " " << i.N_max << " "
              << i.materialID << std::endl;
  }
}

}  // namespace test

}  // namespace mss
