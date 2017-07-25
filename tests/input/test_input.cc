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
#include "../../src/pre/Input.h"

namespace mss {

namespace test {

typedef input::Solution* CreateInputFunc();

input::Solution* ReadInput() {
  return new input::Solution(testDataPath(__FILE__) + "input.txt");
}

input::Solution* ReadWritten() {
  input::Solution s(testDataPath(__FILE__) + "input.txt");
  std::ofstream file(testDataPath(__FILE__) + "output.txt");
  s.Print(file);
  file.close();
  return new input::Solution(testDataPath(__FILE__) +
                             std::string("output.txt"));
}

class InputTest : public testing::TestWithParam<CreateInputFunc*> {
 public:
  virtual ~InputTest() { delete s_; }
  virtual void SetUp() { s_ = (*GetParam())(); }
  virtual void TearDown() {
    delete s_;
    s_ = nullptr;
  }

 protected:
  input::Solution* s_;
};

TEST_P(InputTest, Constructors) {
  EXPECT_EQ(s_->material().size(), 4);
  EXPECT_EQ(s_->material()[2].lambda, 68.5472e9);
  EXPECT_EQ(s_->frequency(), 1.23);
  EXPECT_EQ(s_->matrix().material->ID, "rubber");
  EXPECT_EQ(s_->incident().size(), 2);
  EXPECT_EQ(s_->configFiber().size(), 3);
  EXPECT_EQ(s_->configFiber()[2].material->lambda, 68.5472e9);
  EXPECT_EQ(s_->config().fiber.size(), 5);
  EXPECT_EQ(s_->config().fiber[2].config->radius, 8e-3);
}

INSTANTIATE_TEST_CASE_P(ioCtorsTest, InputTest,
                        testing::Values(&ReadInput, &ReadWritten));

}  // namespace test

}  // namespace mss
