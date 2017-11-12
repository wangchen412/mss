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

#include "../test.h"

namespace mss {

namespace test {

typedef input::Solution* CreateInputFunc();

input::Solution* ReadInput() {
  std::string path(__FILE__);
  std::string fn = path.substr(0, path.find_last_of("/\\")) + "/../data/input.txt";
  return new input::Solution(fn);
}

input::Solution* ReadWritten() {
  input::Solution* tmp = ReadInput();
  std::ofstream out("InputTest.txt");
  tmp->Print(out);
  out.close();
  delete tmp;
  return new input::Solution("InputTest.txt");
}

class InputTest : public testing::TestWithParam<CreateInputFunc*> {
 public:
  virtual ~InputTest() { delete s_; }
  virtual void SetUp() { s_ = (*GetParam())(); }
  virtual void TearDown() {
    delete s_;
    s_ = nullptr;
    std::remove("InputTest.txt");
  }

 protected:
  input::Solution* s_;
};

TEST_P(InputTest, Constructors) {
  EXPECT_EQ(s_->material().size(), 4);
  EXPECT_EQ(s_->material()[2].lambda, 68.5472e9);
  EXPECT_EQ(s_->frequency(), 1.25664e6);
  EXPECT_EQ(s_->matrix().material->ID, "rubber");
  EXPECT_EQ(s_->incident().size(), 2);
  EXPECT_EQ(s_->fiber_config().size(), 3);
  EXPECT_EQ(s_->fiber_config()[2].material->lambda, 68.5472e9);
  EXPECT_EQ(s_->assembly_config().size(), 2);
  EXPECT_EQ(s_->config().fiber.size(), 5);
  EXPECT_EQ(s_->config().fiber[2].config->radius, 8e-3);
  EXPECT_EQ(s_->config().width, 80e-3);
  EXPECT_EQ(s_->config().height, 60e-3);
  EXPECT_EQ(s_->assembly_config()[1].assembly[0].config->fiber.size(), 5);
  EXPECT_EQ(s_->assembly_config()[1].fiber[0].config->radius, 10e-3);
}

INSTANTIATE_TEST_CASE_P(ioCtorsTest, InputTest,
                        testing::Values(&ReadInput, &ReadWritten));

}  // namespace test

}  // namespace mss
