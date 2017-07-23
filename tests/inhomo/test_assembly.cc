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

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include "../../src/incident/IncidentInput.h"
#include "../../src/inhomo/ConfigAssembly.h"

namespace mss {

namespace test {

class AssemblyTest : public testing::Test {
 protected:
  const double omega = 1.25664e6;
  Material rubber = {1300, 1.41908e9, 0.832e9};
  Material lead = {11400, 36.32496e9, 8.43e9};
  Matrix matrix = {rubber, omega};

  ConfigFiber<StateIP> c1 = {"c1", 30, 300, 1e-3, lead, &matrix};
  ConfigFiber<StateAP> c2 = {"c2", 30, 300, 1e-3, lead, &matrix};
};

}  // namespace test

}  // namespace mss
