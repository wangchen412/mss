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

#include <fstream>
#include <string>
#include "../src/post/Output.h"
#include "../src/post/check/Continuity.h"

using namespace mss;

int main(int argc, char* argv[]) {
  if (argc != 2) exit_error_msg({"Input required."});

  Solution<AP> s{input::Solution(argv[1])};
  delete s.incident_[0];
  s.incident_[0] = new IncidentCylinSH(*s.Matrix(), {-0.3, 0.3});
  s.Solve();

  std::cout << mss_msg({"Wave length: ",
                        std::to_string(2 * pi / s.Matrix()->KT())})
            << std::endl;

  post::CC_Solution<AP> cc{&s};
  std::cout << mss_msg({"Maximum mismatch: ", std::to_string(cc.Max())})
            << std::endl;

  post::Output<AP> oo{&s};
  oo.Write();

  return 0;
}
