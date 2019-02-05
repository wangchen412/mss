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

#include "homo.h"

using namespace mss;

int main(int argc, char* argv[]) {
  if (argc != 9) exit_error_msg({"Material parameters required."});
  dcomp rho1{atof(argv[1]), atof(argv[2])}, mu1{atof(argv[3]), atof(argv[4])};
  dcomp rho2{atof(argv[5]), atof(argv[6])}, mu2{atof(argv[7]), atof(argv[8])};

  Scan({rho1, mu1, mu1}, {rho2, mu2, mu2}, 16);

  return 0;
}
