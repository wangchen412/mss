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

int main() {
  Eigen::VectorXcd w(1000), t(1000);
  ReadBv(w, t);

  Mismatch f(16576.2, w, t, {{8550, 356.25}, 0, {39515625000, 2634375000}});
  Eigen::Vector4d x0{2, 10, 2, 10};
  Eigen::Vector4d mm = GradientDescent(f, x0);
  return 0;
}
