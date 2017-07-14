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

#include "Input.h"

namespace mss {

namespace input {

using namespace std;

Material::Material(stringstream data) {
  data >> ID >> rho >> lambda >> mu;
}

Matrix::Matrix(stringstream data) {
  data >> frequency >> delta >> 
}

Solution::Solution(const string& file) : fn_(file) {
  add(material_);
  add(matrix_);
  add(incident_);
}

template <typename T>
void Solution::add(vector<T*>& vec) {
  ifstream file(fn_);
  skipUntil(file, T::key, 2);
  string tmp;
  while (getline(file, tmp)) {
    if (isWhiteSpace(tmp)) break;
    if (tmp[0] == '#') continue;
    vec.push_back(new T(stringstream(tmp)));
  }
  file.close();
}

}  // namespace input

}  // namespace mss
