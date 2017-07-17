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

void operator>>(istream& is, Material& m) {
  is >> m.ID >> m.rho >> m.mu >> m.lambda;
}
ostream& operator<<(ostream& os, const Material& m) {
  return os << m.ID << "\t" << m.rho << "\t" << m.mu << "\t" << m.lambda;
}

void operator>>(istream& is, Matrix& m) {
  is >> m.frequency >> m.delta >> m.materialID;
}
ostream& operator<<(ostream& os, const Matrix& m) {
  return os << m.frequency << "\t" << m.delta << "\t" << m.materialID;
}

void operator>>(istream& is, IncidentPlane& i) {
  is >> i.type >> i.angle >> i.amplitude >> i.phase;
}
ostream& operator<<(ostream& os, const IncidentPlane& i) {
  return os << i.type << "\t" << i.angle << "\t" << i.amplitude << "\t"
            << i.phase;
}

void operator>>(istream& is, ConfigFiber& c) {
  is >> c.ID >> c.radius >> c.N_max >> c.materialID;
}
ostream& operator<<(ostream& os, const ConfigFiber& c) {
  return os << c.ID << "\t" << c.radius << "\t" << c.N_max << "\t"
            << c.materialID;
}

void operator>>(istream& is, Fiber& f) {
  is >> f.position >> f.configID;
}
ostream& operator<<(ostream& os, const Fiber& f) {
  return os << f.position << "\t" << f.configID;
}

Solution::Solution(const string& fn) : fn_(fn) {
  add(material_, "[Materials]");
  add(matrix_, "[Matrix]");
  add(incident_, "[Incident Waves]");
  add(configFiber_, "[Fiber Configurations]");
  add(fiber_, "[Fibers]");
}

template <typename T>
void Solution::add(vector<T>& vec, const string& key) {
  ifstream file(fn_);
  skipUntil(file, key, 2);
  string tmp;
  while (getline(file, tmp)) {
    if (isWhiteSpace(tmp)) break;
    if (tmp[0] == '#') continue;
    vec.emplace_back();
    stringstream(tmp) >> vec.back();
    //    vec.emplace_back(stringstream(tmp));
  }
  file.close();
}

}  // namespace input

}  // namespace mss
