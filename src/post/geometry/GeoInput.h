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

#ifndef MSS_GEOINPUT_H
#define MSS_GEOINPUT_H

#include "Area.h"
#include "Circle.h"
#include "Line.h"
#include "Point.h"

namespace mss {

namespace post {

template <typename T>
class GeoInput {
 public:
  GeoInput(const Solution<T>* solution) : sol_(solution) { add_map(); }
  Geometry<T>* operator()(std::stringstream& ss) {
    std::string type;
    ss >> type;
    return funcMap[type](ss);
  }

 private:
  typedef std::function<Geometry<T>*(std::stringstream&)> funcType;
  std::map<std::string, funcType, ci_comp> funcMap;
  const Solution<T>* sol_;

  void add_map();
};

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T>
void GeoInput<T>::add_map() {
  funcMap["Point"] = [this](std::stringstream& ss) {
    std::string ID;
    PosiVect position;
    double angle;
    ss >> ID >> position >> angle;
    return new Point<T>(sol_, position, angle, ID);
  };

  funcMap["Line"] = [this](std::stringstream& ss) {
    std::string ID;
    PosiVect p1, p2;
    size_t N;
    ss >> ID >> p1 >> p2 >> N;
    return new Line<T>(sol_, p1, p2, N, ID);
  };

  funcMap["Circle"] = [this](std::stringstream& ss) {
    std::string ID;
    PosiVect center;
    double R;
    size_t N;
    ss >> ID >> center >> R >> N;
    return new Circle<T>(sol_, center, R, N, ID);
  };

  funcMap["Area"] = [this](std::stringstream& ss) {
    std::string ID;
    PosiVect p1, p2;
    size_t Nx, Ny;
    ss >> ID >> p1 >> p2 >> Nx >> Ny;
    return new Area<T>(sol_, p1, p2, Nx, Ny, ID);
  };
}

}  // namespace post

}  // namespace mss

#endif
