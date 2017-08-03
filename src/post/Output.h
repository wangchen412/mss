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

#ifndef MSS_OUTPUT_H
#define MSS_OUTPUT_H

#include "../core/Solution.h"
#include "geometry/GeoInput.h"

namespace mss {

namespace post {

template <typename T>
class Output {
 public:
  Output(const Solution<T>* solution, const std::string& fn)
      : sol_(solution), fn_(fn) {
    add_geo();
  }
  Output(const Solution<T>* solution)
      : Output(solution, solution->InputFN()) {}

  void Write() const;
  const GeoCPtrs<T>& Geo() const { return geo_; }
  const Geometry<T>* Geo(const size_t& sn) const { return geo_[sn]; }

 private:
  const Solution<T>* sol_;
  std::string fn_;
  GeoCPtrs<T> geo_;

  void add_geo();
};

typedef Output<StateIP> OutputIP;
typedef Output<StateAP> OutputAP;

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T>
void Output<T>::add_geo() {
  GeoInput<T> f(sol_);
  std::ifstream file(fn_);
  skipUntil(file, "[Output]", 2);

  std::string tmp;
  while (std::getline(file, tmp)) {
    if (isWhiteSpace(tmp)) break;
    if (tmp[0] == '#') continue;
    std::stringstream ss(tmp);
    geo_.push_back(f(ss));
  }

  file.close();
}

template <typename T>
void Output<T>::Write() const {
  for (auto& i : geo_) i->Write();
}

}  // namespace post

}  // namespace mss

#endif
