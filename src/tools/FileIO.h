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

#ifndef MSS_FILEIO_H
#define MSS_FILEIO_H

#include <fstream>
#include <iostream>
#include <string>

namespace mss {

// Case insensitive comparison of two std::strings.
struct ci_comp {
  bool operator()(const std::string &lhs, const std::string &rhs) const {
    return strcasecmp(lhs.c_str(), rhs.c_str()) < 0;
  }
};

template <typename T>
T *FindID(const std::vector<T> &vec, const std::string &name) {
  for (const T *t : vec)
    if (ci_comp()(name, t->ID()))
      return &*t;

  std::cout << "[mss]: Error. ID: " << name << " not found." << std::endl;
  exit(EXIT_FAILURE);
}

} // namespace mss

#endif
