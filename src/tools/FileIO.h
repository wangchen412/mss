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

#include <algorithm>
#include <string>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace mss {

// Case insensitive comparison of two std::strings.
inline bool iequals(const std::string& a, const std::string& b) {
  size_t sz = a.size();
  if (b.size() != sz) return false;
  for (size_t i = 0; i < sz; i++)
    if (tolower(a[i]) != tolower(b[i])) return false;
  return true;
}

// Case insensitive comparison functor for std::map.
struct ci_comp {
  bool operator()(std::string lhs, std::string rhs) const {
    std::transform(lhs.begin(), lhs.end(), lhs.begin(), ::tolower);
    std::transform(rhs.begin(), rhs.end(), rhs.begin(), ::tolower);
    return lhs.compare(rhs) < 0;
  }
};

inline std::string mss_msg(std::initializer_list<std::string> msg) {
  std::string rst("[mss]: ");
  for (auto& i : msg) rst += i;
  return rst;
}

inline std::string error_msg(std::initializer_list<std::string> msg) {
  std::string rst(mss_msg({"Error. "}));
  for (auto& i : msg) rst += i;
  return rst;
}

inline void print_msg(std::initializer_list<std::string> msg) {
  std::cout << mss_msg(msg) << std::endl;
}

inline void print_error_msg(std::initializer_list<std::string> msg) {
  std::cout << error_msg(msg) << std::endl;
}

inline void exit_error_msg(std::initializer_list<std::string> msg) {
  print_error_msg(msg);
  exit(EXIT_FAILURE);
}

// Find the element of which the ID matches "name" in the std::vector "vec".
template <typename T>
const T* FindID(const std::vector<T>& vec, const std::string& name) {
  for (const auto& t : vec)
    if (iequals(name, t.ID)) return &t;
  print_error_msg({"ID: ", name, " not found."});
  exit(EXIT_FAILURE);
}

template <typename T>
const T* FindCPtrID(const std::vector<const T*>& vec,
                    const std::string& name) {
  for (const auto& t : vec)
    if (iequals(name, t->ID())) return t;
  return nullptr;
}

template <typename T>
T* FindPtrID(const std::vector<T*>& vec, const std::string& name) {
  for (const auto& t : vec)
    if (iequals(name, t->ID())) return t;
  return nullptr;
}

// Generate file name that follows the existing files with same format.
inline std::string NewFileName(const std::string& fileName,
                               const std::string& extension) {
  std::string name(fileName);
  int i = 1;
  while (std::ifstream(name + extension))
    name = fileName + std::to_string((long long)i++);
  return name + extension;
}

// Find the last created file with the same format.
inline std::string FindCurrentFile(const std::string& fileName,
                                   const std::string& extension) {
  std::string name(fileName);
  int i = 1;
  while (std::ifstream(name + extension))
    name = fileName + std::to_string((long long)i++);
  return fileName + std::to_string((long long)--i) + extension;
}

// Determine if the string is "white".
inline bool isWhiteSpace(const std::string input) {
  if (input.empty()) return true;
  for (auto it = input.begin(); it != input.end(); it++)
    if (!isspace(*it)) return false;
  return true;
}

// Skip n lines of the input file stream.
inline void skip(std::ifstream& inputFile, int n,
                 std::string* last = nullptr) {
  std::string tmp;
  for (int i = 0; i < n; i++) std::getline(inputFile, tmp);
  if (last) *last = tmp;
}

// Skip lines of the input file stream until find the objective string, then,
// optionally, skip extra n lines.
inline bool skipUntil(std::ifstream& inputFile, const std::string& obj,
                      int n = 0, std::string* last = nullptr) {
  std::string tmp;
  while (std::getline(inputFile, tmp))
    if (iequals(tmp, obj)) {
      skip(inputFile, n, last);
      return true;
    }
  return false;
}

// Copy file.
inline bool copyFile(const std::string& inFileName, std::ofstream& outFile,
                     const std::string& endStr) {
  std::ifstream inFile(inFileName);
  std::string tmp;
  bool notblank = true;
  while (std::getline(inFile, tmp)) {
    notblank = true;
    if (tmp == endStr) return notblank;
    if (isWhiteSpace(tmp)) notblank = false;
    outFile << tmp << std::endl;
  }
  return notblank;
}

// Convert double to std::string.
inline std::string d2string(double x) {
  std::stringstream s;
  s << x;
  return s.str();
}

class separator {
 public:
  separator(const std::string& s, size_t n = 75) : s_(s), n_(n) {}
  friend std::ostream& operator<<(std::ostream& os, const separator& sep) {
    for (size_t i = 0; i < sep.n_; i++) os << sep.s_;
    return os << std::endl;
  }

 private:
  std::string s_;
  size_t n_;
};

}  // namespace mss

#endif
