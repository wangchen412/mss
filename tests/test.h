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

#ifndef MSS_TEST_H
#define MSS_TEST_H

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include "../src/core/Solution.h"

namespace mss {

namespace test {

class Test : public testing::Test {
 protected:
  Test(const std::string& path = __FILE__, const std::string& subDir = "") {
    path_ = DataPath(path, subDir);
    if (subDir != "") path_ += "/";
  }

  virtual ~Test() {
    for (auto& i : sp_)
      for (auto& j : *i) delete j;
  }

  // Read sample points and states from data file.
  template <typename T>
  void ReadSample(const std::string& fn, std::vector<T>& ref);

  // Read computed coefficients from data file.
  void ReadCoeff(const std::string& fn, Eigen::VectorXcd& ref);

  std::string path(const std::string& fn) const { return path_ + fn; }
  const CSCPtrs& SamplePts(int i) const { return *sp_[i]; }
  const CSCPtrs& SamplePtsBack() const { return *sp_.back(); }

  std::string path_;          // The path to the data directory.
  std::vector<CSCPtrs*> sp_;  // The sample points.

 public:
  static std::string DataPath(const std::string& path,
                              const std::string& fn = "") {
    return path.substr(0, path.rfind("/")) + std::string("/data/") + fn;
  }
};

template <typename T>
inline void Test::ReadSample(const std::string& fn, std::vector<T>& ref) {
  sp_.push_back(new CSCPtrs);
  std::ifstream file(path(fn));
  std::string ts;
  while (std::getline(file, ts)) {
    std::stringstream tss(ts);
    PosiVect r;
    T st;
    tss >> r >> st;
    sp_.back()->push_back(new CS(r));
    ref.emplace_back(st);
  }
  file.close();
}

inline void Test::ReadCoeff(const std::string& fn, Eigen::VectorXcd& ref) {
  std::ifstream file(path(fn));
  for (int i = 0; i < ref.size(); i++) file >> ref(i);
  file.close();
}

}  // namespace test

}  // namespace mss

#endif
