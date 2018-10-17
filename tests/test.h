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
  void ReadSample(const std::string& fn, std::vector<T>& ref, int s = 0);
  template <typename T>
  void ReadSample(std::ifstream& file, std::vector<T>& ref, int s = 0);

  // Read computed coefficients from data file.
  void ReadCoeff(const std::string& fn, VectorXcd& ref);

  std::string path(const std::string& fn) const { return path_ + fn; }
  const CSCPtrs& SamplePts(int i = 0) const { return *sp_[i]; }
  const CSCPtrs& SamplePtsBack() const { return *sp_.back(); }

  template <typename T>
  MatrixXcd Extract(const std::vector<T>& states) const;

  std::string path_;          // The path to the data directory.
  std::vector<CSCPtrs*> sp_;  // The sample points.

 public:
  static std::string DataPath(const std::string& path,
                              const std::string& fn = "") {
    return path.substr(0, path.find_last_of("/\\")) + std::string("/data/") +
           fn;
  }
};

template <typename T>
void Test::ReadSample(const std::string& fn, std::vector<T>& ref, int s) {
  std::ifstream file(path(fn));
  ReadSample(file, ref, s);
}

template <typename T>
void Test::ReadSample(std::ifstream& file, std::vector<T>& ref, int s) {
  // Read sample points from a file, then:
  // 1. Create a new empty vector at the end of sp_ (sp_ is a 2D vector of
  //    sample points), and put the points read from the file into it.
  // 2. Put the states read from the file into the vector "ref".

  skip(file, s);
  sp_.push_back(new CSCPtrs);
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

inline void Test::ReadCoeff(const std::string& fn, VectorXcd& ref) {
  std::ifstream file(path(fn));
  for (int i = 0; i < ref.size(); i++) file >> ref(i);
  file.close();
}

template <>
MatrixXcd Test::Extract(const std::vector<StateIP>& s) const {
  MatrixXcd rst(s.size(), 6);
  for (size_t i = 0; i < s.size(); i++) {
    rst(i, 0) = s[i].AngleGLB();
    rst(i, 1) = s[i].Displacement().x;
    rst(i, 2) = s[i].Displacement().y;
    rst(i, 3) = s[i].Stress().xx;
    rst(i, 4) = s[i].Stress().yy;
    rst(i, 5) = s[i].Stress().xy;
  }
  return rst;
}
template <>
MatrixXcd Test::Extract(const std::vector<StateAP>& s) const {
  MatrixXcd rst(s.size(), 4);
  for (size_t i = 0; i < s.size(); i++) {
    rst(i, 0) = s[i].AngleGLB();
    rst(i, 1) = s[i].Displacement().x;
    rst(i, 2) = s[i].Stress().x;
    rst(i, 3) = s[i].Stress().y;
  }
  return rst;
}

}  // namespace test

}  // namespace mss

#endif
