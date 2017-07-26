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

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include "../../src/incident/IncidentInput.h"
#include "../../src/inhomo/ConfigFiber.h"
#include "../../src/inhomo/Fiber.h"

namespace mss {

namespace test {

class FiberTest : public testing::Test {
 protected:
  ~FiberTest() {
    for (auto& i : sampleCS_) delete i;
  }
  CSCPtrs sampleCS_;

  const double omega = 1.25664e6;
  Material rubber    = {1300, 1.41908e9, 0.832e9};
  Material lead      = {11400, 36.32496e9, 8.43e9};
  Matrix matrix      = {rubber, omega};

  input::Solution s2{testDataPath(__FILE__) + "fiber/input.txt"};
  Matrix matrix2{s2};

  ConfigFiber<StateIP> c1 = {"c1", 30, 300, 1e-3, lead, &matrix};
  ConfigFiber<StateAP> c2 = {"c2", 30, 300, 1e-3, lead, &matrix};
  ConfigFiber<StateAP> c3 = {s2.configFiber()[0], &matrix2};
  IncidentPlaneSH inSH{matrix2, s2.incident()[0]};

  Fiber<StateIP> f1 = {&c1};
  Fiber<StateAP> f2 = {&c2, {1, 2}};
  Fiber<StateAP> f3 = {&c3};
};

void FiberTest_ReadCoeff(const std::string& fn, Eigen::VectorXcd& sc,
                         Eigen::VectorXcd& in) {
  std::ifstream file(testDataPath(__FILE__) + fn);
  for (int i = 0; i < sc.size(); i++) file >> sc(i);
  for (int i = 0; i < in.size(); i++) file >> in(i);
  file.close();
}

template <typename T>
void FiberTest_ReadSample(const Fiber<T>& f, const std::string& fn,
                          std::vector<T>& ref, std::vector<T>& com,
                          CSCPtrs& sampleCS) {
  std::ifstream file(testDataPath(__FILE__) + fn);
  std::string ts;
  while (std::getline(file, ts)) {
    std::stringstream tss(ts);
    PosiVect r;
    T st;
    tss >> r >> st;
    sampleCS.push_back(new CS(r));
    ref.emplace_back(st);
    com.emplace_back(f.Scatter(sampleCS.back()));
  }
}

TEST_F(FiberTest, ConfigCtor) {
  EXPECT_EQ(c1.TransMatrix().rows(), 1200);
  EXPECT_EQ(c1.TransMatrix().cols(), 122);
  EXPECT_EQ(c1.NoN(), 300);
  EXPECT_EQ(c1.TopOrder(), 30);
  EXPECT_EQ(c1.ID(), "c1");
  EXPECT_EQ(c1.Radius(), 1e-3);
  EXPECT_EQ(c1.Material(), lead);
  EXPECT_EQ(c1.KL(), omega / lead.CL());
  EXPECT_EQ(c1.KT(), omega / lead.CT());
  EXPECT_EQ(c1.Material_m(), rubber);
  EXPECT_EQ(c1.KL_m(), omega / rubber.CL());
  EXPECT_EQ(c1.KT_m(), omega / rubber.CT());
  EXPECT_EQ(c1.Node().size(), 300);

  EXPECT_EQ(c2.TransMatrix().rows(), 600);
  EXPECT_EQ(c2.TransMatrix().cols(), 61);
  EXPECT_EQ(c2.NoN(), 300);
  EXPECT_EQ(c2.TopOrder(), 30);
  EXPECT_EQ(c2.ID(), "c2");
  EXPECT_EQ(c2.Radius(), 1e-3);
  EXPECT_EQ(c2.Material(), lead);
  EXPECT_EQ(c2.KL(), omega / lead.CL());
  EXPECT_EQ(c2.KT(), omega / lead.CT());
  EXPECT_EQ(c2.Material_m(), rubber);
  EXPECT_EQ(c2.KL_m(), omega / rubber.CL());
  EXPECT_EQ(c2.KT_m(), omega / rubber.CT());
  EXPECT_EQ(c2.Node().size(), 300);
}
TEST_F(FiberTest, FiberCtor) {
  EXPECT_EQ(f1.Position(), PosiVect(0, 0));
  EXPECT_EQ(f1.Basis(), nullptr);
  EXPECT_EQ(f1.NoN(), 300);
  EXPECT_EQ(f1.NoC(), 122);
  EXPECT_EQ(f1.NoE(), 1200);
  EXPECT_EQ(f1.Node().size(), 300);

  EXPECT_EQ(f2.Position(), PosiVect(1, 2));
  EXPECT_EQ(f2.Basis(), nullptr);
  EXPECT_EQ(f2.NoN(), 300);
  EXPECT_EQ(f2.NoC(), 61);
  EXPECT_EQ(f2.NoE(), 600);
  EXPECT_EQ(f1.Node().size(), 300);
}
TEST_F(FiberTest, Contain) {
  const double r = 1e-3;
  PosiVect s(1, 2);
  CS cs0;
  CS cs10(r + epsilon, 0), cs20(r - epsilon, 0);
  CS cs01(0, -r - epsilon), cs02(0, -r + epsilon);

  EXPECT_TRUE(f1.Contain(&cs0));
  EXPECT_FALSE(f1.Contain(&cs10));
  EXPECT_TRUE(f1.Contain(&cs20));
  EXPECT_FALSE(f1.Contain(&cs01));
  EXPECT_TRUE(f1.Contain(&cs02));

  EXPECT_TRUE(f2.Contain(&(cs0 += s)));
  EXPECT_FALSE(f2.Contain(&(cs10 += s)));
  EXPECT_TRUE(f2.Contain(&(cs20 += s)));
  EXPECT_FALSE(f2.Contain(&(cs01 += s)));
  EXPECT_TRUE(f2.Contain(&(cs02 += s)));
}
TEST_F(FiberTest, TT) {
  // Compare with the results computed by previous version of mss.
  // The acceptable relative error is set as 1e-9.

  Eigen::VectorXcd sc(61), in(61);
  FiberTest_ReadCoeff("fiber/Coeff_SH1.dat", sc, in);
  for (int i = 0; i < 61; i++) sc(i) *= f3.Config()->TT(i - 30);
  EXPECT_TRUE(ApproxVectRV(sc, in, 1e-9));
}
TEST_F(FiberTest, Solve) {
  // The acceptable relative error is set as 1e-4.

  Eigen::VectorXcd sc(61), in(61);
  FiberTest_ReadCoeff("fiber/Coeff_SH1.dat", sc, in);
  EXPECT_TRUE(ApproxVectRV(sc, f3.Solve({&inSH}), 1e-4));
}
TEST_F(FiberTest, DISABLED_Scatter) {
  f3.SetCoeff(f3.Solve({&inSH}));
  std::vector<StateAP> ref, com;
  FiberTest_ReadSample(f3, "fiber/Line_sc_SH1.dat", ref, com, sampleCS_);

  EXPECT_EQ(ref.size(), 1000);
  EXPECT_THAT(ref, testing::Eq(com));
}

}  // namespace test

}  // namespace mss
