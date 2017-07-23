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
  const double omega = 1.25664e6;
  Material rubber = {1300, 1.41908e9, 0.832e9};
  Material lead = {11400, 36.32496e9, 8.43e9};
  Matrix matrix = {rubber, omega};

  ConfigFiber<StateIP> c1 = {"c1", 30, 300, 1e-3, lead, &matrix};
  ConfigFiber<StateAP> c2 = {"c2", 30, 300, 1e-3, lead, &matrix};

  Fiber<StateIP> f1 = {&c1};
  Fiber<StateAP> f2 = {&c2, {1, 2}};
  Fiber<StateAP> f3 = {&c2};
};

void FiberTest_ReadFile(const std::string& fn, Eigen::VectorXcd& sc,
                        Eigen::VectorXcd& in) {
  std::ifstream file(testDataPath(__FILE__) + fn);
  skip(file, 4);
  for (int i = 0; i < sc.size(); i++) file >> sc(i);
  for (int i = 0; i < in.size(); i++) file >> in(i);
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
  // Only -12 ~ +12 order coefficients are tested and the acceptable relative
  // error is set as 1e-6. The reference data is computed by the previous
  // version mss.

  int k = 12;
  const double re = 1e-6;

  Eigen::VectorXcd sc(61), in(61);
  FiberTest_ReadFile("Single_SH1.dat", sc, in);
  for (int i = 30 - k; i < 30 + k; i++)
    EXPECT_TRUE(near(in(i), sc(i) * f3.Config()->TT(i - 30), re));
}
TEST_F(FiberTest, SingleScattering) {
  // -14 ~ +14 order coefficients are tested to be near the results obtained
  // by previous program. The acceptable relative error is 1e-4.

  IncidentPlaneSH inSH(matrix, 1.2, 2.3e-6, 3.4);
  Eigen::VectorXcd vsh = inSH.EffectBV(f3.Node());
  Eigen::MatrixXcd Q3 = f3.ModeMatrix(&f3);
  auto svd3 = Q3.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::VectorXcd coeff3 = svd3.solve(vsh);

  int k = 14;
  const double re = 1e-4;

  Eigen::VectorXcd sc(61), in(61);
  FiberTest_ReadFile("Single_SH1.dat", sc, in);
  for (int i = 30 - k; i < 30 + k; i++)
    EXPECT_TRUE(near(sc(i), coeff3(i), re));
}

}  // namespace test

}  // namespace mss
