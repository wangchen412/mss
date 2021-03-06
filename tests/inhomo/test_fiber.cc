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

#include "../test.h"

namespace mss {

namespace test {

class FiberTest : public Test {
 protected:
  FiberTest() : Test(__FILE__, "fiber") {}

  const double omega = 1.25664e6;
  Material rubber = {1300, 1.41908e9, 0.832e9};
  Material lead = {11400, 36.32496e9, 8.43e9};
  Matrix matrix = {rubber, omega};

  input::Solution s2{path("input.txt")};
  Matrix matrix2{s2};

  FiberConfig<IP> c1 = {"c1", 30, 300, 1e-3, lead, &matrix};
  FiberConfig<AP> c2 = {"c2", 30, 300, 1e-3, lead, &matrix};
  FiberConfig<AP> c3 = {s2.fiber_config()[0], &matrix2};
  IncidentPlaneSH inSH{matrix2, s2.incident()[0]};

  Fiber<IP> f1 = {&c1};
  Fiber<AP> f2 = {&c2, {1, 2}};
  Fiber<AP> f3 = {&c3};
};

TEST_F(FiberTest, ConfigCtor) {
  // EXPECT_EQ(c1.ColloMat().rows(), 1200);
  // EXPECT_EQ(c1.ColloMat().cols(), 122);
  EXPECT_EQ(c1.NumNode(), 300);
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

  EXPECT_EQ(c2.ColloMat().rows(), 600);
  EXPECT_EQ(c2.ColloMat().cols(), 61);
  EXPECT_EQ(c2.NumNode(), 300);
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
  EXPECT_EQ(f1.NumNode(), 300);
  EXPECT_EQ(f1.NumCoeff(), 122);
  EXPECT_EQ(f1.NumBv(), 1200);
  EXPECT_EQ(f1.Node().size(), 300);

  EXPECT_EQ(f2.Position(), PosiVect(1, 2));
  EXPECT_EQ(f2.Basis(), nullptr);
  EXPECT_EQ(f2.NumNode(), 300);
  EXPECT_EQ(f2.NumCoeff(), 61);
  EXPECT_EQ(f2.NumBv(), 600);
  EXPECT_EQ(f1.Node().size(), 300);
}
TEST_F(FiberTest, Contains) {
  const double r = 1e-3;
  PosiVect s(1, 2);
  CS cs0;
  CS cs10(r + epsilon, 0), cs20(r - epsilon, 0);
  CS cs01(0, -r - epsilon), cs02(0, -r + epsilon);

  EXPECT_TRUE(f1.Contains(&cs0));
  EXPECT_FALSE(f1.Contains(&cs10));
  EXPECT_TRUE(f1.Contains(&cs20));
  EXPECT_FALSE(f1.Contains(&cs01));
  EXPECT_TRUE(f1.Contains(&cs02));

  EXPECT_TRUE(f2.Contains(&(cs0 += s)));
  EXPECT_FALSE(f2.Contains(&(cs10 += s)));
  EXPECT_TRUE(f2.Contains(&(cs20 += s)));
  EXPECT_FALSE(f2.Contains(&(cs01 += s)));
  EXPECT_TRUE(f2.Contains(&(cs02 += s)));
}
TEST_F(FiberTest, DSolve) {
  VectorXcd ref(61);
  ReadCoeff("Coeff_SH1.dat", ref);
  EXPECT_TRUE(ApproxVectRv(ref, f3.Solve({&inSH}, DFT), 1e-3, 10));
  EXPECT_TRUE(ApproxVectRv(ref, c3.Solve({&inSH}, DFT), 1e-3, 10));
}
TEST_F(FiberTest, TT) {
  // Compare with the results computed by previous version of mss.
  // The acceptable relative error is set as 1e-9.

  VectorXcd ref(122);
  ReadCoeff("Coeff_SH1.dat", ref);
  for (int i = 0; i < 61; i++) ref(i) *= f3.Config()->TT(i - 30);
  EXPECT_TRUE(ApproxVectRv(VectorXcd(ref.segment(0, 61)),
                           VectorXcd(ref.segment(61, 61)), 1e-9));
}
TEST_F(FiberTest, PseudoInciCoeff) {
  f2.SetCoeff(f2.CSolve({&inSH}));

  VectorXcd inner(600), outer(600);
  for (size_t i = 0; i < f2.NumNode(); i++) {
    inner.segment<2>(2 * i) = f2.Inner(f2.Node(i)).Bv();
    outer.segment<2>(2 * i) =
        f2.Scatter(f2.Node(i)).Bv() + f2.Pseudo(f2.Node(i)).Bv();
  }
  EXPECT_TRUE(ApproxVectRv(inner, outer, 1e-13));
}
TEST_F(FiberTest, RMatrix) {
  VectorXcd ref(122);
  ReadCoeff("Coeff_SH1.dat", ref);
  VectorXcd com = f3.Config()->RefraMat() * inSH.EffectBv(f3.Node());
  EXPECT_TRUE(ApproxVectRv(VectorXcd(ref.segment(61, 61)), com, 1e-3, 10));
}
TEST_F(FiberTest, CSolve) {
  // The acceptable relative error is set as 1e-4.

  VectorXcd ref(61);
  ReadCoeff("Coeff_SH1.dat", ref);
  EXPECT_TRUE(ApproxVectRv(ref, f3.CSolve({&inSH}), 1e-4));
}
TEST_F(FiberTest, Modes) {
  for (int n = -30; n <= 30; n++) {
    // Specific point:
    double x = 23e-3, y = 12e-3;
    CS p(x, y);
    // The center of the fiber f3 is at the origin. So the r and theta:
    double r = sqrt(x * x + y * y);
    double t = atan(y / x);
    // The wave number of the matrix:
    double km = matrix2.KT();
    // State:
    dcomp w = Hn(n, km * r) * exp(ii * double(n) * t);
    dcomp tzr = matrix2.Material().Mu() * km * 0.5 *
                (Hn(n - 1, km * r) - Hn(n + 1, km * r)) *
                exp(ii * double(n) * t);
    dcomp tzt = matrix2.Material().Mu() * ii * double(n) * w / r;
    Vector<dcomp> tt(tzr, tzt);
    // Rotate CS:
    tt.RotateInPlace(-t);
    // Normalizer:
    dcomp norm = Hn(n, km * 10e-3);
    // Result:
    StateAP st(w / norm, tt / norm);
    EXPECT_EQ(f3.ScatterMode(&p, 30 + n), st);
  }
}
TEST_F(FiberTest, DScatter) {
  std::vector<AP> ref, com;
  ReadSample("Line_sc_SH1.dat", ref);
  EXPECT_EQ(ref.size(), 100);

  f3.SetCoeff(f3.DSolve({&inSH}));
  for (auto& i : SamplePts()) com.emplace_back(f3.Scatter(i));
  EXPECT_EQ(com.size(), 100);

  for (size_t i = 0; i < 100; i++) EXPECT_TRUE(ref[i].isApprox(com[i], 1e-3));
}
TEST_F(FiberTest, CScatter) {
  std::vector<AP> ref, com;
  ReadSample("Line_sc_SH1.dat", ref);
  EXPECT_EQ(ref.size(), 100);

  f3.SetCoeff(f3.CSolve({&inSH}));
  for (auto& i : SamplePts()) com.emplace_back(f3.Scatter(i));
  EXPECT_EQ(com.size(), 100);

  for (size_t i = 0; i < 100; i++) EXPECT_TRUE(ref[i].isApprox(com[i], 1e-5));
}
TEST_F(FiberTest, ScatterBvMat) {
  std::vector<AP> ref_states;
  ReadSample("Line_sc_SH1.dat", ref_states);
  EXPECT_EQ(ref_states.size(), 100);

  VectorXcd ref(200);
  for (size_t i = 0; i < 100; i++) ref.segment<2>(i * 2) = ref_states[i].Bv();
  VectorXcd com = f3.ScatterBvMat(*sp_[0]) * f3.DSolve({&inSH});
  EXPECT_TRUE(ApproxVectRv(ref, com, 1e-3));
}
TEST_F(FiberTest, PsInBvMat) {
  f2.SetCoeff(f2.CSolve({&inSH}));

  VectorXcd ref(600);
  for (size_t i = 0; i < f2.NumNode(); i++)
    ref.segment<2>(2 * i) =
        f2.Inner(f2.Node(i)).Bv() - f2.Scatter(f2.Node(i)).Bv();
  VectorXcd com = f2.PsInBvMat(f2.Node()) * f2.PsInCoeff();
  EXPECT_TRUE(ApproxVectRv(ref, com, 1e-12));
}
TEST_F(FiberTest, PsInBvMatT) {
  f2.SetCoeff(f2.CSolve({&inSH}));

  VectorXcd ref(600);
  for (size_t i = 0; i < f2.NumNode(); i++)
    ref.segment<2>(2 * i) =
        f2.Inner(f2.Node(i)).Bv() - f2.Scatter(f2.Node(i)).Bv();
  VectorXcd com = f2.PsInBvMatT(f2.Node()) * f2.ScatterCoeff();
  EXPECT_TRUE(ApproxVectRv(ref, com, 1e-12));
}

}  // namespace test

}  // namespace mss
