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

class PeriodicTest : public Test {
 protected:
  PeriodicTest() : Test(__FILE__, "periodic") {}

  input::Solution s{path("input.txt")};
  Matrix matrix{s};
  AssemblyConfig<AP> c{s.config(), &matrix};
  IncidentPlaneSH inSH{matrix, s.incident()[0]};
};

TEST_F(PeriodicTest, BoundaryModeMat) {
  c.Solve({&inSH}, DFT);

  // Reference:
  VectorXcd ref(c.NumBv() * 2);
  for (size_t i = 0; i < c.NumNode() * 2; i++)
    ref.segment(i * 2, 2) =
        c.Resultant(c.Boundary().DNode()[i], {&inSH}).Bv();

  // Computed by using BoundaryModeMat times scattering coefficients.
  VectorXcd coeff(c.NumCoeff());
  size_t n = 0;
  for (size_t i = 0; i < c.inhomo().size(); i++) {
    coeff.segment(n, c.NumCoeff(i)) = c.inhomo(i)->ScatterCoeff();
    n += c.NumCoeff(i);
  }

  VectorXcd com =
      c.BoundaryModeMat() * coeff + inSH.EffectBv(c.Boundary().DNode());

  EXPECT_TRUE(ApproxVectRv(ref, com));


  // TransMat
  VectorXcd coeff_com = c.TransMat() * inSH.EffectBv(c.Node());
  EXPECT_TRUE(ApproxVectRv(coeff, coeff_com, 1e-3, 15, true));
}

TEST_F(PeriodicTest, DtN_Map) {
  c.Solve({&inSH}, DFT);

  // Reference:
  VectorXcd ref_w(c.NumBv() / 2), ref_t(c.NumBv() / 2);
  for (size_t i = 0; i < c.NumNode(); i++) {
    auto tmp = c.Resultant(c.Node(i), {&inSH}).Bv();
    ref_w(i) = tmp(0);
    ref_t(i) = tmp(1);
  }

  VectorXcd in     = inSH.EffectBv(c.Node());
  VectorXcd com_ww = c.z1_mat() * in, com_tt = c.z2_mat() * in;
  VectorXcd com_w(c.NumNode()), com_t(c.NumNode());
  for (size_t i = 0; i < c.NumNode(); i++) {
    com_w(i) = com_ww(2 * i);
    com_t(i) = com_tt(2 * i);
  }

  std::cout << "Node number: " << c.NumNode() << std::endl;
  std::cout << "ref_w vector size: " << ref_w.size() << std::endl;
  std::cout << "com_w vector size: " << com_w.size() << std::endl;
  std::cout << "com_ww size: " << com_ww.size() << std::endl;

  std::ofstream file1("com_ww.txt"), file2("com_w.txt"), file3("com_www.txt");
  file1 << com_ww << std::endl;
  file2 << com_w << std::endl;
  file1.close();
  file2.close();

  for (long i = 0; i < com_ww.size(); i++) {
    dcomp tmp = com_ww(i);
    file3 << tmp << std::endl;
  }
  file3.close();

  // std::cout << com_ww(2 * c.NumNode() - 4) << std::endl;
  // std::cout << com_tt(2 * c.NumNode() - 4) << std::endl;

  // std::cout << com_w(c.NumNode() - 2) << std::endl;
  // std::cout << com_t(c.NumNode() - 2) << std::endl;

  // std::cout << com_ww(2 * c.NumNode() - 1) << std::endl;
  // std::cout << com_tt(2 * c.NumNode() - 1) << std::endl;

  // std::cout << com_ww << std::endl;

  // std::cout << c.z1_mat().row(c.NumNode() - 1) << std::endl;
  // std::cout << c.z2_mat().row(c.NumNode() - 1) << std::endl;

  // EXPECT_TRUE(ApproxVectRv(ref_w, com_w, 1e-2, 0, true));
  // EXPECT_TRUE(ApproxVectRv(ref_t, com_t, 1e-2, 0, true));
}

TEST_F(PeriodicTest, InToRstMap) {
  c.Solve({&inSH}, DFT);

  // Reference:
  //   The effects of incident wave at the complimentary nodes are given by
  //   the real incident.
  VectorXcd ref(c.NumBv());
  for (size_t i = 0; i < c.NumNode(); i++)
    ref.segment(i * 2, 2) =
        c.Resultant(c.Boundary().DNode()[i], {&inSH}).Bv();

  // Computed by using BoundaryModeMat times TransMat times incident.
  //   The effects of incident wave at the complimentary nodes are given by
  //   the linear interpolation.
  MatrixXcd m = c.BoundaryModeMat() * c.TransMat();

  VectorXcd com = m * inSH.EffectBv(c.Node());

  EXPECT_TRUE(ApproxVectRv(ref, com, 1e-3, 0, true));
  EXPECT_FALSE(ApproxVectRv(ref, com));
}

TEST_F(PeriodicTest, DISABLED_CharPoly) {
  // std::cout << c.CharPoly(exp(ii * k * c.Width()), 1) << std::endl;
  // std::cout << c.CharPoly(1, 1) << std::endl;

  // std::cout << c.Z_mat(1, 1).size() << std::endl;
  // std::cout << c.CharPoly(1, 1) << std::endl;

  std::ofstream file("exp.dat");

  int n     = 20;
  double ll = 250, hh = 260;
  double dk = (hh - ll) / n;

  for (int i = 0; i < n; i++) {
    double k  = i * dk + ll;
    dcomp psx = exp(ii * k * c.Width());
    file << k << "\t" << DeterExpon(c.Y_mat(psx, 1)) << std::endl;
  }

  file.close();
}

TEST_F(PeriodicTest, DISABLED_InvMat2) {
  // MatrixXcd m = c.Y_mat(exp(ii * pi2), 1);
  // MatrixXcd m = c.z1_mat();  // Singular
  MatrixXcd m = c.BoundaryModeMat() * c.TransMat();

  MatrixXcd I(m.rows(), m.rows());
  I.setIdentity();
  std::cout << (m * m.inverse() - I).norm() << std::endl;
}

TEST_F(PeriodicTest, DISABLED_InvMat) {
  MatrixXcd m = c.Y_mat(1, 1);

  auto lu = m.partialPivLu();

  MatrixXcd P = lu.permutationP();
  MatrixXcd L = lu.matrixLU().triangularView<Eigen::Upper>();
  MatrixXcd U = lu.matrixLU().triangularView<Eigen::UnitLower>();

  MatrixXcd Li = L.inverse();
  MatrixXcd Ui = U.inverse();

  MatrixXcd I(Li.rows(), Li.rows());
  I.setIdentity();

  // std::cout << L*Li << std::endl;
  std::cout << (L * Li - I).norm() << std::endl;
}

}  // namespace test

}  // namespace mss
