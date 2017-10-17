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

  // // TransMat
  // VectorXcd coeff_com = c.TransMat() * inSH.EffectBv(c.Node());
  // EXPECT_TRUE(ApproxVectRv(coeff, coeff_com, 1e-3, 15, true));
}

TEST_F(PeriodicTest, DtN_Map) {
  c.Solve({&inSH}, DFT);

  /// At the nodes:
  // Reference:
  VectorXcd ref_w(c.NumBv() / 2), ref_t(c.NumBv() / 2);
  for (size_t i = 0; i < c.NumNode(); i++) {
    auto tmp = c.Resultant(c.Node(i), {&inSH}).Bv();
    ref_w(i) = tmp(0);
    ref_t(i) = tmp(1);
  }

  // Computed:
  VectorXcd in     = inSH.EffectBv(c.Node());
  VectorXcd com_ww = c.z1_mat() * in, com_tt = c.z2_mat() * in;
  VectorXcd com_w(c.NumNode()), com_t(c.NumNode());
  for (size_t i = 0; i < c.NumNode(); i++) {
    com_w(i) = com_ww(2 * i);
    com_t(i) = com_tt(2 * i);
  }

  EXPECT_TRUE(ApproxVectRv(ref_w, com_w, 1e-6));
  EXPECT_TRUE(ApproxVectRv(ref_t, com_t, 1e-6));

  /// At the doubled nodes (include the complimentary nodes):
  // Reference:
  VectorXcd ref_w2(c.NumBv()), ref_t2(c.NumBv());
  for (size_t i = 0; i < c.NumBv(); i++) {
    auto tmp  = c.Resultant(c.Boundary().DNode()[i], {&inSH}).Bv();
    ref_w2(i) = tmp(0);
    ref_t2(i) = tmp(1);
  }

  // Computed:
  VectorXcd com_w2 = c.z1_mat() * in, com_t2 = c.z2_mat() * in;
  EXPECT_TRUE(ApproxVectRv(ref_w, com_w, 1e-6));
  EXPECT_TRUE(ApproxVectRv(ref_t, com_t, 1e-6));

  // // Inverse:
  // MatrixXcd I(c.NumBv(), c.NumBv());
  // I.setIdentity();
  // std::cout << (c.z1_mat() * c.z1_mat().inverse() - I).norm() << std::endl;
  // std::cout << (c.z2_mat() * c.z2_mat().inverse() - I).norm() << std::endl;

  // std::cout
  //     << (c.z1_mat() * c.z1_mat().inverse() - I).array().abs().maxCoeff()
  //     << std::endl;

  // MatrixXcd z1_inv = c.z1_mat().fullPivLu().inverse();
  // MatrixXcd DtN    = c.z2_mat() * z1_inv;
  // VectorXcd com_t3 = DtN * com_w2;

  VectorXcd com_t3 = c.z1_mat().fullPivLu().solve(c.z2_mat() * com_w2);

  ApproxVectRv(com_t2, com_t3, 1e-2, 0, true);
}

TEST_F(PeriodicTest, DISABLED_InToRstMap) {
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
