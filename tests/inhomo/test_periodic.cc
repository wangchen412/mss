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

TEST_F(PeriodicTest, DISABLED_InvBdMat) {
  VectorXcd ref_in = inSH.EffectBv(c.inhomo(0)->Node());

  VectorXcd in0    = inSH.EffectBv(c.Boundary().Node());
  VectorXcd com_in = c.BdIntMatT() * in0;

  ApproxVectRv(ref_in, com_in, 1e-4, 0, true);

  std::cout << c.BdIntMatT().rows() << "  " << c.BdIntMatT().cols()
            << std::endl;

  MatrixXcd bd_inv  = PseudoInverse(c.BdIntMatT());
  VectorXcd com_in0 = bd_inv * com_in;

  ApproxVectRv(in0, com_in0, 1e-4, 0, true);
}

TEST_F(PeriodicTest, Matrices) {
  c.Solve({&inSH}, DFT);

  // Reference:
  VectorXcd ref(c.Boundary().NumDBv());
  for (size_t i = 0; i < c.Boundary().NumDNode(); i++)
    ref.segment(i * 2, 2) =
        c.Resultant(c.Boundary().DNode()[i], {&inSH}).Bv();

  VectorXcd in_d = inSH.EffectBv(c.Boundary().DNode());
  VectorXcd in   = inSH.EffectBv(c.Boundary().Node());

  VectorXcd coeff(c.NumCoeff());
  size_t n = 0;
  for (size_t i = 0; i < c.inhomo().size(); i++) {
    coeff.segment(n, c.NumCoeff(i)) = c.inhomo(i)->ScatterCoeff();
    n += c.NumCoeff(i);
  }

  // BoundaryModeMat:
  // Computed by using BoundaryModeMat times scattering coefficients.
  VectorXcd com1 = c.BoundaryModeMat() * coeff + in_d;
  EXPECT_TRUE(ApproxVectRv(ref, com1));

  // BoundaryModeMat * TransMat:
  VectorXcd com2 = c.BoundaryModeMat() * c.TransMat() * in + in_d;
  EXPECT_TRUE(ApproxVectRv(ref, com2, 1e-6));

  // BoundaryModeMat * TransMat + Interpolation:
  VectorXcd com3 = c.InToRstMat() * in;
  EXPECT_TRUE(ApproxVectRv(ref, com3, 1e-5));
}

// The SVD of matrix B, the boundary integral matrix.
TEST_F(PeriodicTest, DISABLED_MatrixB) {
  MatrixXcd A = c.BdIntMatT();
  MatrixXcd B(A.rows(), A.cols() / 2);
  MatrixXcd C(A.rows() / 2, A.cols() / 2);

  for (long u = 0; u < A.cols() / 2; u += 2) B.col(u) = A.col(u * 2);
  for (long u = 0; u < B.rows() / 2; u += 2) C.row(u) = B.row(u * 2);


  std::cout << A.rows() << "  " << A.cols() << std::endl;
  std::cout << B.rows() << "  " << B.cols() << std::endl;
  std::cout << C.rows() << "  " << C.cols() << std::endl;


  Eigen::VectorXd svda = A.bdcSvd().singularValues();
  Eigen::VectorXd svdb = B.bdcSvd().singularValues();
  Eigen::VectorXd svdc = C.bdcSvd().singularValues();

  std::cout << svda << std::endl;
  std::cout << "++++++++++++++++++++++" << std::endl;
  std::cout << svdb << std::endl;
  std::cout << "++++++++++++++++++++++" << std::endl;
  std::cout << svdc << std::endl;
}

TEST_F(PeriodicTest, DISABLED_DtN_Map) {
  c.Solve({&inSH}, DFT);

  VectorXcd in_d = inSH.EffectBv(c.Boundary().DNode());
  VectorXcd in   = inSH.EffectBv(c.Boundary().Node());

  /// At the doubled nodes (include the complimentary nodes):
  // Reference:
  VectorXcd ref_w(c.Boundary().NumDBv() / 2);
  VectorXcd ref_t(c.Boundary().NumDBv() / 2);
  for (size_t i = 0; i < c.Boundary().NumDNode(); i++) {
    auto tmp = c.Resultant(c.Boundary().DNode()[i], {&inSH}).Bv();
    ref_w(i) = tmp(0);
    ref_t(i) = tmp(1);
  }

  // Computed:
  VectorXcd com_w = c.z1_mat() * in, com_t = c.z2_mat() * in;
  EXPECT_TRUE(ApproxVectRv(ref_w, com_w, 1e-5));
  EXPECT_TRUE(ApproxVectRv(ref_t, com_t, 1e-5));

  /// PseudoInverse:
  // MatrixXcd z1_inv = PseudoInverse(c.z1_mat());
  // MatrixXcd DtN = c.z2_mat() * z1_inv;
  // VectorXcd com_t2 = DtN * com_w;
  // ApproxVectRv(com_t, com_t2, 1, 0, true);

  /// Least square:
  auto svd = c.z1_mat().jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);
  VectorXcd com_in = svd.solve(com_w);
  ApproxVectRv(in, com_in, 1, 0, true);
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
