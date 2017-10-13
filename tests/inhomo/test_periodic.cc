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

  Eigen::VectorXcd rst1(c.NumBv() * 2);
  for (size_t i = 0; i < c.NumNode() * 2; i++)
    rst1.segment(i * 2, 2) =
        c.Resultant(c.Boundary().DNode()[i], {&inSH}).Bv();

  Eigen::VectorXcd coeff(c.NumCoeff());
  size_t n = 0;
  for (size_t i = 0; i < c.inhomo().size(); i++) {
    coeff.segment(n, c.NumCoeff(i)) = c.inhomo(i)->ScatterCoeff();
    n += c.NumCoeff(i);
  }

  Eigen::VectorXcd rst2 =
      c.BoundaryModeMat() * coeff + inSH.EffectBv(c.Boundary().DNode());

  EXPECT_TRUE(ApproxVectRv(rst1, rst2));
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

TEST_F(PeriodicTest, InvMat2) {
  MatrixXcd m = c.Y_mat(exp(ii * pi2), 1);
  MatrixXcd I(m.rows(), m.rows());
  I.setIdentity();
  std::cout << (m*m.inverse() - I).norm() << std::endl;
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

  //std::cout << L*Li << std::endl;
  std::cout << (L*Li - I).norm() << std::endl;
}

}  // namespace test

}  // namespace mss
