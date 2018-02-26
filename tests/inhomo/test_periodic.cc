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

// Pseudo-incident wave field expansion.
class PeriodicTest : public Test {
 protected:
  Material rubber{1300, 1.41908e9, 0.832e9}, lead{11400, 36.32496e9, 8.43e9};
  Matrix m{rubber, 8e5};
  IncidentPlaneSH inc{m, pi_2 / 2, 1e-5};
  double R{3e-3};
  FiberConfig<AP> fc{"1", 35, 1000, R, lead, &m};
  Fiber<AP> f{&fc};
  double A1{6 * R}, A2{2 * R};
  Boundary<AP> b1{10 * m.KT(), {{-A1, A1}, {A1, -A1}}, &m};
  Boundary<AP> b2{10 * m.KT(), {{-A2, A2}, {A2, -A2}}, &m};
};

// Using Bessel J expansion.
TEST_F(PeriodicTest, DtN_single_Cylindrical) {
  // z is the transformation from scattering wave expansion coefficients to
  // the boundary values of resultant wave field. The incident wave is
  // represented by the Bessel J expansion.
  MatrixXcd z(f.ScatterBvMat(b1.Node()) + f.PsiBvMatT(b1.Node()));
  MatrixXcd z1(z.rows() / 2, z.cols());
  MatrixXcd z2(z.rows() / 2, z.cols());

  for (long i = 0; i < z.rows() / 2; i++) {
    z1.row(i) = z.row(i * 2);
    z2.row(i) = z.row(i * 2 + 1);
  }

  MatrixXcd zz(z1.transpose() * z1);
  MatrixXcd dtn(z2 * zz.inverse() * z1.transpose());

  f.SetCoeff(f.CSolve({&inc}));
  VectorXcd u(b1.NumNode()), t(b1.NumNode());
  for (size_t i = 0; i < b1.NumNode(); i++) {
    Vector2cd tmp = f.Scatter(b1.Node(i)).Bv() + inc.Effect(b1.Node(i)).Bv();
    u(i) = tmp(0);
    t(i) = tmp(1);
  }
  VectorXcd tt = dtn * u;
  EXPECT_TRUE(ApproxVectRv(t, tt, 1e-3, 0, true));
}

// Using plane wave expansion.
TEST_F(PeriodicTest, ColloDMat) {
  VectorXcd ref(inc.EffectDv(f.Node()));
  VectorXcd com(f.ColloDMat() * f.CSolve(inc.EffectBv(f.Node())));
  EXPECT_TRUE(ApproxVectRv(ref, com, 1e-10, 0, true));
}
TEST_F(PeriodicTest, ExPoDBMat) {
  VectorXcd ref(inc.EffectBv(b2.Node()));
  VectorXcd com(b2.ExPoDBMat(f.Node()) * inc.EffectDv(f.Node()));
  EXPECT_TRUE(ApproxVectRv(ref, com, 1e-4, 0, true));
}
TEST_F(PeriodicTest, DtN_single_Plane) {
  MatrixXcd z(f.ScatterBvMat(b2.Node()) +
              b2.ExPoDBMat(f.Node()) * f.ColloDMat());
  MatrixXcd z1(z.rows() / 2, z.cols());
  MatrixXcd z2(z.rows() / 2, z.cols());

  for (long i = 0; i < z.rows() / 2; i++) {
    z1.row(i) = z.row(i * 2);
    z2.row(i) = z.row(i * 2 + 1);
  }

  MatrixXcd dtn(z2 * PseudoInverse(z1));

  f.SetCoeff(f.CSolve({&inc}));
  VectorXcd u(b2.NumNode()), t(b2.NumNode());
  for (size_t i = 0; i < b2.NumNode(); i++) {
    Vector2cd tmp = f.Scatter(b2.Node(i)).Bv() + inc.Effect(b2.Node(i)).Bv();
    u(i) = tmp(0);
    t(i) = tmp(1);
  }
  VectorXcd tt = dtn * u;
  EXPECT_TRUE(ApproxVectRv(t, tt, 1e-4, 0, true));
}

// TEST_F(PeriodicTest, DISABLED_PBC) {
//   MatrixXcd z(f.ScatterBvMat(b.Node()) + f.PsiBvMatT(b.Node()));
//   MatrixXcd z1(z.rows() / 2, z.cols());
//   MatrixXcd z2(z.rows() / 2, z.cols());

//   size_t N = b.NumNode() / 4;
//   z1 = z.block(0, 0, 4 * N, z.cols());
//   for (size_t i = 0; i < N; i++) {
//     z2.block(2 * i, 0, 2, z.cols()) =
//         z.block(2 * (3 * N - 1 - i), 0, 2, z.cols());
//     z2.block(2 * (i + N), 0, 2, z.cols()) =
//         z.block(2 * (4 * N - 1 - i), 0, 2, z.cols());
//   }

//   MatrixXcd zz(z2 - z1), uz(z1);
//   uz.block(uz.rows() / 2, 0, uz.rows() / 2, uz.cols()).setZero();

//   // std::cout << zz << std::endl;

//   MatrixXcd Z(zz.transpose() * zz);
//   MatrixXcd UZ(uz.transpose() * zz + zz.transpose() * uz);
//   MatrixXcd U(uz.transpose() * uz);

//   // std::cout << Z << std::endl;

//   dcomp eta_x = 1;
//   MatrixXcd ZZZ(Z - UZ * eta_x + U * eta_x * eta_x);
//   // std::cout << ZZZ.determinant() << std::endl;
//   // std::cout << ZZZ << std::endl;
// }

}  // namespace test

}  // namespace mss
