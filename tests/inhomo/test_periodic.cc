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
  PeriodicTest() : Test(__FILE__, "periodic") {}

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
  MatrixXcd z(f.ScatterBvMat(b1.Node()) + f.PsInBvMatT(b1.Node()));
  MatrixXcd z1(z.rows() / 2, z.cols());
  MatrixXcd z2(z.rows() / 2, z.cols());

  for (long i = 0; i < z.rows() / 2; i++) {
    z1.row(i) = z.row(i * 2);
    z2.row(i) = z.row(i * 2 + 1);
  }

  MatrixXcd zz(z1.adjoint() * z1);
  MatrixXcd dtn(z2 * zz.inverse() * z1.adjoint());

  f.SetCoeff(f.CSolve({&inc}));
  VectorXcd u(b1.NumNode()), t(b1.NumNode());
  for (size_t i = 0; i < b1.NumNode(); i++) {
    Vector2cd tmp = f.Scatter(b1.Node(i)).Bv() + inc.Effect(b1.Node(i)).Bv();
    u(i) = tmp(0);
    t(i) = tmp(1);
  }
  VectorXcd tt = dtn * u;
  EXPECT_TRUE(ApproxVectRv(t, tt, 1e-4, 0, true));
}

// Using plane wave expansion.
TEST_F(PeriodicTest, ColloDMat) {
  VectorXcd ref(inc.EffectDv(f.Node()));
  VectorXcd com(f.ColloDMat() * f.CSolve(inc.EffectBv(f.Node())));
  EXPECT_TRUE(ApproxVectRv(ref, com, 1e-10, 0, true));
}
TEST_F(PeriodicTest, PlaneEBMat) {
  VectorXcd ref(inc.EffectBv(b2.Node()));
  VectorXcd com(b2.PlaneEBMat(f.Node()) * inc.EffectDv(f.Node()));
  EXPECT_TRUE(ApproxVectRv(ref, com, 1e-4, 0, true));
}
TEST_F(PeriodicTest, DtN_single_plane) {
  MatrixXcd z(f.ScatterBvMat(b2.Node()) +
              b2.PlaneEBMat(f.Node()) * f.ColloDMat());
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
TEST_F(PeriodicTest, DtN_single_cylindrical) {
  input::Solution input{path("input.txt")};
  Matrix matrix(input);
  IncidentPlaneSH inc(input);
  AssemblyConfig<AP> ac(input.config(), &matrix);

  MatrixXcd z(ac.ResBvMat(ac.Node()));
  // z << ac.CylinEBMat(0), ac.CylinEBMat(1), ac.CylinEBMat(2),
  // ac.CylinEBMat(3); z += ac.inhomo(0)->ScatterBvMat(ac.Node());

  MatrixXcd z1(z.rows() / 2, z.cols()), z2(z.rows() / 2, z.cols());
  for (long i = 0; i < z.rows() / 2; i++) {
    z1.row(i) = z.row(i * 2);
    z2.row(i) = z.row(i * 2 + 1);
  }

  MatrixXcd dtn(z2 * PseudoInverse(z1));
  ac.DSolve({&inc});
  VectorXcd u(ac.NumNode()), t(ac.NumNode());
  for (size_t i = 0; i < ac.NumNode(); i++) {
    Vector2cd tmp = ac.Resultant(ac.Node(i), {&inc}).Bv();
    u(i) = tmp(0);
    t(i) = tmp(1);
  }
  VectorXcd tt = dtn * u;
  EXPECT_TRUE(ApproxVectRv(t, tt, 1e-4, 0, true));
}

// Eigenvalue problem
TEST_F(PeriodicTest, DISABLED_Eigenvalue_DtN_check) {
  input::Solution input{path("input2.txt")};
  Matrix matrix(input);
  IncidentPlaneSH inc(input);
  AssemblyConfig<AP> ac(input.config(), &matrix);

  MatrixXcd z(ac.ResBvMat(ac.Node()));
  // z << ac.CylinEBMat(0), ac.CylinEBMat(1), ac.CylinEBMat(2),
  // ac.CylinEBMat(3); z += ac.inhomo(0)->ScatterBvMat(ac.Node());

  MatrixXcd z1(z.rows() / 2, z.cols()), z2(z.rows() / 2, z.cols());
  for (long i = 0; i < z.rows() / 2; i++) {
    z1.row(i) = z.row(i * 2);
    z2.row(i) = z.row(i * 2 + 1);
  }

  MatrixXcd dtn(z2 * PseudoInverse(z1));
  ac.DSolve({&inc});
  VectorXcd u(ac.NumNode()), t(ac.NumNode());
  for (size_t i = 0; i < ac.NumNode(); i++) {
    Vector2cd tmp = ac.Resultant(ac.Node(i), {&inc}).Bv();
    u(i) = tmp(0);
    t(i) = tmp(1);
  }
  VectorXcd tt = dtn * u;
  EXPECT_TRUE(ApproxVectRv(t, tt, 1e-4, 0, true));
}
TEST_F(PeriodicTest, Eigenvalue_single) {
  input::Solution input{path("input2.txt")};
  std::ifstream data(path("BlochK_45.dat"));
  std::string tmp;
  Eigen::VectorXd omega(14), ref_k(14), com_k(14);
  skipUntil(data, "omega:");
  for (int i = 0; i < 14; i++) {
    getline(data, tmp);
    std::stringstream(tmp) >> omega(i);
  }
  skipUntil(data, "k:");
  for (int i = 0; i < 14; i++) {
    getline(data, tmp);
    std::stringstream(tmp) >> ref_k(i);
  }

  Material nickle(input.material()[0]);
  Material aluminum(input.material()[1]);

  for (int n = 0; n < 14; n++) {
    Matrix matrix(aluminum, omega(n));
    AssemblyConfig<AP> ac(input.config(), &matrix);
    ac.Boundary().ReverseEdge();

    MatrixXcd z1(2 * (ac.Edge(0).size() + ac.Edge(1).size()), ac.NumCoeff());
    MatrixXcd z2(2 * (ac.Edge(2).size() + ac.Edge(3).size()), ac.NumCoeff());
    z1 << ac.ResBvMat(ac.Edge(0)), ac.ResBvMat(ac.Edge(1));
    z2 << ac.ResBvMat(ac.Edge(2)), ac.ResBvMat(ac.Edge(3));
    for (long i = 1; i < z1.rows(); i += 2) z1.row(i) *= -1;
    for (long i = 0; i < z1.rows(); i++) {
      dcomp p = z1.row(i).array().mean();
      // dcomp p = GeometricMean(z1.row(i).array());
      z1.row(i) /= p;
      z2.row(i) /= p;
    }

    MatrixXcd A(PseudoInverse(z1) * z2);
    Eigen::ComplexEigenSolver<MatrixXcd> ces;
    ces.compute(A);
    VectorXcd ev = ces.eigenvalues();
    VectorXcd mv = ev.array().abs();
    auto ue = FindUnitEigenvalue(ev, 0.01);
    if (ue.empty()) exit_error_msg({"No unit eigenvalues."});

    if ((log(ev(ue[0])) / ii).real() > 0)
      com_k(n) = (log(ev(ue[0])) / ii / pi).real();
    else
      com_k(n) = (log(ev(ue[1])) / ii / pi).real();
  }

  EXPECT_TRUE(ApproxVectRv(ref_k, com_k, 2e-2));
}

}  // namespace test

}  // namespace mss
