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
TEST_F(PeriodicTest, DISABLED_DtN_single_Cylindrical) {
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
TEST_F(PeriodicTest, DISABLED_ColloDMat) {
  VectorXcd ref(inc.EffectDv(f.Node()));
  VectorXcd com(f.ColloDMat() * f.CSolve(inc.EffectBv(f.Node())));
  EXPECT_TRUE(ApproxVectRv(ref, com, 1e-10, 0, true));
}
// TEST_F(PeriodicTest, DISABLED_PlaneEBMat) {
//   VectorXcd ref(inc.EffectBv(b2.Node()));
//   VectorXcd com(b2.PlaneEBMat(f.Node()) * inc.EffectDv(f.Node()));
//   EXPECT_TRUE(ApproxVectRv(ref, com, 1e-4, 0, true));
// }
// TEST_F(PeriodicTest, DISABLED_DtN_single_plane) {
//   MatrixXcd z(f.ScatterBvMat(b2.Node()) +
//               b2.PlaneEBMat(f.Node()) * f.ColloDMat());
//   MatrixXcd z1(z.rows() / 2, z.cols());
//   MatrixXcd z2(z.rows() / 2, z.cols());

//   for (long i = 0; i < z.rows() / 2; i++) {
//     z1.row(i) = z.row(i * 2);
//     z2.row(i) = z.row(i * 2 + 1);
//   }

//   MatrixXcd dtn(z2 * PseudoInverse(z1));

//   f.SetCoeff(f.CSolve({&inc}));
//   VectorXcd u(b2.NumNode()), t(b2.NumNode());
//   for (size_t i = 0; i < b2.NumNode(); i++) {
//     Vector2cd tmp = f.Scatter(b2.Node(i)).Bv() + inc.Effect(b2.Node(i)).Bv();
//     u(i) = tmp(0);
//     t(i) = tmp(1);
//   }
//   VectorXcd tt = dtn * u;
//   EXPECT_TRUE(ApproxVectRv(t, tt, 1e-4, 0, true));
// }
TEST_F(PeriodicTest, DISABLED_DtN_single_cylindrical) {
  input::Solution input{path("input.txt")};
  Matrix matrix(input);
  IncidentPlaneSH inc(matrix, input.incident()[0]);
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
  EXPECT_TRUE(ApproxVectRv(t, tt, 1e-3, 0, true));
}

// Eigenvalue problem
TEST_F(PeriodicTest, DISABLED_Eigenvalue_DtN_check) {
  input::Solution input{path("input2.txt")};
  Matrix matrix(input);
  IncidentPlaneSH inc(matrix, input.incident()[0]);
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
TEST_F(PeriodicTest, DISABLED_Eigenvalue_single) {
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
    AssemblyConfig<AP> ac(input.assembly_config()[0], &matrix);
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
TEST_F(PeriodicTest, DISABLED_ResMat_multiple) {
  // A preliminary check for the ResMat.

  input::Solution input{path("input2.txt")};
  Matrix matrix(input);
  IncidentPlaneSH inc(matrix, input.incident()[0]);
  AssemblyConfig<AP> ac(input.assembly_config()[1], &matrix);
  CSCPtrs node = ac.EdgeNode();
  MatrixXcd z(ac.ResBvMat(node));

  ac.DSolve({&inc});
  VectorXcd com = z * ac.ScatterCoeff();
  VectorXcd ref(2 * node.size());
  for (size_t i = 0; i < node.size(); i++)
    ref.segment<2>(2 * i) = ac.Resultant(node[i], {&inc}).Bv();

  EXPECT_TRUE(ApproxVectRv(ref, com, 1e-4, 0, true));
}
TEST_F(PeriodicTest, DISABLED_ResMat_multiple_DtN) {
  // Check if the relation derived can be used in BEM solving.

  input::Solution input{path("input2.txt")};
  Matrix matrix(input);
  IncidentPlaneSH inc(matrix, input.incident()[0]);
  AssemblyConfig<AP> ac(input.assembly_config()[1], &matrix);
  AssemblyConfig<AP> ac_hi(input.assembly_config()[2], &matrix);
  CSCPtrs node = ac.EdgeNode();
  MatrixXcd z(ac.ResBvMat(node));

  MatrixXcd zw(z.rows() / 2, z.cols()), zt(z.rows() / 2, z.cols());
  for (long i = 0; i < zw.rows(); i++)
    zw.row(i) = z.row(2 * i), zt.row(i) = z.row(2 * i + 1);

  ac.DSolve({&inc});
  ac_hi.DSolve({&inc});
  VectorXcd w(node.size()), t(node.size());
  VectorXcd w_hi(w), t_hi(t);
  for (size_t i = 0; i < node.size(); i++) {
    Vector2cd tmp = ac.Resultant(node[i], {&inc}).Bv();
    w(i) = tmp(0);
    t(i) = tmp(1);
    tmp = ac_hi.Resultant(node[i], {&inc}).Bv();
    w_hi(i) = tmp(0);
    t_hi(i) = tmp(1);
  }

  // ApproxVectRv(w, w_hi, 1e-8, 0, true);
  // ApproxVectRv(t, t_hi, 1e-8, 0, true);

  // VectorXcd coeff_r = ac.ScatterCoeff();
  // VectorXcd coeff_w = zw.jacobiSvd(40).solve(w);
  // VectorXcd coeff_t = zt.jacobiSvd(40).solve(t);
  // EXPECT_TRUE(ApproxVectRv(coeff_r, coeff_w, 1e-2, 0, true));
  // EXPECT_TRUE(ApproxVectRv(coeff_r, coeff_w, 1e-2, 0, true));

  MatrixXcd DtN = zt * PseudoInverse(zw);
  VectorXcd tt = DtN * w;
  EXPECT_TRUE(ApproxVectRv(t, tt, 3e-4));

  // long n = DtN.rows() / 2;
  // MatrixXcd z11 = DtN.block(0, 0, n, n);
  // MatrixXcd z12 = DtN.block(0, n, n, n);
  // MatrixXcd z21 = DtN.block(n, 0, n, n);
  // MatrixXcd z22 = DtN.block(n, n, n, n);
  // MatrixXcd iii = MatrixXcd::Identity(n, n);
  // MatrixXcd zzz = MatrixXcd::Zero(n, n);

  // MatrixXcd A(2 * n, 2 * n), B(2 * n, 2 * n);
  // A << z11 - z22, -z21, -iii, zzz;
  // B << -z12, zzz, zzz, -iii;

  // writeMatrix(A, "A");
  // writeMatrix(B, "B");
}
TEST_F(PeriodicTest, DISABLED_ResMat_Larger_DtN) {
  // Check if the DtN map derived can be derived for with larger cells.

  input::Solution input{path("input2.txt")};
  Matrix matrix(input);
  IncidentPlaneSH inc(matrix, input.incident()[0]);
  AssemblyConfig<AP> ac(input.assembly_config()[3], &matrix);
  AssemblyConfig<AP> ac_hi(input.assembly_config()[4], &matrix);
  CSCPtrs node = ac.EdgeNode();
  MatrixXcd z(ac.ResBvMat(node));
  MatrixXcd zw(z.rows() / 2, z.cols()), zt(z.rows() / 2, z.cols());

  for (long i = 0; i < zw.rows(); i++)
    zw.row(i) = z.row(2 * i), zt.row(i) = z.row(2 * i + 1);

  ac.DSolve({&inc});
  ac_hi.DSolve({&inc});
  VectorXcd w(node.size()), t(node.size());
  VectorXcd w_hi(w), t_hi(t);
  for (size_t i = 0; i < node.size(); i++) {
    Vector2cd tmp = ac.Resultant(node[i], {&inc}).Bv();
    w(i) = tmp(0);
    t(i) = tmp(1);
    tmp = ac_hi.Resultant(node[i], {&inc}).Bv();
    w_hi(i) = tmp(0);
    t_hi(i) = tmp(1);
  }

  // // VectorXcd coeff_r = ac.ScatterCoeff();
  // // VectorXcd coeff_w = zw.jacobiSvd(40).solve(w);
  // // VectorXcd coeff_t = zt.jacobiSvd(40).solve(t);
  // // EXPECT_TRUE(ApproxVectRv(coeff_r, coeff_w, 1e-2, 0, true));
  // // EXPECT_TRUE(ApproxVectRv(coeff_r, coeff_w, 1e-2, 0, true));

  EXPECT_TRUE(ApproxVectRv(w, w_hi, 1e-3, 0, true));
  EXPECT_TRUE(ApproxVectRv(t, t_hi, 1e-3, 0, true));

  VectorXcd tt = zt * PseudoInverse(zw) * w;
  // VectorXcd tt = zt * zw.jacobiSvd(40).solve(w);

  EXPECT_TRUE(ApproxVectRv(t, tt, 2e-2));
}
TEST_F(PeriodicTest, DISABLED_Eigenvalue_multiple_DtN) {
  input::Solution input{path("input2.txt")};
  Matrix matrix(input);
  AssemblyConfig<AP> ac(input.assembly_config()[1], &matrix);
  ac.Boundary().ReverseEdge();

  MatrixXcd z1(2 * (ac.Edge(0).size() + ac.Edge(1).size()), ac.NumCoeff());
  MatrixXcd z2(2 * (ac.Edge(2).size() + ac.Edge(3).size()), ac.NumCoeff());
  z1 << ac.ResBvMat(ac.Edge(0)), ac.ResBvMat(ac.Edge(1));
  z2 << ac.ResBvMat(ac.Edge(2)), ac.ResBvMat(ac.Edge(3));
  for (long i = 1; i < z1.rows(); i += 2) z1.row(i) *= -1;
  MatrixXcd z(z1.rows() + z2.rows(), z1.cols());
  z << z1, z2;

  MatrixXcd zw(z.rows() / 2, z.cols()), zt(z.rows() / 2, z.cols());
  for (long i = 0; i < zw.rows(); i++)
    zw.row(i) = z.row(2 * i), zt.row(i) = z.row(2 * i + 1);
  MatrixXcd DtN = zt * PseudoInverse(zw);

  long n = DtN.rows() / 2;
  MatrixXcd z11 = DtN.block(0, 0, n, n), z12 = DtN.block(0, n, n, n);
  MatrixXcd z21 = DtN.block(n, 0, n, n), z22 = DtN.block(n, n, n, n);
  MatrixXcd iii = MatrixXcd::Identity(n, n), zzz = MatrixXcd::Zero(n, n);

  MatrixXcd A(2 * n, 2 * n), B(2 * n, 2 * n);
  A << z11 - z22, -z21, -iii, zzz;
  B << -z12, zzz, zzz, -iii;

  writeMatrix(A, "A");
  writeMatrix(B, "B");
}
TEST_F(PeriodicTest, DISABLED_Eigenvalue_multiple) {
  input::Solution input{path("input2.txt")};
  Matrix matrix(input);
  AssemblyConfig<AP> ac(input.assembly_config()[1], &matrix);
  ac.Boundary().ReverseEdge();

  MatrixXcd z1(2 * (ac.Edge(0).size() + ac.Edge(1).size()), ac.NumCoeff());
  MatrixXcd z2(2 * (ac.Edge(2).size() + ac.Edge(3).size()), ac.NumCoeff());
  z1 << ac.ResBvMat(ac.Edge(0)), ac.ResBvMat(ac.Edge(1));
  z2 << ac.ResBvMat(ac.Edge(2)), ac.ResBvMat(ac.Edge(3));
  for (long i = 1; i < z1.rows(); i += 2) z1.row(i) *= -1;
  MatrixXcd z(z1.rows() + z2.rows(), z1.cols());

  writeMatrix(z1, "z1");
  writeMatrix(z2, "z2");
}

TEST_F(PeriodicTest, ResMat_improve) {
  input::Solution input{path("input2.txt")};
  Matrix matrix(input);
  IncidentPlaneSH inc(matrix, input.incident()[0]);
  AssemblyConfig<AP> ac(input.assembly_config()[0], &matrix);
  ac.Boundary().ReverseEdge();
  CSCPtrs node = ac.EdgeNode();
  MatrixXcd z(ac.ResBvMat(node));
  MatrixXcd zw(z.rows() / 2, z.cols()), zt(z.rows() / 2, z.cols());
  for (long i = 0; i < zw.rows(); i++)
    zw.row(i) = z.row(2 * i), zt.row(i) = z.row(2 * i + 1);

  ac.DSolve({&inc});
  // ac.CSolve({&inc});
  VectorXcd w(node.size()), t(node.size());
  for (size_t i = 0; i < node.size(); i++) {
    Vector2cd tmp = ac.Resultant(node[i], {&inc}).Bv();
    w(i) = tmp(0);
    t(i) = tmp(1);
  }
  VectorXcd coeff(ac.ScatterCoeff());
  VectorXcd ww = zw * coeff, tt = zt * coeff;

  std::ofstream disp_file("disp.dat"), trac_file("trac.dat");
  ApproxVectRv(w, ww, 1e-40, 0, true, disp_file);
  ApproxVectRv(t, tt, 1e-40, 0, true, trac_file);
}
TEST_F(PeriodicTest, DISABLED_Eigen_Disp) {
  input::Solution input{path("input2.txt")};
  Matrix matrix(input);
  AssemblyConfig<AP> ac(input.assembly_config()[0], &matrix);
  ac.Boundary().ReverseEdge();

  MatrixXcd z1(2 * (ac.Edge(0).size() + ac.Edge(1).size()), ac.NumCoeff());
  MatrixXcd z2(2 * (ac.Edge(2).size() + ac.Edge(3).size()), ac.NumCoeff());

  // MatrixXcd z1((ac.Edge(0).size() + ac.Edge(1).size()), ac.NumCoeff());
  // MatrixXcd z2((ac.Edge(2).size() + ac.Edge(3).size()), ac.NumCoeff());
  // z1 << ac.CylinEDMat(ac.Edge(0)), ac.CylinEDMat(ac.Edge(1));
  // z2 << ac.CylinEDMat(ac.Edge(2)), ac.CylinEDMat(ac.Edge(3));

  z1 << ac.ResBvMat(ac.Edge(0)), ac.ResBvMat(ac.Edge(1));
  z2 << ac.ResBvMat(ac.Edge(2)), ac.ResBvMat(ac.Edge(3));
  // z1 << ac.ScatterBvMat(ac.Edge(0)), ac.ScatterBvMat(ac.Edge(1));
  // z2 << ac.ScatterBvMat(ac.Edge(2)), ac.ScatterBvMat(ac.Edge(3));
  // z1 << ac.CylinEBMat(ac.Edge(0)), ac.CylinEBMat(ac.Edge(1));
  // z2 << ac.CylinEBMat(ac.Edge(2)), ac.CylinEBMat(ac.Edge(3));
  for (long i = 1; i < z1.rows(); i += 2) z1.row(i) *= -1;

  MatrixXcd za1(z1), zg1(z1), za2(z2), zg2(z2);
  for (long i = 0; i < z1.rows(); i++) {
    dcomp p1 = z1.row(i).array().mean();
    za1.row(i) /= p1;
    za2.row(i) /= p1;
    dcomp p2 = GeometricMean(z1.row(i).array());
    zg1.row(i) /= p2;
    zg2.row(i) /= p2;
  }
  std::cout << z1.rows() << "  " << z1.cols() << std::endl;
  std::cout << "z1: " << ConditionNum(z1) << std::endl;
  std::cout << z1.jacobiSvd(40).singularValues() << std::endl;
  std::cout << "z2: " << ConditionNum(z2) << std::endl;
  std::cout << "Befo: " << ConditionNum(z1) << std::endl;
  std::cout << "Alge: " << ConditionNum(za1) << std::endl;
  std::cout << "Geom: " << ConditionNum(zg1) << std::endl;

  MatrixXcd z(z1.rows() + z2.rows(), z1.cols());
  z << z1, z2;
  writeMatrix(z, "z.dat");
  writeMatrix(z1, "z1");
  writeMatrix(z2, "z2");

  // MatrixXcd z11 = z.block(0, 0, z.rows() / 2, z.cols() / 2);
  // MatrixXcd z22 =
  //     z.block(z.rows() / 2, z.cols() / 2, z.rows() / 2, z.cols() / 2);

  // writeMatrix(z11, "z11.dat");
  // writeMatrix(z22, "z22.dat");

  std::cout << "z: " << ConditionNum(z) << std::endl;
  // std::cout << "z11: " << ConditionNum(z11) << std::endl;
  // std::cout << "z22: " << ConditionNum(z22) << std::endl;


  // MatrixXcd z(z1.rows() + z2.rows(), z1.cols());
  // z << z1, z2;
  // MatrixXcd zw(z.rows() / 2, z.cols()), zt(z.rows() / 2, z.cols());
  // for (long i = 0; i < zw.rows(); i++)
  //   zw.row(i) = z.row(2 * i), zt.row(i) = z.row(2 * i + 1);
  // MatrixXcd zaw(zw), zgw(zw), zat(zt), zgt(zt);
  // for (long i = 0; i < zw.rows(); i++) {
  //   dcomp p1 = zw.row(i).array().mean();
  //   zaw.row(i) /= p1;
  //   zat.row(i) /= p1;
  //   dcomp p2 = GeometricMean(zw.row(i).array());
  //   zgw.row(i) /= p2;
  //   zgt.row(i) /= p2;
  // }
  // std::cout << "zw: " << std::endl;
  // std::cout << zw.rows() << "  " << zw.cols() << std::endl;
  // std::cout << "Befo: " << ConditionNum(zw) << std::endl;
  // std::cout << "Alge: " << ConditionNum(zaw) << std::endl;
  // std::cout << "Geom: " << ConditionNum(zgw) << std::endl;
}

}  // namespace test

}  // namespace mss
