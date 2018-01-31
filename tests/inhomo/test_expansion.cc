#include "../test.h"

namespace mss {

namespace test {

class ExpansionTest : public Test {
 protected:
  ExpansionTest() {
    for (size_t i = 0; i < np; i++)
      for (size_t j = 0; j < np; j++)
        points.push_back(new CS(-a + j * d, a - i * d));
    f2.SetCoeff(f2.CSolve({&in1}));
    // DFT();
    collocation();
  }
  virtual ~ExpansionTest() {
    for (auto i : points) delete i;
  }

  size_t np{300};
  double r{3e-3};
  double a{8 * r};
  double d = 2 * a / np;
  int N{300};
  CSCPtrs points;
  VectorXcd c{N};

  Material rubber{1300, 1.41908e9, 0.832e9}, lead{11400, 36.32496e9, 8.43e9};
  Matrix matrix{rubber, 1e6};
  IncidentPlaneSH in1{matrix, 0, 1e-5, 2};
  FiberConfig<AP> fc{"1", 40, 500, r, lead, &matrix};
  Fiber<AP> fiber{&fc, {0, 0}};
  Fiber<AP> f2{&fc, {a + r, 0}};

  StateAP mode(const CS* objCS, int n) const {
    return ModeT<AP>(fiber.LocalCS(), objCS,
                     EigenFunctor(Jn, n, matrix.KT(), fiber.Radius()),
                     matrix.Material());
  }
  // MatrixXcd colloMat(CSCPtrs objCSs, int N) const {
  //   MatrixXcd rst(objCSs.size(), 2 * N + 1);
  //   for (size_t i = 0; i < objCSs.size(); i++)
  //     for (int j = -N; j <= N; j++)
  //       rst(i, j + N) = mode(objCSs[i], j).Displacement().x;
  //   return rst;
  // }
  MatrixXcd colloMat(CSCPtrs objCSs) const {
    size_t P = objCSs.size();
    MatrixXcd rst(P, N);
    for (size_t i = 0; i < P; i++)
      for (int j = 0; j < N; j++)
        rst(i, j) =
            _planeWaveAP(objCSs[i], pi2 / N * j, &matrix).Displacement().x;
    // rst.block<2, 1>(i * 2, j) =
    //     _planeWaveAP(objCSs[i], pi2 / N * j, &matrix).Bv();

    return rst;
  }
  void collocation() {
    VectorXcd b(fiber.NumNode());
    for (size_t i = 0; i < fiber.NumNode(); i++)
      // b(i) = (in1.Effect(fiber.Node(i))).Displacement().x;
      b(i) = (in1.Effect(fiber.Node(i)) + f2.Scatter(fiber.Node(i)))
                 .Displacement()
                 .x;

    // VectorXcd b(fiber.NumBv());
    // for (size_t i = 0; i < fiber.NumNode(); i++)
    //   b.segment<2>(i * 2) =
    //       (in1.Effect(fiber.Node(i)) + f2.Scatter(fiber.Node(i))).Bv();

    // c = colloMat(fiber.Node()).lu().solve(b);

    c = PseudoInverse(colloMat(fiber.Node())) * b;

    // c = colloMat(fiber.Node())
    //         .jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV)
    //         .solve(b);
  }
  void DFT() {
    for (int n = -N; n <= N; n++) c(n + N) = coeff_DFT(n);
  }
  dcomp coeff_DFT(int N) {
    double k = matrix.KT();
    double r = fiber.Radius();
    double mu = matrix.Material().Mu();

    VectorXcd in(fiber.NumNode());
    for (size_t i = 0; i < fiber.NumNode(); i++)
      in(i) = (in1.Effect(fiber.Node(i)) + f2.Scatter(fiber.Node(i))).Bv()(0);
    return DFT_v(fiber.NumNode(), N).dot(in);
  }

  StateAP compute(const CS* objCS) const {
    // StateAP rst(objCS);
    // if (objCS->PositionGLB().Length() < r) return rst;
    // for (int n = -N; n <= N; n++) rst += mode(objCS, n) * c(n + N);
    // return rst;

    StateAP rst(objCS);
    for (long i = 0; i < c.size(); i++)
      rst += _planeWaveAP(objCS, pi2 / c.size() * i, &matrix) * c(i);
    return rst;
  }
};

TEST_F(ExpansionTest, Incident) {
  std::ofstream file("incident.txt");
  for (auto p : points) file << in1.Effect(p) + f2.Scatter(p) << std::endl;
  // for (auto p : points) file << in1.Effect(p) << std::endl;
  file.close();
}
TEST_F(ExpansionTest, Residue) {
  std::ofstream file("residue.txt");
  for (auto p : points)
    file << compute(p) - in1.Effect(p) - f2.Scatter(p) << std::endl;
  file.close();
}
TEST_F(ExpansionTest, Expansion) {
  std::ofstream file("expansion.txt");
  for (auto p : points) file << compute(p) << std::endl;
  file.close();
  // //std::cout << c << std::endl;
  MatrixXcd m = colloMat(fiber.Node());
  std::cout << "kr: " << matrix.KT() * r << std::endl;
  std::cout << m.rows() << "  " << m.cols() << std::endl;
  std::cout << m.jacobiSvd().rank() << std::endl;
  // std::cout << m.jacobiSvd().singularValues() << std::endl;
  // MatrixXcd mm = PseudoInverse(m);
  // std::cout << mm * m << std::endl;

  VectorXcd ref_w(fiber.NumNode()), com_w(fiber.NumNode());
  VectorXcd ref_t(fiber.NumNode()), com_t(fiber.NumNode());

  for (size_t i = 0; i < fiber.NumNode(); i++) {
    const CS* p = fiber.Node(i);
    Vector2cd ref = (in1.Effect(p) + f2.Scatter(p)).Bv();
    Vector2cd com = compute(p).Bv();
    ref_w(i) = ref(0);
    ref_t(i) = ref(1);
    com_w(i) = com(0);
    com_t(i) = com(1);
  }

  EXPECT_TRUE(ApproxVectRv(ref_w, com_w, 1e-3, 0, true));
  EXPECT_TRUE(ApproxVectRv(ref_t, com_t, 1e-3, 0, true));
}
TEST_F(ExpansionTest, DISABLED_BesselJ) {
  double kr = matrix.KT() * fiber.Radius();
  std::cout << "kr = \t" << kr << std::endl;

  for (int NN = 10; NN < 40; NN++) {
    std::cout << "N = " << NN << std::endl;
    std::cout << "Jn(kr) = \t" << Jn(NN, kr) << std::endl;
    std::cout << "Jn(6kr) = \t" << Jn(NN, 6 * kr) << std::endl;
    std::cout << "Ratio = \t" << Jn(NN, 6 * kr) / Jn(NN, kr) << std::endl;
    std::cout << "Coeff = \t" << coeff_DFT(NN) << std::endl;
    std::cout << "Contribution = \t"
              << coeff_DFT(NN) * Jn(NN, 6 * kr) / Jn(NN, kr) << std::endl;
    std::cout << "++++++++++++++++++++++++++++++++++++" << std::endl;
  }
}

}  // namespace test

}  // namespace mss
