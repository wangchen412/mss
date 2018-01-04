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
  int N{25};
  CSCPtrs points;
  VectorXcd c{2 * N + 1};

  Material rubber{1300, 1.41908e9, 0.832e9}, lead{11400, 36.32496e9, 8.43e9};
  Matrix matrix{rubber, 1e6};
  IncidentPlaneSH in1{matrix, 0, 1, 2};
  FiberConfig<AP> fc{"1", 40, 1000, r, lead, &matrix};
  Fiber<AP> fiber{&fc, {0, 0}};
  Fiber<AP> f2{&fc, {a + r, 0}};

  StateAP mode(const CS* objCS, int n) const {
    return ModeT<AP>(fiber.LocalCS(), objCS,
                     EigenFunctor(Jn, n, matrix.KT(), fiber.Radius()),
                     matrix.Material());
  }
  MatrixXcd colloMat(CSCPtrs objCSs, int N) const {
    MatrixXcd rst(objCSs.size(), 2 * N + 1);
    for (size_t i = 0; i < objCSs.size(); i++)
      for (int j = -N; j <= N; j++)
        rst(i, j + N) = mode(objCSs[i], j).Displacement().x;
    return rst;
  }
  void collocation() {
    VectorXcd b(fiber.NumNode());
    for (size_t i = 0; i < fiber.NumNode(); i++)
      b(i) = (in1.Effect(fiber.Node(i)) + f2.Scatter(fiber.Node(i)))
                 .Displacement()
                 .x;
    c = colloMat(fiber.Node(), N)
            .jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV)
            .solve(b);
  }
  void DFT() {
    for (int n = -N; n <= N; n++) c(n + N) = coeff_DFT(n);
  }
  dcomp coeff_DFT(int N) {
    double k  = matrix.KT();
    double r  = fiber.Radius();
    double mu = matrix.Material().Mu();

    VectorXcd in(fiber.NumNode());
    for (size_t i = 0; i < fiber.NumNode(); i++)
      in(i) = (in1.Effect(fiber.Node(i)) + f2.Scatter(fiber.Node(i))).Bv()(0);
    return DFT_v(fiber.NumNode(), N).dot(in);
  }

  StateAP compute(const CS* objCS) const {
    StateAP rst(objCS);
    if (objCS->PositionGLB().Length() < r) return rst;
    for (int n = -N; n <= N; n++) rst += mode(objCS, n) * c(n + N);
    return rst;
  }
};

TEST_F(ExpansionTest, Incident) {
  std::ofstream file("incident.txt");
  for (auto p : points) file << in1.Effect(p) + f2.Scatter(p) << std::endl;
  file.close();
}
TEST_F(ExpansionTest, Expansion) {
  std::ofstream file("expansion.txt");
  for (auto p : points) file << compute(p) << std::endl;
  file.close();
  std::cout << c << std::endl;
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
