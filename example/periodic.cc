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

#include "mismatch.h"

using namespace mss;

class FiberRes {
 public:
  FiberRes(const Fiber<AP>& fiber) : fiber(fiber) {}

  StateAP Resultant(const CS* cs) const {
    if (fiber.Contains(cs))
      return fiber.Inner(cs);
    else
      return fiber.Scatter(cs) + fiber.Pseudo(cs);
  }

 private:
  const Fiber<AP>& fiber;
};

void bv(double omega, double theta, VectorXcd& w, VectorXcd& t) {
  const Material steel(7670, 116e9, 84.3e9), lead(11400, 36e9, 8.43e9);

  Matrix matrix(steel, omega);
  auto fc = new FiberConfig<AP>("1", 14, 200, 0.06, lead, &matrix);
  InhomoPtrs<AP> fibers;
  fibers.push_back(new Fiber<AP>(fc, {0.1, 0.1}));
  auto ac = new AssemblyConfig<AP>("1", fibers, 0.2, 0.2, 80000, &matrix);
  ac->Boundary().ReverseEdge();

  std::cout << ac->Edge(0).size() << std::endl;

  MatrixXcd z1(2 * (ac->Edge(0).size() + ac->Edge(1).size()), ac->NumCoeff());
  MatrixXcd z2(2 * (ac->Edge(2).size() + ac->Edge(3).size()), ac->NumCoeff());
  z1 << ac->ResBvMat(ac->Edge(0)), ac->ResBvMat(ac->Edge(1));
  z2 << ac->ResBvMat(ac->Edge(2)), ac->ResBvMat(ac->Edge(3));
  for (long i = 1; i < z1.rows(); i += 2) z1.row(i) *= -1;
  for (long i = 0; i < z1.rows(); i++) {
    dcomp p = z1.row(i).array().mean();
    // dcomp p = GeometricMean(z1.row(i).array());
    z1.row(i) /= p;
    z2.row(i) /= p;
  }
  dcomp ee = MinDet(z2, z1, theta);
  std::cout << "Ka (pi): " << ee << std::endl;

  VectorXcd xx =
      NewtonEigen(z2 - PhaseShift(exp(ii * ee * pi), theta, z1.rows()) * z1);

  Fiber<AP> f(fc);
  f.SetCoeff(xx);

  // FiberRes fr(f);
  // post::Area<AP>(&fr, {-0.1, 0.1}, {0.1, -0.1}, 400, 400,
  //                std::to_string(xx.dot(yy).real()))
  //     .Write();

  Boundary<AP, 4> box(w.rows() / 0.8, {{-0.1, 0.1}, {0.1, -0.1}}, &matrix);
  std::vector<StateAP> v(box.Node().size());
#ifdef NDEBUG
#pragma omp parallel for
#endif
  for (size_t k = 0; k < box.NumNode(); k++) {
    v[k] = f.Scatter(box.Node(k)) + f.Pseudo(box.Node(k));
    w(k) = v[k].Bv()(0);
    t(k) = v[k].Bv()(1);
  }

  // std::ofstream file("wt_com.txt");
  // for (int i = 0; i < 400; i++) file << w(i) << "\t" << t(i) << std::endl;
  // file.close();
}
double read(int n, VectorXcd& w, VectorXcd& t) {
  auto file = std::ifstream(std::string("wt/") + std::to_string(n));
  std::string tmp;

  std::getline(file, tmp);
  std::stringstream s(tmp);
  double f;
  s >> f;

  int i = 0;
  while (std::getline(file, tmp)) {
    std::stringstream s(tmp);
    double x, y, wr, wi, tr, ti;
    s >> x >> y >> wr >> wi >> tr >> ti;
    w(i) = dcomp(wr, wi);
    t(i) = dcomp(tr, ti);
    ++i;
  }
  if (i != 400) std::cout << "error: " << i << std::endl;
  return f;
}
Eigen::VectorXd homo(double omega, const VectorXcd& w, const VectorXcd& t) {
  const Material norm_mat({11400, 11400}, 0, {84e9, 84e9});
  Eigen::VectorXd x0(4);
  x0.setOnes();
  Mismatch f(omega, w, t, norm_mat, 0.2, 0.2, 500);
  // return BasinHopping(2, 2, f, x0);
  return NelderMead(f, x0, &std::cout);
}

int main() {
  std::ifstream freq_file("freq.txt");
  std::ofstream out_file("mss_pbc.txt");
  std::vector<double> freq(80);
  std::string tmp;
  for (int i = 0; i < 80; i++) {
    std::getline(freq_file, tmp);
    std::stringstream s(tmp);
    s >> freq[i];
  }
  for (int i = 0; i < 80; i++) {
    VectorXcd w(400), t(400);
    w.setZero();
    t.setZero();

    for (int j = 0; j < 45; j++) {
      VectorXcd ww(400), tt(400);
      bv(freq[i] * pi2, pi / 4 / 45 * j, ww, tt);
      w += ww;
      t += tt;
    }

    out_file << freq[i] << "\t" << homo(freq[i], w, t).transpose()
             << std::endl;
  }
  out_file.close();
  return 0;
}
