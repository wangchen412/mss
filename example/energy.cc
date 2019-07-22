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

#include <fstream>
#include <string>
#include "../src/post/Output.h"
#include "../src/post/check/Continuity.h"

using namespace mss;

class solution {
 public:
  solution(double omega, double theta) : matrix_(matrix_mat_, omega) {
    fc = new FiberConfig<AP>("1", 14, 200, 0.06, inhomo_mat_, &matrix_);
    fibers.push_back(new Fiber<AP>(fc, {0.1, 0.1}));
    ac = new AssemblyConfig<AP>("1", fibers, 0.2, 0.2, 80000, &matrix_);
    ac->Boundary().ReverseEdge();

    MatrixXcd z1(2 * (ac->Edge(0).size() + ac->Edge(1).size()),
                 ac->NumCoeff());
    MatrixXcd z2(2 * (ac->Edge(2).size() + ac->Edge(3).size()),
                 ac->NumCoeff());
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
    MatrixXcd z = z2 - PhaseShift(ee * pi, theta, z1.rows()) * z1;
    MatrixXcd zp = z2 - PhaseShiftDiff(ee * pi, theta, z1.rows()) * z1;
    VectorXcd xx = NewtonEigen(z, zp);

    f = new Fiber<AP>(fc);
    f->SetCoeff(xx);
  }
  StateAP Resultant(const CS* cs) const {
    if (f->Contains(cs))
      return f->Inner(cs);
    else
      return f->Scatter(cs) + f->Pseudo(cs);
  }

  const Material& material(const CS* cs) const {
    if (f->Contains(cs))
      return inhomo_mat_;
    else
      return matrix_mat_;
  }
  double Frequency() const { return matrix_.Frequency(); }
  const Matrix* matrix() const { return &matrix_; }

 private:
  Material matrix_mat_{7670, 116e9, 84.3e9};
  Material inhomo_mat_{11400, 36e9, 8.43e9};
  Matrix matrix_;
  FiberConfig<AP>* fc;
  InhomoPtrs<AP> fibers;
  AssemblyConfig<AP>* ac;
  Fiber<AP>* f;
};

class expo {
 public:
  expo(double omega, const post::Area<AP>& area)
      : omega(omega), area(area), N(area.Points().size()) {
    MatrixXcd A(N, 3);
    VectorXcd b(N);
    for (long i = 0; i < N; i++) {
      PosiVect p = area.Point(i)->PositionGLB();
      A(i, 0) = p.x;
      A(i, 1) = p.y;
      A(i, 2) = 1;
      b(i) = log(area.Point(i)->State().Displacement().x);
    }
    c = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
  }

  dcomp Displacement(const PosiVect& p) const {
    Eigen::VectorXcd A(3);
    A << p.x, p.y, 1;
    return exp(A.dot(c));
  }

  StateAP Resultant(const CS* cs) const {
    PosiVect p = cs->PositionGLB();
    return StateAP(Displacement(p), dcomp(0), dcomp(0));
  }
  double Frequency() const { return 0; }
  Material material(const CS*) const { return matrix_mat_; }

  double rho() {
    double ww = 0;
    for (long i = 0; i < N; i++)
      ww += pow(std::abs(Displacement(area.Point(i)->PositionGLB())), 2);
    ww /= N;
    return area.KineticEnergyDensity() / (ww / 2 * omega * omega);
  }

  double mu() {
    double ww = 0;
    for (long i = 0; i < N; i++)
      ww += pow(std::abs(Displacement(area.Point(i)->PositionGLB())), 2);
    ww /= N;
    ww *= pow(std::abs(c(0)), 2) + pow(std::abs(c(1)), 2);
    ww /= 2;

    return area.StrainEnergyDensity() / ww;
  }

 private:
  double omega;
  const post::Area<AP>& area;
  long N;
  VectorXcd c;
  Material matrix_mat_{7670, 116e9, 84.3e9};
};
class plane {
 public:
  template <typename S>
  plane(const S* sol, const post::Area<AP>& area)
      : omega(sol->Frequency()), area(area), N(area.Points().size()) {
    Eigen::MatrixXcd A(N, 3);
    Eigen::VectorXcd b(N);
    for (long i = 0; i < N; i++) {
      PosiVect p = area.Point(i)->PositionGLB();
      A(i, 0) = p.x;
      A(i, 1) = p.y;
      A(i, 2) = 1;
      b(i) = area.Point(i)->State().Displacement().x;
    }
    c = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

    bb = new Boundary<AP>(1500, {{-0.1, 0.1}, {0.1, -0.1}}, sol->matrix());
    w.resize(bb->NumNode());
    t.resize(bb->NumNode());
    for (size_t i = 0; i < bb->NumNode(); i++) {
      StateAP s = sol->Resultant(bb->Node(i));
      w(i) = s.Displacement().x;
      t(i) = s.Stress().x;
    }
  }

  dcomp Displacement(const PosiVect& p) const {
    Eigen::VectorXcd A(3);
    A << p.x, p.y, 1;
    return A.dot(c);
  }
  StateAP Resultant(const CS* cs) const {
    PosiVect p = cs->PositionGLB();
    return StateAP(Displacement(p), dcomp(0), dcomp(0));
  }
  double Frequency() const { return 0; }
  Material material(const CS*) const { return matrix_mat_; }

  double rho() {
    double ww = 0;
    for (long i = 0; i < N; i++)
      ww += pow(std::abs(Displacement(area.Point(i)->PositionGLB())), 2);
    ww /= N;
    return area.KineticEnergyDensity() / (ww / 2 * omega * omega);
  }
  double mu() {
    VectorXcd w2(bb->NumNode());
    for (size_t i = 0; i < bb->NumNode(); i++)
      w2(i) = Resultant(bb->Node(i)).Displacement().x;

    VectorXcd dw = w2 - w;
    double wt = dw.dot(t).real() / t.rows() * 20;

    return (area.StrainEnergyDensity() + wt) /
           (pow(std::abs(c(0)), 2) + pow(std::abs(c(1)), 2)) * 2;
  }

 private:
  double omega;
  const post::Area<AP>& area;
  Boundary<AP>* bb;
  VectorXcd w, t;
  long N;
  VectorXcd c;
  Material matrix_mat_{7670, 116e9, 84.3e9};
};
class box {
 public:
  template <typename S>
  box(const S* sol, const post::Area<AP>& area)
      : omega(sol->Frequency()), area(area) {
    bb = new Boundary<AP>(500, {{-0.1, 0.1}, {0.1, -0.1}}, sol->matrix());
    MatrixXcd A(bb->Node().size(), 3);
    VectorXcd b(bb->Node().size());

    for (size_t i = 0; i < bb->Node().size(); i++) {
      PosiVect p = bb->Node(i)->PositionGLB();
      A(i, 0) = p.x;
      A(i, 1) = p.y;
      A(i, 2) = 1;
      b(i) = sol->Resultant(bb->Node(i)).Displacement().x;
    }
    c = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
  }

  dcomp Displacement(const PosiVect& p) const {
    Eigen::VectorXcd A(3);
    A << p.x, p.y, 1;
    return A.dot(c);
  }
  StateAP Resultant(const CS* cs) const {
    PosiVect p = cs->PositionGLB();
    return StateAP(Displacement(p), dcomp(0), dcomp(0));
  }
  double Frequency() const { return 0; }
  Material material(const CS*) const { return matrix_mat_; }

  double rho() {
    double ww = 0;
    for (size_t i = 0; i < area.Points().size(); i++)
      ww += pow(std::abs(Displacement(area.Point(i)->PositionGLB())), 2);
    ww /= area.Points().size();
    return area.KineticEnergyDensity() / (ww / 2 * omega * omega);
  }
  double mu() {
    return area.StrainEnergyDensity() /
           (pow(std::abs(c(0)), 2) + pow(std::abs(c(1)), 2)) * 2;
  }

 private:
  double omega;
  const post::Area<AP>& area;
  Boundary<AP>* bb;
  Eigen::VectorXcd c;
  Material matrix_mat_{7670, 116e9, 84.3e9};
};

int main() {
  std::ofstream file("energy_expo.txt");

  int N = 20;
  double fmax = 4800;
  double fmin = 500;
  double df = (fmax - fmin) / N;

  for (int i = 0; i < 20; i++) {
    std::cout << i << std::endl;
    double omega = (fmin + df * i) * pi2;
    solution s(omega, 0);
    post::Area<AP> area(&s, {-0.1, 0.1}, {0.1, -0.1}, 100, 100, "eigen");
    // plane a(&s, area);
    expo a(omega, area);
    file << fmin + df * i << "\t" << a.rho() << "\t" << a.mu() << std::endl;
  }

  file.close();

  return 0;
}

// int main() {
//   double omega = 16576.24319112025;
//   solution s(omega, pi / 5);
//   post::Area<AP> area(&s, {-0.1, 0.1}, {0.1, -0.1}, 300, 300, "eigen");
//   area.Write();
//   expo ex(omega, area);
//   post::Area<AP>(&ex, {-0.1, 0.1}, {0.1, -0.1}, 300, 300, "expo").Write();
//   return 0;
// }
