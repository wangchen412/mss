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

class periodic {
 public:
  periodic(double omega, double theta) : matrix_(matrix_mat_, omega) {
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
    double ee = MinDet(z2, z1, theta);
    k = ee * pi / 0.2;
    kx = k * cos(theta);
    ky = k * sin(theta);

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

  double k, kx, ky;

 private:
  Material matrix_mat_{7670, 116e9, 84.3e9};
  Material inhomo_mat_{11400, 36e9, 8.43e9};
  Matrix matrix_;
  FiberConfig<AP>* fc;
  InhomoPtrs<AP> fibers;
  AssemblyConfig<AP>* ac;
  Fiber<AP>* f;
};
class ms {
 public:
  ms(double omega, double angle) {
    input::Solution in("input.txt");
    in.update_frequency(omega);
    in.update_incident_angle(angle);
    s = new Solution<AP>(in);
    s->Solve();

    int P = 100;
    post::AreaAP area(s, {-0.5, 0.5}, {0.5, -0.5}, P, P);
    MatrixXcd A(P * P, 3);
    VectorXcd b(P * P);
    for (long i = 0; i < P * P; i++) {
      PosiVect p = area.Point(i)->PositionGLB();
      A(i, 0) = p.x;
      A(i, 1) = p.y;
      A(i, 2) = 1;
      b(i) = log(area.Point(i)->State().Displacement().x);
    }
    VectorXcd c =
        A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
    kx = (c(0) / ii).real();
    ky = (c(1) / ii).real();
    k = sqrt(kx * kx + ky * ky);
    std::cout << "kx from fitting: " << c(0) << std::endl;
    std::cout << "ky from fitting: " << c(1) << std::endl;
    std::cout << "k from fitting: " << k << std::endl;
    std::cout << "Amplitude from fitting: " << exp(c(2)) << std::endl;
  }

  Material material(const CS* cs) const { return s->material(cs); }
  StateAP Resultant(const CS* cs) const { return s->Resultant(cs); }
  double Frequency() const { return s->Frequency(); }

  double k, kx, ky;

 private:
  Solution<AP>* s;
};

class homo {
 public:
  template <typename S>
  homo(const S* sol, const post::Area<AP>& area)
      : omega(sol->Frequency()),
        area(area),
        N(area.Points().size()),
        k(sol->k),
        kx(sol->kx),
        ky(sol->ky) {
    for (auto& i : area.Points())
      A += sol->Resultant(i->LocalCS()).Displacement().x *
           exp(-ii * kx * i->PositionGLB().x - ii * ky * i->PositionGLB().y);
    A /= N;
    std::cout << "Amp from averaging: " << A << std::endl;
  }

  dcomp Displacement(const PosiVect& p) const {
    return A * exp(ii * kx * p.x + ii * ky * p.y);
  }

  double rho() {
    return area.KineticEnergyDensity() * 2 / pow(std::abs(A) * omega, 2);
  }

  double mu() { return pow(omega / k, 2) * rho(); }

 private:
  double omega;
  const post::Area<AP>& area;
  long N;
  dcomp A{0};
  double k, kx, ky;
};

template <typename T>
Eigen::Vector2d homo_ang(double omega, double angle) {
  T s(omega, angle);
  post::Area<AP> area(&s, {-0.1, 0.1}, {0.1, -0.1}, 100, 100, "eigen");
  homo p(&s, area);
  return Eigen::Vector2d(p.rho(), p.mu());
}

template <typename T>
Eigen::Vector2d homo_iso(double omega, int N = 45) {
  Eigen::MatrixXd rst(2, N);
  for (long i = 0; i < N; i++)
    rst.col(i) = homo_ang<T>(omega, pi / 4 / N * i);
  return Eigen::Vector2d(rst.row(0).mean(), rst.row(1).mean());
}

int main() {
  std::ofstream file("energy_test.txt");

  int N = 1;
  double fmax = 4878;
  double fmin = 4800;
  double df = (fmax - fmin) / N;

  for (int i = 0; i < N; i++) {
    std::cout << i << std::endl;
    double omega = (fmin + df * i) * pi2;
    file << fmin + df * i << "\t" << homo_iso<ms>(omega, 1).transpose()
         << std::endl;
  }
  file.close();

  // solution ss(4792.8 * pi2, pi / 4);
  // post::Area<AP> aa(&ss, {-0.1, 0.1}, {0.1, -0.1}, 300, 300, "eigen");
  // plane p(&ss, aa);
  // std::cout << 4792.8 * pi2 << "\t" << p.rho() << "\t" << p.mu() <<
  // std::endl;

  return 0;
}
