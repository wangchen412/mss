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
#include <sstream>
#include "../src/core/Solution.h"
#include "../src/post/Output.h"

using namespace mss;

template <typename T, int N>
class MultiBEM {
 public:
  MultiBEM(double omega, const Material& inhomo_mat,
           const Material& matrix_mat, const std::vector<PosiVect>& ps)
      : omega_(omega),
        matrix_(matrix_mat, omega),
        inhomo_(inhomo_mat, omega),
        n_(ps.size() / 4),
        b_(n_ + 1),
        nn_(n_ + 1) {
    in_ = new IncidentPlaneSH(matrix_);

    b_[0] =
        new Boundary<T, N>(10 * matrix_.KT(), ps, &matrix_, INCL_RECT, true);
    nn_[0] = b_[0]->NumNode();

    std::vector<PosiVect> positions{{-1, 1}, {1, -1}};

    for (int i = 0; i < n_; i++) {
      b_[i + 1] = new Boundary<T, N>(
          10 * matrix_.KT(),
          std::vector<PosiVect>(ps.begin() + i * 4, ps.begin() + i * 4 + 4),
          &inhomo_, INCL_RECT);
      nn_[i + 1] = b_[i + 1]->NumNode();
    }

    Eigen::MatrixXcd dtn(nn_[0], nn_[0]);
    dtn.setZero();
    long mm = 0;
    for (int i = 1; i <= n_; i++) {
      dtn.block(mm, mm, nn_[i], nn_[i]) = b_[i]->DtN();
      mm += nn_[i];
    }

    Eigen::VectorXcd w = (b_[0]->MatrixH() + b_[0]->MatrixG() * dtn)
                             .lu()
                             .solve(in_->EffectDv(b_[0]->Node()));

    Eigen::VectorXcd bv(nn_[0] * 2);
    mm = 0;
    for (int i = 1; i <= n_; i++) {
      Eigen::VectorXcd v = b_[i]->DispToEffect() * w.segment(mm, nn_[i]);
      b_[i]->SetBv(v);
      bv.segment(mm * 2, nn_[i] * 2) = v;
      mm += nn_[i];
    }

    for (long i = 1; i < bv.size(); i += 2) bv[i] *= -1;
    b_[0]->SetBv(bv);
  }

  ~MultiBEM() {
    delete in_;
    for (auto& i : b_) delete i;
  }

  T Resultant(const CS* cs) const {
    for (int i = 1; i <= n_; i++) {
      if (b_[i]->Contains(cs) == 1) return b_[i]->Effect(cs);
      if (b_[i]->Contains(cs) == 0) return b_[i]->Effect(cs, true);
    }
    return in_->Effect(cs) + b_[0]->Effect(cs);
  }

  double Frequency() const { return omega_; }
  Material material(const CS*) const { return matrix_.Material(); }
  const Matrix* matrix() const { return &matrix_; }

 private:
  double omega_;
  Matrix matrix_, inhomo_;
  int n_;
  std::vector<Boundary<T, N>*> b_;
  std::vector<long> nn_;
  Incident<T>* in_;
};

class FiberRes {
 public:
  FiberRes(const Fiber<AP>& fiber) : fiber(fiber) {}

  StateAP Resultant(const CS* cs) const {
    if (fiber.Contains(cs))
      return fiber.Inner(cs);
    else
      return fiber.Scatter(cs) + fiber.Pseudo(cs);
  }

  double Frequency() const { return 1; }
  Material material(const CS* a) const { return fiber.material(a); }

 private:
  const Fiber<AP>& fiber;
};

void Recover(const VectorXcd& bv) {
  const Material steel(7670, 116e9, 84.3e9), lead(11400, 36e9, 8.43e9);
  double omega = 16576.24319112025;
  Matrix matrix(steel, omega);
  auto fc = new FiberConfig<AP>("1", 14, 200, 0.06, lead, &matrix);
  InhomoPtrs<AP> fibers;
  fibers.push_back(new Fiber<AP>(fc, {0.1, 0.1}));
  auto ac = new AssemblyConfig<AP>("1", fibers, 0.2, 0.2, 500, &matrix);

  // std::cout << ac->Boundary().NumNode() << std::endl;
  // std::ofstream file("ac.txt");
  // for (auto& i : ac->Boundary().Node()) file << i->PositionGLB() <<
  // std::endl; file.close();

  VectorXcd coeff = ac->ResBvMat(ac->Node()).jacobiSvd(40).solve(bv);
  // std::cout << coeff << std::endl;

  Fiber<AP> f(fc);
  f.SetCoeff(coeff);
  FiberRes fr(f);
  post::Area<AP>(&fr, {-0.1, 0.1}, {0.1, -0.1}, 500, 500, "rec").Write();
  post::Circle<AP>(&fr, {0, 0}, 0.0601, 1000, "rec").Write();
  post::Line<AP>(&fr, {-0.1, 0}, {0.1, 0}, 1000, "rec").Write();
}

int main() {
  const Material steel(7670, 116e9, 84.3e9), lead(11400, 36e9, 8.43e9);
  Material eff_mat(0.803696 * 11400, 1, 0.609922 * 84.3e9);

  Eigen::MatrixXd ps1(2, 4);
  double R = 1, a = 19 * 0.2 + 0.08, b = 9 * 0.2 + 0.08, d = 0.04;
  ps1 << R - a - d, R - a - d, R + d, R + d, d - R, -R - b - d, -R - b - d,
      d - R;

  Eigen::MatrixXd rot(2, 2);
  rot << cos(pi2 / 3), -sin(pi2 / 3), sin(pi2 / 3), cos(pi2 / 3);
  Eigen::MatrixXd ps2 = rot * ps1;
  Eigen::MatrixXd ps3 = rot * ps2;

  std::vector<PosiVect> ps(12);
  for (int i = 0; i < 4; i++) {
    ps[i].x = ps1(0, i);
    ps[i].y = ps1(1, i);
    ps[i + 4].x = ps2(0, i);
    ps[i + 4].y = ps2(1, i);
    ps[i + 8].x = ps3(0, i);
    ps[i + 8].y = ps3(1, i);
  }

  double omega = 16576.24319112025;

  MultiBEM<AP, 10> s(omega, eff_mat, steel, ps);
  post::Area<AP>(&s, {-0.94, -1.74}, {-0.74, -1.94}, 500, 500, "bem").Write();
  post::Circle<AP>(&s, {-0.84, -1.84}, 0.0601, 1000, "bem").Write();
  post::Line<AP>(&s, {-0.94, -1.84}, {-0.74, -1.84}, 1000, "bem").Write();

  // post::Circle<AP>(&s, {0, 0}, 0.5155, 1000, "1").Write();
  // post::Circle<AP>(&s, {0, 0}, 2.5773, 1000, "2").Write();
  // post::Circle<AP>(&s, {0, 0}, 4.6392, 1000, "3").Write();

  Boundary<AP, 4> box(500.00001, {{-0.94, -1.74}, {-0.74, -1.94}},
                      s.matrix());

  // std::ofstream file("box.txt");
  // for (auto& i : box.Node()) file << i->PositionGLB() << std::endl;
  // file.close();

  VectorXcd bv(800);
  for (int i = 0; i < 400; i++)
    bv.segment<2>(i * 2) = s.Resultant(box.Node(i)).Bv();

  Recover(bv);

  return 0;
}
