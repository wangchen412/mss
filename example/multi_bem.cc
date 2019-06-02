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

 private:
  double omega_;
  Matrix matrix_, inhomo_;
  int n_;
  std::vector<Boundary<T, N>*> b_;
  std::vector<long> nn_;
  Incident<T>* in_;
};

int main() {
  Material steel(7670, 116e9, 84.3e9), lead(11400, 36e9, 8.43e9);

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

  MultiBEM<AP, 10> s(1657.624 * 2, steel, lead, ps);
  post::Area<AP>(&s, {-6, 6}, {6, -6}, 300, 300, "bem").Write();

  return 0;
}
