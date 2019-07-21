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

template <typename T, int N>
class MultiBEM {
 public:
  MultiBEM(double omega, const Material& matrix_mat,
           const MatrixXcd& inhomo_DtN, const std::vector<PosiVect>& ps)
      : s_(inhomo_DtN.rows()),
        omega_(omega),
        matrix_(matrix_mat, omega),
        inhomo_DtN_(inhomo_DtN),
        inhomo_DtE_(s_ * 2, s_),
        n_(ps.size() / 4) {
    MatrixXcd I = MatrixXcd::Identity(s_, s_);
    for (long i = 0; i < s_; i++) {
      inhomo_DtE_.row(i * 2) = I.row(i);
      inhomo_DtE_.row(i * 2 + 1) = inhomo_DtN_.row(i);
    }

    in_ = new IncidentPlaneSH(matrix_);

    b_ = new Boundary<T, N>(300, ps, &matrix_, INCL_RECT, true);
    nn_ = b_->NumNode();

    std::cout << "Total node number: " << nn_ << std::endl;

    Eigen::MatrixXcd dtn(nn_, nn_);
    dtn.setZero();
    long mm = 0;
    for (int i = 1; i <= n_; i++) {
      dtn.block(mm, mm, s_, s_) = inhomo_DtN_;
      mm += s_;
    }

    Eigen::VectorXcd w = (b_->MatrixH() + b_->MatrixG() * dtn)
                             .lu()
                             .solve(in_->EffectDv(b_->Node()));

    Eigen::VectorXcd bv(nn_ * 2);
    mm = 0;
    for (int i = 1; i <= n_; i++) {
      Eigen::VectorXcd v = inhomo_DtE_ * w.segment(mm, s_);
      bv.segment(mm * 2, s_ * 2) = v;
      mm += s_;
    }

    for (long i = 1; i < bv.size(); i += 2) bv[i] *= -1;
    b_->SetBv(bv);
  }

  ~MultiBEM() {
    delete in_;
    delete b_;
  }

  T Resultant(const CS* cs) const { return in_->Effect(cs) + b_->Effect(cs); }

 private:
  long s_;
  double omega_;
  Matrix matrix_;
  MatrixXcd inhomo_DtN_, inhomo_DtE_;
  int n_;
  long nn_;
  Boundary<T, N>* b_;
  Incident<T>* in_;
};

MatrixXcd DtN() {
  input::Solution input("input3.txt");
  Matrix matrix(input);
  IncidentPlaneSH inc(matrix, input.incident()[0]);

  AssemblyConfig<AP> ac(input.assembly_config("low"), &matrix);
  CSCPtrs node = ac.EdgeNode();

  std::cout << node.size() << std::endl;

  MatrixXcd z(ac.ResBvMat(node));
  MatrixXcd zw(z.rows() / 2, z.cols()), zt(z.rows() / 2, z.cols());
  for (long i = 0; i < zw.rows(); i++) {
    zw.row(i) = z.row(2 * i);
    zt.row(i) = z.row(2 * i + 1);
  }
  MatrixXcd dtn = zt * PseudoInverse(zw);

  AssemblyConfig<AP> ac2(input.assembly_config("high"), &matrix);
  ac2.DSolve({&inc});
  VectorXcd w(node.size()), t(node.size());
  for (size_t i = 0; i < node.size(); i++) {
    Vector2cd tmp = ac2.Resultant(node[i], {&inc}).Bv();
    w(i) = tmp(0);
    t(i) = tmp(1);
  }

  VectorXcd tt = dtn * w;

  std::ofstream file("miss.dat");
  ApproxVectRv(t, tt, 1e-30, 0, true, file);
  file.close();

  std::cout << "DtN computed." << std::endl;

  return dtn;
}

MatrixXcd DtN2(double omega, const std::vector<PosiVect>& ps) {
  Matrix m({0.803696 * 11400, 1, 0.609922 * 84.3e9}, omega);
  Boundary<AP, 10> b(300, ps, &m, INCL_RECT);
  return b.DtN();
}

int main() {
  // DtN();
  std::vector<PosiVect> ps(4);
  ps[0].x = -0.5;
  ps[0].y = 0.5;
  ps[1].x = -0.5;
  ps[1].y = -0.5;
  ps[2].x = 0.5;
  ps[2].y = -0.5;
  ps[3].x = 0.5;
  ps[3].y = 0.5;

  Material steel(7670, 116e9, 84.3e9);
  // MultiBEM<AP, 10> s(34557.5192, steel, DtN2(34557.5192, ps), ps);
  MultiBEM<AP, 4> s(34557.5192, steel, DtN(), ps);

  std::cout << "Solved. Preparing output." << std::endl;

  post::Area<AP> a(&s, {-4, 4}, {4, -4}, 300, 300, "bem");
  a.Write();

  return 0;
}
