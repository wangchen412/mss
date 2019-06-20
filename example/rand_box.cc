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

#include "../src/post/Output.h"
#include "../src/post/check/Continuity.h"

using namespace mss;

class Mismatch {
 public:
  Mismatch(double r, double omega, const Eigen::VectorXcd& w,
           const Eigen::VectorXcd& t, const Material& m0)
      : omega_(omega), r_(r), w_(w), t_(t), m0_(m0) {}

  double operator()(const Eigen::Vector4d& r) const {
    Matrix matrix(m0_.mul_comp(r), omega_);
    Boundary<AP, 4> b(100, {{-r_, r_}, {r_, -r_}}, &matrix);
    return (b.MatrixH() * w_ - b.MatrixG() * t_).norm();
  }

 private:
  double omega_, r_;
  Eigen::VectorXcd w_, t_;
  const Material m0_;
};

Eigen::Vector4d box_homo(double r, double x, double y, const Solution<AP>& s,
                         const Eigen::Vector4d& x0) {
  Boundary<AP, 4> b(100, {{x - r, y + r}, {x + r, y - r}}, s.Matrix());
  Eigen::VectorXcd w(b.NumNode()), t(b.NumNode());
  std::vector<StateAP> v(b.Node().size());

#ifdef NDEBUG
#pragma omp parallel for
#endif
  for (size_t k = 0; k < b.NumNode(); k++) {
    v[k] = s.Resultant(b.Node(k));
    w(k) = v[k].Bv()(0);
    t(k) = v[k].Bv()(1);
  }

  Mismatch f(r, s.Frequency(), w, t, {{11400, 11400}, 0, {84e9, 84e9}});
  return BasinHopping(0, 0, f, x0);
}

int main() {
  Solution<AP> s{input::Solution("input.txt")};
  std::ifstream coeff_in("coeff.txt");
  s.ReadCoeff(coeff_in);
  coeff_in.close();

  Eigen::Vector4d x0;
  x0 << 0.8, 0, 0.6, 0;
  std::ofstream file("boxes.dat");
  for (int m = 0; m <= 35; m++) {
    for (int i = -2; i <= 2; i += 2)
      for (int j = -2; j <= 2; j += 2) {
        file << box_homo(0.3 + 0.02 * m, j, i, s, x0).transpose() << "\t";
        file.flush();
        std::cout << i << "  " << j << std::endl;
      }
    file << std::endl;
  }
  file.close();

  return 0;
}
