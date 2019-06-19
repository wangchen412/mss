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
    Boundary<AP, 4> b(500, {{0, 0}, {0, r_}}, &matrix, CIRCULAR);
    return (b.MatrixH() * w_ - b.MatrixG() * t_).norm();
  }

 private:
  double omega_, r_;
  Eigen::VectorXcd w_, t_;
  const Material m0_;
};

Eigen::Vector4d circ_homo(double r, const Solution<AP>& s) {
  Boundary<AP, 4> b(500, {{0, 0}, {0, r}}, s.Matrix(), CIRCULAR);
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
  return BasinHopping(2, 2, f, Eigen::Vector4d::Ones());
}

int main() {
  Solution<AP> s{input::Solution("input.txt")};
  std::ifstream coeff_in("coeff.txt");
  s.ReadCoeff(coeff_in);
  coeff_in.close();

  std::ofstream file("circ_homo.dat");
  for (int i = 0; i < 41; i++)
    file << circ_homo(0.1 + 0.005 * i, s) << std::endl;
  file.close();

  return 0;
}
