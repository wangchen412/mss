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
#include "../src/post/check/Continuity.h"

using namespace mss;

class Mismatch {
 public:
  Mismatch(double omega, const Eigen::VectorXcd& w, const Eigen::VectorXcd& t,
           const Material& m0)
      : omega_(omega), w_(w), t_(t), m0_(m0) {}

  double operator()(const Eigen::Vector4d& r) const {
    Matrix matrix(m0_ * r, omega_);
    Boundary<AP, 4> b{500, {{-0.3, 0.3}, {0.3, -0.3}}, &matrix};
    return (b.MatrixH() * w_ - b.MatrixG() * t_).norm();
  }

 private:
  double omega_;
  Eigen::VectorXcd w_, t_;
  const Material m0_;
};

void compute_bv(double omega, Eigen::VectorXcd& w, Eigen::VectorXcd& t) {
  input::Solution in("input.txt");
  in.update_frequency(omega);
  Solution<AP> s(in);
  s.Solve();

  post::CC_Solution<AP> cc{&s};
  std::cout << mss_msg({"Maximum mismatch: ", std::to_string(cc.Max())})
            << std::endl;

  Boundary<AP, 4> b{500, {{-0.3, 0.3}, {0.3, -0.3}}, s.Matrix()};
#ifdef NDEBUG
#pragma omp parallel for
#endif
  for (size_t i = 0; i < 1200; i++) {
    VectorXcd bv = s.Resultant(b.Node(i)).Bv();
    w(i) = bv(0);
    t(i) = bv(1);
  }

  std::ofstream file("bv.dat");
  for (size_t i = 0; i < b.NumNode(); i++)
    file << setMaxPrecision << b.Node(i)->PositionGLB() << "\t"
         << b.Node(i)->AngleGLB() << "\t" << w(i) << "\t" << t(i)
         << std::endl;
  file.close();

  std::cout << mss_msg({"Boundary values computed."}) << std::endl;
}
void read_bv(const std::string& fn, Eigen::VectorXcd& w,
             Eigen::VectorXcd& t) {
  std::ifstream file(fn);
  std::string tmp;
  double x, y, ang;
  for (int i = 0; i < 1200; i++) {
    getline(file, tmp);
    std::stringstream ss(tmp);
    ss >> x >> y >> ang >> w(i) >> t(i);
  }
  std::cout << mss_msg({"Boundary values read from ", fn}) << std::endl;
}

int main(int argc, char** argv) {
  if (argc != 6) exit_error_msg({"Initial values and step size needed."});
  double omega = 16576.243191120248 * 1.95;

  Eigen::VectorXcd w(1200), t(1200);
  // if (argc == 2) compute_bv(omega, w, t);
  // if (argc == 3) read_bv(argv[2], w, t);
  read_bv("bv.dat", w, t);

  Mismatch f(omega, w, t, {{11400, 11400}, 0, {84e9, 84e9}});
  std::ofstream file("iterations.dat");
  GradientDescent(f, &file, {atof(argv[1]), atof(argv[2]), atof(argv[3]), atof(argv[4])}, atof(argv[5]));
  file.close();
  return 0;
}
