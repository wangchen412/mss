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

#ifndef MSS_EXAMPLE_MISMATCH
#define MSS_EXAMPLE_MISMATCH

#include <fstream>
#include <sstream>
#include "../src/core/Solution.h"
#include "../src/post/Output.h"
#include "../src/post/check/Continuity.h"

using namespace mss;

class Mismatch {
 public:
  Mismatch(double omega, const Eigen::VectorXcd& w, const Eigen::VectorXcd& t,
           const Material& m0, double width = 0.6, double height = 0.6,
           double density = 500)
      : omega_(omega),
        width_(width),
        height_(height),
        density_(density),
        w_(w),
        t_(t),
        m0_(m0) {}

  double operator()(const Eigen::Vector4d& r) const {
    Matrix matrix(m0_ * r, omega_);
    Boundary<AP, 4> b{density_, {{0, height_}, {width_, 0}}, &matrix};
    return t_.dot(b.MatrixH() * w_ - b.MatrixG() * t_).real();
  }

  Material material(const Eigen::Vector4d& r) const { return m0_ * r; }

 private:
  double omega_, width_, height_, density_;
  Eigen::VectorXcd w_, t_;
  const Material m0_;
};
void read_bv(const std::string& fn, Eigen::VectorXcd& w,
             Eigen::VectorXcd& t) {
  std::ifstream file(fn);
  std::string tmp;
  double x, y, ang;
  dcomp t1, t2;
  std::vector<dcomp> ww, tt;
  while (getline(file, tmp)) {
    std::stringstream(tmp) >> x >> y >> ang >> t1 >> t2;
    ww.push_back(t1);
    tt.push_back(t2);
  }
  w = Eigen::Map<Eigen::VectorXcd>(&ww[0], ww.size());
  t = Eigen::Map<Eigen::VectorXcd>(&tt[0], tt.size());
}
void compute_bv(double omega, Eigen::VectorXcd& w, Eigen::VectorXcd& t,
                double x1, double y1, double x2, double y2,
                double density = 500) {
  input::Solution in("input.txt");
  in.update_frequency(omega);
  Solution<AP> s(in);
  s.Solve();

  post::CC_Solution<AP> cc{&s};
  std::cout << mss_msg({"Maximum mismatch: ", std::to_string(cc.Max())})
            << std::endl;

  Boundary<AP, 4> b{density, {{x1, y1}, {x2, y2}}, s.Matrix()};

  w.resize(b.NumNode());
  t.resize(b.NumNode());
#ifdef NDEBUG
#pragma omp parallel for
#endif
  for (size_t i = 0; i < b.NumNode(); i++) {
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

  std::cout << mss_msg({std::to_string(b.NumNode()),
                        " pairs of boundary values computed."})
            << std::endl;
}

#endif
