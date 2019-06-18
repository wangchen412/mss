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

class res {
 public:
  res(IncidentPlaneSH* in, Boundary<AP, 14>* b0, Boundary<AP, 14>* b1)
      : in(in), b0(b0), b1(b1) {}

  StateAP Resultant(const CS* cs) const {
    switch (b1->Contains(cs)) {
      case 1:
        return b1->Effect(cs);
      case -1:
        return in->Effect(cs) + b0->Effect(cs);
    }
    return b1->Effect(cs, true);
  }

 private:
  IncidentPlaneSH* in;
  Boundary<AP, 14>* b0;
  Boundary<AP, 14>* b1;
};

class Mismatch {
 public:
  Mismatch(double omega, const Eigen::VectorXcd& w, const Eigen::VectorXcd& t,
           const Material& m0)
      : omega_(omega), w_(w), t_(t), m0_(m0) {}
  double operator()(const Eigen::Vector4d& r) const {
    Matrix matrix(m0_.mul_comp(r), omega_);
    Boundary<AP, 4> b{0, {}, &matrix, INPUT};
    return (b.MatrixH() * w_ - b.MatrixG() * t_).norm();
  }

  Material material(const Eigen::Vector4d& r) const {
    return m0_.mul_comp(r);
  }

 private:
  double omega_;
  Eigen::VectorXcd w_, t_;
  const Material m0_;
};

int main() {
  // 1. Solve and output the wave field of multiple scattering.
  Solution<AP> s{input::Solution("input.txt")};
  s.Solve();
  post::CC_Solution<AP> cc{&s};
  std::cout << mss_msg({"Maximum mismatch: ", std::to_string(cc.Max())})
            << std::endl;
  // post::Output<AP> o{&s};
  // o.Write();

  std::ofstream coeff_file("coeff.txt");
  for (auto& i : s.inhomo()) coeff_file << i->ScatterCoeff();
  coeff_file.close();

  // 2. Compute boundary values around the RVE.
  //   Boundary<AP, 4> b(0, {}, s.Matrix(), INPUT);
  //   Eigen::VectorXcd w(b.NumNode()), t(b.NumNode());
  //   std::vector<StateAP> v(b.Node().size());
  // #ifdef NDEBUG
  // #pragma omp parallel for
  // #endif
  //   for (size_t i = 0; i < b.NumNode(); i++) {
  //     v[i] = s.Resultant(b.Node(i));
  //     w(i) = v[i].Bv()(0);
  //     t(i) = v[i].Bv()(1);
  //   }
  //   std::ofstream bv_out("bv.dat");
  //   for (size_t i = 0; i < b.NumNode(); i++)
  //     bv_out << setMaxPrecision << v[i].Basis()->PositionGLB() << "\t"
  //            << v[i].Basis()->AngleGLB() << "\t" << w(i) << "\t" << t(i)
  //            << std::endl;
  //   bv_out.close();

  Boundary<AP, 4> box(500, {{-0.3, 0.3}, {0.3, -0.3}}, s.Matrix());
  Eigen::VectorXcd w(box.NumNode()), t(box.NumNode());
  std::vector<StateAP> v(box.Node().size());

  std::cout << box.Node().size() << std::endl;

#ifdef NDEBUG
#pragma omp parallel for
#endif
  for (size_t k = 0; k < box.NumNode(); k++) {
    v[k] = s.Resultant(box.Node(k));
    w(k) = v[k].Bv()(0);
    t(k) = v[k].Bv()(1);
  }

  // 3. Homogenization.
  Mismatch f(s.Frequency(), w, t, {{11400, 11400}, 0, {84e9, 84e9}});
  std::ofstream file("iterations.dat");
  Eigen::VectorXd p = BasinHopping(10, 10, f, Eigen::Vector4d::Ones(), &file);
  file.close();
  std::cout << mss_msg({"Effective properties: "}) << p.transpose()
            << std::endl;

  // 4. Solve and output the wave field of effective homogeneous scatterer.
  Material steel(7670, 116e9, 84.3e9);
  mss::Matrix m{steel, s.Frequency()}, ff{f.material(p), s.Frequency()};

  IncidentPlaneSH* in = new IncidentPlaneSH(m);
  Boundary<AP, 14>* b0 = new Boundary<AP, 14>(50 * m.KT(), {{-2, 2}, {2, -2}},
                                              &m, RECTANGULAR, true);
  Boundary<AP, 14>* b1 =
      new Boundary<AP, 14>(50 * m.KT(), {{-2, 2}, {2, -2}}, &ff);

  Eigen::VectorXcd bv =
      b1->DispToEffect() * (b0->MatrixH() + b0->MatrixG() * b1->DtN())
                               .lu()
                               .solve(in->EffectDv(b0->Node()));
  b1->SetBv(bv);
  for (long i = 0; i < bv.size(); i++)
    if (i % 2) bv(i) *= -1;
  b0->SetBv(bv);

  // res* rst = new res(in, b0, b1);
  // post::Area<AP> a2(rst, {-6, 6}, {6, -6}, 1200, 1200, "bem");
  // a2.Write();

  delete in;
  delete b0;
  delete b1;
  return 0;
}
