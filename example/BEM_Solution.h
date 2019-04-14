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

#ifndef MSS_BEM_SOLUTION_H
#define MSS_BEM_SOLUTION_H

#include <fstream>
#include <sstream>
#include "../src/core/Solution.h"
#include "../src/post/Output.h"

using namespace mss;

template <typename T, int N>
class BEM_Solution {
 public:
  BEM_Solution(double omega, const Material& inhomo_mat,
               const Material& matrix_mat = {7670, 116e9, 84.3e9},
               const std::vector<PosiVect>& positions = {{-2, 2}, {2, -2}})
      : omega_(omega),
        matrix_(matrix_mat, omega),
        inhomo_(inhomo_mat, omega) {
    in_ = new IncidentPlaneSH(matrix_);
    b0_ = new Boundary<T, N>(10 * matrix_.KT(), positions, &matrix_,
                             RECTANGULAR, true);
    b1_ = new Boundary<T, N>(10 * matrix_.KT(), positions, &inhomo_,
                             RECTANGULAR);

    Eigen::VectorXcd bv =
        b1_->DispToEffect() * (b0_->MatrixH() + b0_->MatrixG() * b1_->DtN())
                                  .lu()
                                  .solve(in_->EffectDv(b0_->Node()));

    b1_->SetBv(bv);
    for (long i = 0; i < bv.size(); i++)
      if (i % 2) bv(i) *= -1;
    b0_->SetBv(bv);
  }

  ~BEM_Solution() {
    delete in_;
    delete b0_;
    delete b1_;
  }

  T Resultant(const CS* cs) const {
    switch (b1_->Contains(cs)) {
      case 1:
        return b1_->Effect(cs);
      case -1:
        return in_->Effect(cs) + b0_->Effect(cs);
    }
    return b1_->Effect(cs, true);
  }

 private:
  double omega_;
  Matrix matrix_, inhomo_;
  Boundary<T, N>* b0_;
  Boundary<T, N>* b1_;
  Incident<T>* in_;
};

#endif
