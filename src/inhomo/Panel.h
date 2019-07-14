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

#ifndef MSS_PANEL_H
#define MSS_PANEL_H

#include "../core/Modes.h"
#include "Inhomo.h"

namespace mss {

template <typename T, int N = 2>
class Panel {
 public:
  Panel(const CS* center, double length, const Matrix* matrix)
      : center_(center), hl_(length / 2.0), matrix_(matrix) {
    add_gauss_point();
  }
  virtual ~Panel() {
    for (auto& i : point_) delete i;
  }

  Eigen::Matrix<dcomp, T::NumV, T::NumBv> InfStateMatT(const CS* objCS) const;
  MatrixNcd<T> InfMatT(const CS* objCS) const;
  MatrixNcd<T> InfMat_L(const CS* objCS) const;
  MatrixNcd<T> InfMat_T(const CS* objCS) const;

 private:
  // The center CS of the panel, with the out-going normal
  // direction along the x-positive.
  const CS* center_;
  CSCPtrs point_;  // The Gaussian Integration points.
  double hl_;      // Half length of the panel.
  const Matrix* matrix_;
  static LegendreRoot<N> l_;

  void add_gauss_point();
};

template <typename T, int N>
LegendreRoot<N> Panel<T, N>::l_;

template <typename T, int N = 2>
using PanelCPtrs = std::vector<const Panel<T, N>*>;

template <typename T, int N = 2>
using PanelPtrs = std::vector<Panel<T, N>*>;

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T, int N>
Eigen::Matrix<dcomp, T::NumV, T::NumBv> Panel<T, N>::InfStateMatT(
    const CS* objCS) const {
  Eigen::Matrix<dcomp, T::NumV, T::NumBv> rst;
  rst.setZero();
  for (int i = 0; i < N; i++)
    rst += GreenStateT<T>(point_[i], objCS, matrix_) * l_.weight(i);
  return rst * hl_;
}

// TODO Need renaming. InfBvMatT. (Not all the state components.)
template <typename T, int N>
MatrixNcd<T> Panel<T, N>::InfMatT(const CS* objCS) const {
  MatrixNcd<T> rst;
  rst.setZero();
  for (int i = 0; i < N; i++)
    rst += GreenT<T>(point_[i], objCS, matrix_) * l_.weight(i);
  return rst * hl_;
}

// TEMP In-plane problem test.
template <typename T, int N>
MatrixNcd<T> Panel<T, N>::InfMat_L(const CS* objCS) const {
  MatrixNcd<T> rst;
  rst.setZero();
  for (int i = 0; i < N; i++)
    rst += Green(point_[i], objCS, matrix_->KL_comp()) * l_.weight(i);
  return rst * hl_;
}
template <typename T, int N>
MatrixNcd<T> Panel<T, N>::InfMat_T(const CS* objCS) const {
  MatrixNcd<T> rst;
  rst.setZero();
  for (int i = 0; i < N; i++)
    rst += Green(point_[i], objCS, matrix_->KT_comp()) * l_.weight(i);
  return rst * hl_;
}

template <typename T, int N>
void Panel<T, N>::add_gauss_point() {
  for (int i = 0; i < N; i++) {
    CS* pp = new CS(*center_);
    pp->Translate(0, l_.root(i) * hl_);
    point_.push_back(pp);
  }
}

}  // namespace mss
#endif
