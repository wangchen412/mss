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

template <typename T, int N = 4>
class Panel {
 public:
  Panel(const CS* center, double length, const Matrix* matrix)
      : center_(center), hl_(length / 2.0), matrix_(matrix) {
    add_gauss_point();
  }
  virtual ~Panel() {
    for (auto& i : point_) delete i;
  }

  auto InfMatT(const CS* objCS) const;

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

template <typename T, int N = 4>
using PanelCPtrs = std::vector<const Panel<T, N>*>;

template <typename T, int N = 4>
using PanelPtrs = std::vector<Panel<T, N>*>;

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T, int N>
auto Panel<T, N>::InfMatT(const CS* objCS) const {
  Eigen::Matrix<dcomp, T::NumBv, T::NumBv> rst;
  rst.setZero();
  for (int i = 0; i < N; i++)
    rst += GreenT<T>(point_[i], objCS, matrix_) * l_.root(i);
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
