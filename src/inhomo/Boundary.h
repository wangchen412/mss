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

#ifndef MSS_BOUNDARY_H
#define MSS_BOUNDARY_H

#include "Panel.h"

namespace mss {

template <typename T, int N = 4>
class Boundary {
 public:
  Boundary(double density, const std::vector<PosiVect>& positions,
           const Matrix* matrix, Tessellation type = RECTANGULAR)
      : density_(density), matrix_(matrix) {
    switch (type) {
      case RECTANGULAR:
        assert(positions.size() == 2);
        add_rect(positions[0], positions[1]);
        break;
      default:
        exit_error_msg({"Wrong boundary type."});
    }
    P_ = node_.size();
  }
  virtual ~Boundary() {
    for (auto& i : node_) delete i;
    for (auto& i : panel_) delete i;
  }

  const CSCPtrs& Node() const { return node_; }
  MatrixXcd EffectMatT(const CS* objCS) const;
  MatrixXcd EffectMatT(const CSCPtrs& objCSs) const;
  VectorXcd EffectBvT(const Inhomo<T>* obj, const VectorXcd& psi) const;

 private:
  double density_;
  CSCPtrs node_;
  PanelCPtrs<T, N> panel_;
  size_t P_;
  const Matrix* matrix_;
  int n_{T::NumBv};

  void add_rect(const PosiVect& p1, const PosiVect& p2);
  void add_line(const PosiVect& p1, const PosiVect& p2);
};

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T, int N>
MatrixXcd Boundary<T, N>::EffectMatT(const mss::CS* objCS) const {
  MatrixXcd rst(n_, n_ * P_);
  for (size_t i = 0; i < P_; i++)
    rst.block(0, n_ * i, n_, n_) = panel_[i]->InfMatT(objCS);
  return rst;
}

template <typename T, int N>
MatrixXcd Boundary<T, N>::EffectMatT(const CSCPtrs& objCSs) const {
  MatrixXcd rst(n_ * objCSs.size(), n_ * P_);

#ifdef NDEBUG
#pragma omp parallel for
#endif
  for (size_t i = 0; i < objCSs.size(); i++)
    rst.block(n_ * i, 0, n_, n_ * P_) = EffectMatT(objCSs[i]);
  return rst;
}

template <typename T, int N>
VectorXcd Boundary<T, N>::EffectBvT(const Inhomo<T>* obj,
                                    const VectorXcd& psi) const {
  return EffectMatT(obj->Node()) * psi;
}

template <typename T, int N>
void Boundary<T, N>::add_rect(const PosiVect& p1, const PosiVect& p2) {
  add_line({p1.x, p1.y}, {p1.x, p2.y});
  add_line({p1.x, p2.y}, {p2.x, p2.y});
  add_line({p2.x, p2.y}, {p2.x, p1.y});
  add_line({p2.x, p1.y}, {p1.x, p1.y});
}

template <typename T, int N>
void Boundary<T, N>::add_line(const PosiVect& p1, const PosiVect& p2) {
  size_t n   = (p2 - p1).Length() * density_;
  PosiVect d = (p2 - p1) / n;
  double len = d.Length();
  double ang = d.Angle() - pi_2;
  for (size_t i = 0; i < n; i++) {
    node_.push_back(new CS(p1 + d * (i + 0.5), ang));
    panel_.push_back(new Panel<T, N>(node_.back(), len, matrix_));
  }
}

}  // namespace mss

#endif
