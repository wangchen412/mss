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

template <typename T, int N = 2>
class Boundary {
 public:
  Boundary(double density, const std::vector<PosiVect>& positions,
           const Matrix* matrix, BoundaryShape shape = RECTANGULAR)
      : density_(density), matrix_(matrix) {
    switch (shape) {
      case RECTANGULAR:
        assert(positions.size() == 2);
        add_rect(positions[0], positions[1]);
        break;
      default:
        exit_error_msg({"Wrong boundary type."});
    }
    P_ = node_.size();
  }
  Boundary(double density, double height, double width, const Matrix* matrix)
      : Boundary(density, {{0, height}, {width, 0}}, matrix) {}
  virtual ~Boundary() {
    for (auto& i : node_) delete i;
    for (auto& i : node_c_) delete i;
    for (auto& i : panel_) delete i;
  }

  const CSCPtrs& Node() const { return node_; }
  const CSCPtrs& DNode();
  MatrixXcd EffectMatT(const CS* objCS) const;
  MatrixXcd EffectMatT(const CSCPtrs& objCSs) const;
  MatrixXcd EffectMatT(const InhomoCPtrs<T>& objs) const;
  VectorXcd EffectBvT(const Inhomo<T>* obj, const VectorXcd& psi) const;

 private:
  double density_;
  CSCPtrs node_;
  CSCPtrs node_c_;  // Complementary nodes.
  CSCPtrs node_d_;  // Doubled nodes.
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
const CSCPtrs& Boundary<T, N>::DNode() {
  if (!node_d_.empty()) return node_d_;

  for (size_t i = 0; i < node_.size() - 1; i++) {
    PosiVect p = (node_[i]->Position() + node_[i + 1]->Position()) / 2;
    double ang = (node_[i]->Angle() + node_[i + 1]->Angle()) / 2;
    node_c_.push_back(new CS(p, ang));
  }
  PosiVect p = (node_.back()->Position() + node_.begin()->Position()) / 2;
  double ang = (node_.back()->Angle() + node_.begin()->Angle()) / 2;
  node_c_.push_back(new CS(p, ang));

  for (size_t i = 0; i < node_.size(); i++) {
    node_d_.push_back(node_[i]);
    node_d_.push_back(node_c_[i]);
  }

  return node_d_;
}

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
MatrixXcd Boundary<T, N>::EffectMatT(const InhomoCPtrs<T>& objs) const {
  size_t m = 0;
  for (auto& i : objs) m += i->NumBv();

  MatrixXcd rst(m, n_ * P_);
  for (size_t u = 0; u < objs.size(); u++) {
    m = 0;
    for (size_t k = 0; k < u; k++) m += objs[k]->NumBv();
    rst.block(m, 0, objs[u]->NumBv(), n_ * P_) = EffectMatT(objs[u]->Node());
  }

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
