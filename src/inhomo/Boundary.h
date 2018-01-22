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
        r_cc_   = (positions[0] - positions[1]).Length() / 2;
        center_ = CS((positions[0] + positions[1]) / 2);
        break;
      case CIRCULAR:
        assert(positions.size() == 2);
        r_cc_   = (positions[0] - positions[1]).Length() / 2;
        center_ = CS(positions[0]);
        add_circle(positions[0], r_cc_);
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
  const CS* Node(size_t i) const { return node_[i]; }
  const CSCPtrs& DNode() const { return node_d_; }
  size_t NumDNode() const { return node_d_.size(); }
  size_t NumDBv() const { return NumDNode() * T::NumBv; }
  size_t NumNode() const { return P_; }
  size_t NumNode(int i) const { return nn_[i]; }
  size_t NumBv() const { return NumNode() * T::NumBv; }
  size_t NumCoeff() const { return 2 * N_ + 1; }
  MatrixXcd EffectMatT(const CS* objCS) const;
  MatrixXcd EffectMatT(const CSCPtrs& objCSs) const;
  MatrixXcd EffectMatT(const InhomoCPtrs<T>& objs) const;
  VectorXcd EffectBvT(const Inhomo<T>* obj, const VectorXcd& psi) const;

  // MatrixXcd Extrapolation(const CSCPtrs& inner, size_t P) const;

  // This four methods are for the tests which are about expanding the wave
  // field inside the boundary with cylindrical wave modes. For the circular
  // boundary, it works well. But for the rectangular one, the collocation
  // matrix has large condition number, because the cylindrical modes are not
  // orthogonal along the rectangular boundary.
  MatrixXcd ColloMatT();
  MatrixXcd ModeMatT(const CS* objCS) const;
  MatrixXcd ModeMatT(const CSCPtrs& objCSs) const;
  MatrixXcd ModeMatT(const InhomoCPtrs<T>& objs) const;

 private:
  double density_;
  CSCPtrs node_;
  CSCPtrs node_c_;  // Complementary nodes.
  CSCPtrs node_d_;  // Doubled nodes.
  PanelCPtrs<T, N> panel_;
  size_t P_;
  std::vector<size_t> nn_;  // The number of nodes along each edge.
  const Matrix* matrix_;
  const int n_{T::NumBv};
  MatrixXcd c_;
  bool c_computed_{false};
  int N_{20};    // TODO The top order of the incident wave expansion. TEMP
  double r_cc_;  // Radius of the circumscribed circle.
  CS center_;    // Center of the circumscribed circle.

  void add_rect(const PosiVect& p1, const PosiVect& p2);
  size_t add_line(const PosiVect& p1, const PosiVect& p2);

  void add_circle(const PosiVect& p, double r);
};

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T, int N>
MatrixXcd Boundary<T, N>::ModeMatT(const CS* objCS) const {
  MatrixXcd rst(n_, NumCoeff());
  for (int n = -N_; n <= N_; n++) {
    EigenFunctor J(Jn, n, matrix_->KT(), r_cc_);
    StateAP s = ModeT<T>(&center_, objCS, J, matrix_->Material());
    rst.block<2, 1>(0, n + N_) = s.Bv();
  }
  return rst;
}

template <typename T, int N>
MatrixXcd Boundary<T, N>::ModeMatT(const CSCPtrs& objCSs) const {
  MatrixXcd rst(n_ * objCSs.size(), NumCoeff());

#ifdef NDEBUG
#pragma omp parallel for
#endif
  for (size_t i = 0; i < objCSs.size(); i++)
    rst.block(n_ * i, 0, n_, NumCoeff()) = ModeMatT(objCSs[i]);
  return rst;
}

template <typename T, int N>
MatrixXcd Boundary<T, N>::ModeMatT(const InhomoCPtrs<T>& objs) const {
  size_t m = 0;
  for (auto& i : objs) m += i->NumBv();

  MatrixXcd rst(m, NumCoeff());
  for (size_t u = 0; u < objs.size(); u++) {
    m = 0;
    for (size_t k = 0; k < u; k++) m += objs[k]->NumBv();
    rst.block(m, 0, objs[u]->NumBv(), NumCoeff()) = ModeMatT(objs[u]->Node());
  }

  return rst;
}

template <typename T, int N>
MatrixXcd Boundary<T, N>::ColloMatT() {  // TODO: in-plane
  if (c_computed_) return c_;

  c_.resize(NumBv(), NumCoeff());
  for (int n = -N_; n <= N_; n++) {
    EigenFunctor J(Jn, n, matrix_->KT(), r_cc_);
    for (size_t i = 0; i < P_; i++) {
      StateAP s = ModeT<AP>(&center_, node_[i], J, matrix_->Material());
      c_.block<2, 1>(2 * i, n + N_) = s.Bv();
    }
  }

  c_computed_ = true;

  return c_;
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

// template <typename T, int N>
// MatrixXcd Boundary<T, N>::Extrapolation(const CSCPtrs& inner,
//                                         size_t P) const {
//   // Fit P plane waves at given points. Then extrapolate the field at
//   // boundary pionts.

//   MatrixXcd fit_m;
//   MatrixXcd extra_m;
// }

template <typename T, int N>
void Boundary<T, N>::add_rect(const PosiVect& p1, const PosiVect& p2) {
  nn_.push_back(add_line({p1.x, p1.y}, {p1.x, p2.y}));
  nn_.push_back(add_line({p1.x, p2.y}, {p2.x, p2.y}));
  nn_.push_back(add_line({p2.x, p2.y}, {p2.x, p1.y}));
  nn_.push_back(add_line({p2.x, p1.y}, {p1.x, p1.y}));
}

template <typename T, int N>
size_t Boundary<T, N>::add_line(const PosiVect& p1, const PosiVect& p2) {
  size_t n   = (p2 - p1).Length() * density_;
  PosiVect d = (p2 - p1) / n;
  double len = d.Length();
  double ang = d.Angle() - pi_2;

  // Quadral nodes.
  for (size_t i = 0; i < n; i++) {
    if (i > 0) {
      node_c_.push_back(new CS(p1 + d * (i + 0.25), ang));
      node_d_.push_back(node_c_.back());
    }

    node_.push_back(new CS(p1 + d * (i + 0.5), ang));
    node_d_.push_back(node_.back());

    if (i < n - 1) {
      node_c_.push_back(new CS(p1 + d * (i + 0.75), ang));
      node_d_.push_back(node_c_.back());

      node_c_.push_back(new CS(p1 + d * (i + 1), ang));
      node_d_.push_back(node_c_.back());
    }

    panel_.push_back(new Panel<T, N>(node_.back(), len, matrix_));
  }
  return n;
}

template <typename T, int N>
void Boundary<T, N>::add_circle(const PosiVect& p, double r) {
  size_t n = pi2 * r * density_;
  double t = pi2 / n;
  for (size_t i = 0; i < n; i++)
    node_.push_back(new CS(p + PosiVect(r, t * i).Cartesian(), t * i));
}

}  // namespace mss

#endif
