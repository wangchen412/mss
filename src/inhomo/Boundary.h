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

#include <fstream>
#include <sstream>
#include "Panel.h"

namespace mss {

template <typename T, int N = 2>
class Boundary {
 public:
  Boundary(double density, const std::vector<PosiVect>& positions,
           const Matrix* matrix, BoundaryShape shape = RECTANGULAR,
           bool external = false)
      : shape_(shape),
        positions_(positions),
        density_(density),
        matrix_(matrix),
        ext_(external) {
    switch (shape) {
      case INCL_RECT:
        for (size_t i = 0; i < positions.size() / 4; i++)
          add_quad(std::vector<PosiVect>(positions.begin() + i * 4,
                                         positions.begin() + i * 4 + 4));
        r_cc_ = (positions[0] - positions[2]).Length() / 2;
        center_ = CS((positions[0] + positions[2]) / 2);

        height_ = (positions_[0] - positions_[1]).Length();
        width_ = (positions_[2] - positions_[1]).Length();
        angle_ = (positions_[2] - positions_[1]).Angle(width_);
        break;
      case RECTANGULAR:
        for (size_t i = 0; i < positions.size() / 2; i++)
          add_rect(positions[i * 2], positions[i * 2 + 1]);
        r_cc_ = (positions[0] - positions[1]).Length() / 2;
        center_ = CS((positions[0] + positions[1]) / 2);
        break;
      case CIRCULAR:
        assert(positions.size() == 2);
        r_cc_ = (positions[0] - positions[1]).Length();
        center_ = CS(positions[0]);
        add_circle(positions[0], r_cc_);
        break;
      case INPUT:
        add_input_node();
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
    for (auto& i : panel_) delete i;
  }

  const CSCPtrs& Node() const { return node_; }
  const CS* Node(size_t i) const { return node_[i]; }
  const std::vector<CSCPtrs>& Edge() const { return edge_; }
  const CSCPtrs& Edge(size_t i) const { return edge_[i]; }
  size_t NumNode() const { return P_; }
  size_t NumBv() const { return NumNode() * T::NumBv; }
  size_t NumDv() const { return NumNode() * T::NumDv; }
  size_t NumCoeff() const { return 2 * N_ + 1; }
  MatrixXcd EffectStateMatT(const CS* objCS) const;
  MatrixXcd EffectMatT(const CS* objCS) const;
  MatrixXcd EffectMatT(const CSCPtrs& objCSs) const;
  MatrixXcd EffectMatT(const InhomoCPtrs<T>& objs) const;
  VectorXcd EffectBvT(const Inhomo<T>* obj, const VectorXcd& psi) const;

  void ReverseEdge();

  // Boundary element method.
  // Influence matrices.
  MatrixXcd MatrixH();
  MatrixXcd MatrixG();
  MatrixXcd DtN();

  // TEMP Should be in in-plane problem.
  MatrixXcd MatrixHL();
  MatrixXcd MatrixGL();
  MatrixXcd MatrixHT();
  MatrixXcd MatrixGT();
  MatrixXcd MatrixH_IP();
  MatrixXcd MatrixG_IP();

  template <int NN>
  MatrixXcd CombinedMatrix(Boundary<T, NN>& other);

  // Transformation from displacement to displacement and traction.
  MatrixXcd DispToEffect();

  // Representation integral with displacement only.
  MatrixXcd DispMatT(const CSCPtrs& objCSs);

  // These four methods are for the tests which are about expanding the wave
  // field inside the boundary with cylindrical wave modes. For the circular
  // boundary, it works well. But for the rectangular one, the collocation
  // matrix has large condition number, because the cylindrical modes are not
  // orthogonal along the rectangular boundary.
  MatrixXcd ColloMatT();
  MatrixXcd ModeMatT(const CS* objCS) const;
  MatrixXcd ModeMatT(const CSCPtrs& objCSs) const;
  MatrixXcd ModeMatT(const InhomoCPtrs<T>& objs) const;

  const VectorXcd& Bv() const { return bv_; }
  void SetBv(const VectorXcd& bv) { bv_ = bv; }
  T Effect(const CS* objCS, bool overlap = false) const;
  int Contains(const CS* objCS) const;

 private:
  BoundaryShape shape_;
  std::vector<PosiVect> positions_;
  double density_;
  CSCPtrs node_;
  std::vector<CSCPtrs> edge_;
  PanelCPtrs<T, N> panel_;
  size_t P_;
  const Matrix* matrix_;
  const int n_{T::NumBv};
  MatrixXcd c_;
  bool c_computed_{false};
  int N_{20};    // TODO The top order of the incident wave expansion. TEMP
  int pr_{200};  // The ratio between the point densities.
  double r_cc_;  // Radius of the circumscribed circle.
  CS center_;    // Center of the circumscribed circle.

  MatrixXcd H_, G_, DtN_;  // Influence matrices.
  MatrixXcd HL_, GL_, HT_, GT_;
  bool HG_computed_{false}, DtN_computed_{false};
  bool HGL_computed_{false}, HGT_computed_{false};

  bool ext_;
  VectorXcd bv_;
  double epsilon_{1e-4};

  // Counterclockwise
  void add_quad(const std::vector<PosiVect>& ps);
  // top-left -> bottom-right
  void add_rect(const PosiVect& p1, const PosiVect& p2);
  void add_line(const PosiVect& p1, const PosiVect& p2);
  void add_circle(const PosiVect& p, double r);
  void add_input_node();
  void compute_HG();
  void compute_HG_L();
  void compute_HG_T();

  double height_;
  double width_;
  double angle_;
};

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T, int N>
MatrixXcd Boundary<T, N>::MatrixH() {
  if (HG_computed_) return H_;
  compute_HG();
  return H_;
}
template <typename T, int N>
MatrixXcd Boundary<T, N>::MatrixG() {
  if (HG_computed_) return G_;
  compute_HG();
  return G_;
}
template <typename T, int N>
MatrixXcd Boundary<T, N>::MatrixHL() {
  if (HGL_computed_) return HL_;
  compute_HG_L();
  return HL_;
}
template <typename T, int N>
MatrixXcd Boundary<T, N>::MatrixGL() {
  if (HGL_computed_) return GL_;
  compute_HG_L();
  return GL_;
}
template <typename T, int N>
MatrixXcd Boundary<T, N>::MatrixHT() {
  if (HGT_computed_) return HT_;
  compute_HG_T();
  return HT_;
}
template <typename T, int N>
MatrixXcd Boundary<T, N>::MatrixGT() {
  if (HGT_computed_) return GT_;
  compute_HG_T();
  return GT_;
}
template <typename T, int N>
MatrixXcd Boundary<T, N>::MatrixH_IP() {
  MatrixXcd rst(P_ * 2, P_ * 2);
  MatrixXcd zero(P_, P_);
  zero.setZero();
  rst << MatrixHL(), zero, MatrixHT(), zero;
  return rst;
}
template <typename T, int N>
MatrixXcd Boundary<T, N>::MatrixG_IP() {
  MatrixXcd rst(P_ * 2, P_ * 2);
  MatrixXcd zero(P_, P_);
  zero.setZero();
  rst << MatrixGL(), zero, MatrixGT(), zero;
  return rst;
}

template <typename T, int N>
MatrixXcd Boundary<T, N>::DtN() {
  if (DtN_computed_) return DtN_;
  DtN_.resize(P_, P_);
  DtN_ = MatrixG().inverse() * MatrixH();
  DtN_computed_ = true;
  return DtN_;
}
template <typename T, int N>
MatrixXcd Boundary<T, N>::DispToEffect() {
  MatrixXcd rst(NumBv(), NumDv());
  MatrixXcd I = MatrixXcd::Identity(NumDv(), NumDv());
  for (size_t i = 0; i < NumDv(); i++) {
    rst.row(i * 2) = I.row(i);
    rst.row(i * 2 + 1) = DtN().row(i);
  }
  return rst;
}
template <typename T, int N>
void Boundary<T, N>::compute_HG() {
  H_.resize(P_, P_);
  G_.resize(P_, P_);

#ifdef NDEBUG
#pragma omp parallel for
#endif
  for (size_t i = 0; i < P_; i++)
    for (size_t j = 0; j < P_; j++) {
      MatrixNcd<T> m = panel_[j]->InfMatT(node_[i]);
      H_(i, j) = -m(0, 0);
      G_(i, j) = m(0, 1);
      if (i == j) H_(i, j) += 0.5;
    }
  HG_computed_ = true;
}
template <typename T, int N>
void Boundary<T, N>::compute_HG_L() {
  HL_.resize(P_, P_);
  GL_.resize(P_, P_);

#ifdef NDEBUG
#pragma omp parallel for
#endif
  for (size_t i = 0; i < P_; i++)
    for (size_t j = 0; j < P_; j++) {
      MatrixNcd<T> m = panel_[j]->InfMat_L(node_[i]);
      HL_(i, j) = -m(0, 0);
      GL_(i, j) = m(0, 1);
      if (i == j) HL_(i, j) += 0.5;
    }
  HGL_computed_ = true;
}
template <typename T, int N>
void Boundary<T, N>::compute_HG_T() {
  HT_.resize(P_, P_);
  GT_.resize(P_, P_);

#ifdef NDEBUG
#pragma omp parallel for
#endif
  for (size_t i = 0; i < P_; i++)
    for (size_t j = 0; j < P_; j++) {
      MatrixNcd<T> m = panel_[j]->InfMat_T(node_[i]);
      HT_(i, j) = -m(0, 0);
      GT_(i, j) = m(0, 1);
      if (i == j) HT_(i, j) += 0.5;
    }
  HGT_computed_ = true;
}
template <typename T, int N>
template <int NN>
MatrixXcd Boundary<T, N>::CombinedMatrix(Boundary<T, NN>& other) {
  long n = this->NumDv();
  MatrixXcd rst(n * 2, n * 2);
  rst.block(0, 0, n, n) = this->MatrixH();
  rst.block(0, n, n, n) = this->MatrixG();
  rst.block(n, 0, n, n) = other.MatrixH();
  rst.block(n, n, n, n) = -other.MatrixG();
  return rst;
}
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
MatrixXcd Boundary<T, N>::EffectStateMatT(const CS* objCS) const {
  const size_t m_ = T::NumV;
  MatrixXcd rst(m_, n_ * P_);
  for (size_t i = 0; i < P_; i++)
    rst.block(0, n_ * i, m_, n_) = panel_[i]->InfStateMatT(objCS);
  return rst;
}

// TODO Need renaming. EffectBvMatT. (Not all the state components.)
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
MatrixXcd Boundary<T, N>::DispMatT(const CSCPtrs& objCSs) {
  // The displacement transfering matrix
  // Not optimal. TODO: derive the displacement representation directly
  // without computing the traction.

  MatrixXcd tmp = EffectMatT(objCSs) * DispToEffect();
  MatrixXcd rst(tmp.rows() / 2, tmp.cols());
  for (long i = 0; i < rst.rows(); i++) rst.row(i) = tmp.row(i * 2);
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
T Boundary<T, N>::Effect(const CS* objCS, bool overlap) const {
  if (overlap) {
    PosiVect p = center_.PositionIn(objCS);
    CS tmpCS = *objCS + p / p.Length() * epsilon_;
    return Eigen::Matrix<dcomp, T::NumV, 1>(EffectStateMatT(&tmpCS) * bv_);
  }
  return T(Eigen::Matrix<dcomp, T::NumV, 1>(EffectStateMatT(objCS) * bv_));
}
template <typename T, int N>
int Boundary<T, N>::Contains(const CS* objCS) const {
  PosiVect p = objCS->PositionGLB();
  switch (shape_) {
    case INCL_RECT: {
      CS local(positions_[1], angle_);
      PosiVect p = objCS->PositionIn(&local);
      if (p.x > epsilon_ && p.x < width_ - epsilon_ && p.y > epsilon_ &&
          p.y < height_ - epsilon_)
        return 1;
      else if (p.x < -epsilon_ || p.x > width_ + epsilon_ ||
               p.y > height_ + epsilon_ || p.y < -epsilon_)
        return -1;
    }
    case RECTANGULAR:
      if (p.x > positions_[0].x + epsilon_ &&
          p.x < positions_[1].x - epsilon_ &&
          p.y < positions_[0].y - epsilon_ &&
          p.y > positions_[1].y + epsilon_)
        return 1;
      else if (p.x < positions_[0].x - epsilon_ ||
               p.x > positions_[1].x + epsilon_ ||
               p.y > positions_[0].y + epsilon_ ||
               p.y < positions_[1].y - epsilon_)
        return -1;
    case CIRCULAR:
      double r = (center_.PositionGLB() - p).Length();
      if (r < r_cc_ - epsilon)
        return 1;
      else if (r > r_cc_ + epsilon)
        return -1;
  }
  return 0;
}
template <typename T, int N>
void Boundary<T, N>::ReverseEdge() {
  for (size_t i = edge_.size() - 1; i > edge_.size() / 2 - 1; i--)
    std::reverse(edge_[i].begin(), edge_[i].end());
}

template <typename T, int N>
void Boundary<T, N>::add_quad(const std::vector<PosiVect>& ps) {
  add_line(ps[0], ps[1]);
  add_line(ps[1], ps[2]);
  add_line(ps[2], ps[3]);
  add_line(ps[3], ps[0]);
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
  size_t n = (p2 - p1).Length() * density_;
  size_t ne = n / pr_;
  PosiVect d = (p2 - p1) / n;
  PosiVect de = (p2 - p1) / ne;
  double len = d.Length();
  double ang = d.Angle(len) - pi_2 + ext_ * pi;
  edge_.push_back(CSCPtrs());

  for (size_t i = 0; i < n; i++) {
    node_.push_back(new CS(p1 + d * (i + 0.5), ang));
    panel_.push_back(new Panel<T, N>(node_.back(), len, matrix_));
  }
  for (size_t i = 0; i < ne; i++)
    edge_.back().push_back(new CS(p1 + de * (i + 0.5), ang));
}
template <typename T, int N>
void Boundary<T, N>::add_circle(const PosiVect& p, double r) {
  size_t n = pi2 * r * density_;
  double t = pi2 / n;
  for (size_t i = 0; i < n; i++) {
    node_.push_back(
        new CS(p + PosiVect(r, t * i).Cartesian(), t * i + ext_ * pi));
    panel_.push_back(
        new Panel<T, N>(node_.back(), 2 * tan(pi / n) * r, matrix_));
  }
}
template <typename T, int N>
void Boundary<T, N>::add_input_node() {
  std::ifstream file("boundary.txt");
  std::string tmp;
  PosiVect p1, p2;
  getline(file, tmp);
  std::stringstream(tmp) >> p1;
  while (getline(file, tmp)) {
    std::stringstream(tmp) >> p2;
    PosiVect d = p2 - p1;
    double len = d.Length();
    double ang = d.Angle(len) - pi_2 + ext_ * pi;
    node_.push_back(new CS(p1 + d * 0.5, ang));
    panel_.push_back(new Panel<T, N>(node_.back(), len, matrix_));
    p1 = p2;
  }
}

}  // namespace mss

#endif
