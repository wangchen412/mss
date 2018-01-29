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

// Virtual base class of inhomogeneity classes.
// Contains local CS and information about the matrix.

#ifndef MSS_INHOMO_H
#define MSS_INHOMO_H

#include "../incident/Incident.h"

namespace mss {

enum InhomoType { FIBER, ASSEMBLY };

template <typename T>
class Inhomo {
 public:
  explicit Inhomo(const PosiVect& position, InhomoType type, double angle = 0,
                  const CS* basis = nullptr)
      : localCS_(position, angle, basis), type_(type) {}

  virtual ~Inhomo() {}

  // Collocation matrix.
  virtual const MatrixXcd& ColloMat() const {
    print_error_msg({"The Collocation matrix is not available."});
    exit(EXIT_FAILURE);
  }

  // Check if the position of the objective CS is inside the inhomogeneity.
  virtual const Inhomo* Contains(const CS* objCS) const = 0;

  // Transform matrix.
  virtual const MatrixXcd& TransMat() const = 0;

  // The matrix consists the modes from this inhomo to the nodes of the
  // objective inhomo.
  MatrixXcd ModeMat(const Inhomo* obj) const;

  // Resultant states.
  virtual T Scatter(const CS* objCS) const = 0;
  virtual T Inner(const CS* objCS) const = 0;

  // The nth Modes. The n should be the serial number, instead of the order.
  // The transformation from the serial number to the order should be done in
  // the derived class.
  virtual T ScatterMode(const CS* objCS, size_t sn) const = 0;
  VectorXcd ScatterBv(const CSCPtrs& objCSs, size_t sn) const;

  virtual size_t NumNode() const = 0;
  virtual size_t NumCoeff() const = 0;
  virtual size_t NumBv() const = 0;

  VectorXcd IncVec(const InciCPtrs<T>& inc) const;
  VectorXcd TransIncVec(const VectorXcd& incBv) { return TransMat() * incBv; }
  VectorXcd TransIncVec(const InciCPtrs<T>& inc) {
    return TransIncVec(IncVec(inc));
  }

  // Pass the scattering coefficents and set the inner coefficients with
  // T-matrix.
  virtual void SetCoeff(const VectorXcd&) = 0;

  virtual const VectorXcd& ScatterCoeff() const = 0;
  virtual void PrintCoeff(std::ostream& os) const;

  // This inhomogeneity's nodes.
  virtual const CSCPtrs& Node() const = 0;
  virtual const CS* Node(size_t i) const { return Node()[i]; }

  const CS* LocalCS() const { return &localCS_; }
  const CS* Basis() const { return localCS_.Basis(); }
  const PosiVect& Position() const { return localCS_.Position(); }
  PosiVect PositionGLB() const { return localCS_.PositionGLB(); }
  double Angle() const { return localCS_.Angle(); }
  double AngleGLB() const { return localCS_.AngleGLB(); }
  const InhomoType& Type() const { return type_; }

 protected:
  const CS localCS_;
  InhomoType type_;
};

template <typename T>
using InhomoPtrs = std::vector<Inhomo<T>*>;

template <typename T>
using InhomoCPtrs = std::vector<const Inhomo<T>*>;

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T>
VectorXcd Inhomo<T>::IncVec(const InciCPtrs<T>& inc) const {
  VectorXcd rst(NumBv());
  rst.setZero();
  for (auto& i : inc) rst += i->EffectBv(Node());
  return rst;
}
template <typename T>
VectorXcd Inhomo<T>::ScatterBv(const CSCPtrs& objCSs, size_t sn) const {
  VectorXcd rst(objCSs.size() * T::NumBv);
  for (size_t i = 0; i < objCSs.size(); i++)
    rst.segment(i * T::NumBv, T::NumBv) = ScatterMode(objCSs[i], sn).Bv();
  return rst;
}
template <typename T>
MatrixXcd Inhomo<T>::ModeMat(const Inhomo<T>* obj) const {
  MatrixXcd rst(obj->NumBv(), NumCoeff());
  for (size_t sn = 0; sn < NumCoeff(); sn++)
    rst.col(sn) = ScatterBv(obj->Node(), sn);
  return rst;
}
template <typename T>
void Inhomo<T>::PrintCoeff(std::ostream& os) const {
  os << setMaxPrecision << ScatterCoeff() << std::endl;
}

}  // namespace mss

#endif
