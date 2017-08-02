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

#ifndef MSS_CONFIGASSEMBLY_H
#define MSS_CONFIGASSEMBLY_H

#include "../pre/Input.h"
#include "ConfigFiber.h"
#include "Fiber.h"
#include "Inhomogeneity.h"

namespace mss {

template <typename T>
class ConfigAssembly;

template <typename T>
using ConfigAssemPtrs = std::vector<ConfigAssembly<T>*>;

template <typename T>
using ConfigAssemCPtrs = std::vector<const ConfigAssembly<T>*>;

template <typename T>
class ConfigAssembly {
 public:
  ConfigAssembly(const std::string& ID, const input::ConfigAssembly& input,
                 const Matrix* matrix)
      : ID_(ID), matrix_(matrix), input_(input) {
    add_inhomo();
    allocate();
  }

  virtual ~ConfigAssembly() { delete_inhomo(); }

  const Eigen::MatrixXcd& TransMatrix() const;

  const double& CharLength() const { return height_ + width_; }
  size_t NoN() const;  // TODO
  size_t NoE() const;  // TODO
  size_t NoC() const;  // TODO

  const double& Height() const { return height_; }
  const double& Width() const { return width_; }
  const std::string& ID() const { return ID_; }

  void Solve(const InciCPtrs<T>& incident);

  Inhomogeneity<T>* InWhich(const CS* objCS) const;
  T Resultant(const CS* objCS, const Inhomogeneity<T>* inhomo,
              const InciCPtrs<T>& incident) const;
  T Resultant(const CS* objCS, const InciCPtrs<T>& incident) const;

  void PrintCoeff(std::ostream& os) const;

  const InhomoCPtrs<T>& Inhomo() const { return inhomoC_; }
  const Inhomogeneity<T>* Inhomo(const size_t& sn) const {
    return inhomoC_[sn];
  }

 protected:
  const std::string ID_;
  InhomoPtrs<T> inhomo_;
  InhomoCPtrs<T> inhomoC_;
  ConfigFiberCPtrs<T> configFiber_;
  ConfigAssemCPtrs<T> configAssembly_;

  const size_t P_      = {0};                // TODO
  const double height_ = {0}, width_ = {0};  // TODO

  const class Matrix* matrix_;
  const input::ConfigAssembly& input_;

  Eigen::MatrixXcd Q_;
  Eigen::MatrixXcd C_;

  void add_inhomo();
  void add_fiber();
  void add_configFiber();
  void delete_inhomo();
  void delete_configFiber();
  void allocate();
  void compute_MatrixC();
  Eigen::VectorXcd inVect(const InciCPtrs<T>& incident);
  void distSolution(const Eigen::VectorXcd& solution);
};

// ---------------------------------------------------------------------------
// Inline functions:

template <typename T>
void ConfigAssembly<T>::Solve(const InciCPtrs<T>& incident) {
  // Jacobi SVD:
  compute_MatrixC();
  auto svd = C_.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::VectorXcd solution = svd.solve(inVect(incident));
  distSolution(solution);
}

template <typename T>
Eigen::VectorXcd ConfigAssembly<T>::inVect(const InciCPtrs<T>& incident) {
  // The effect vector of incident wave along all the interfaces inside the
  // assembly.

  size_t n = 0, u = 0;
  for (auto& i : inhomo_) n += i->NoE();
  Eigen::VectorXcd rst(n);
  for (auto& i : inhomo_) {
    rst.segment(u, i->NoE()) = i->InciVect(incident);
    u += i->NoE();
  }
  return rst;
}

template <typename T>
void ConfigAssembly<T>::distSolution(const Eigen::VectorXcd& solution) {
  size_t u = 0;
  for (auto& i : inhomo_) {
    i->SetCoeff(solution.segment(u, i->NoC()));
    u += i->NoC();
  }
}

template <typename T>
Inhomogeneity<T>* ConfigAssembly<T>::InWhich(const CS* objCS) const {
  // Return the pointer to the inhomogeneity in which the objCS is.
  // The local CS is considered as the global CS.

  Inhomogeneity<T>* rst = nullptr;
  for (auto& i : inhomo_)
    if (i->Contain(objCS)) rst = i;
  return rst;
}

template <typename T>
T ConfigAssembly<T>::Resultant(const CS* objCS, const Inhomogeneity<T>* in,
                               const InciCPtrs<T>& incident) const {
  T rst(objCS);
  if (in)
    rst = in->Inner(objCS);
  else {
    for (auto& i : incident) rst += i->Effect(objCS);
    for (auto& i : inhomo_) rst += i->Scatter(objCS);
  }
  return rst;
}

template <typename T>
T ConfigAssembly<T>::Resultant(const CS* objCS,
                               const InciCPtrs<T>& incident) const {
  return Resultant(objCS, InWhich(objCS), incident);
}

template <typename T>
void ConfigAssembly<T>::PrintCoeff(std::ostream& os) const {
  for (auto& i : inhomo_) i->PrintCoeff(os);
}

template <typename T>
void ConfigAssembly<T>::add_inhomo() {
  add_fiber();
  for (auto& i : inhomo_) inhomoC_.push_back(i);
}

template <typename T>
void ConfigAssembly<T>::add_fiber() {
  add_configFiber();
  for (auto& i : input_.fiber)
    inhomo_.push_back(
        new Fiber<T>(FindPtrID(configFiber_, i.configID), i.position));
}

template <typename T>
void ConfigAssembly<T>::add_configFiber() {
  for (auto& i : *input_.configFiber)
    configFiber_.push_back(new ConfigFiber<T>(i, matrix_));
}

template <typename T>
void ConfigAssembly<T>::delete_inhomo() {
  for (auto& i : inhomo_) delete i;
  delete_configFiber();
}

template <typename T>
void ConfigAssembly<T>::delete_configFiber() {
  for (auto& i : configFiber_) delete i;
}

template <typename T>
void ConfigAssembly<T>::allocate() {
  size_t noe = 0, noc = 0;
  for (auto& i : inhomo_) {
    noe += i->NoE();
    noc += i->NoC();
  }
  C_.resize(noe, noc);  // The combined transform matrix.
}

template <typename T>
void ConfigAssembly<T>::compute_MatrixC() {
  for (size_t v = 0; v < inhomo_.size(); v++) {
    int i = 0, j = 0;
    for (size_t k = 0; k < v; k++) j += inhomo_[k]->NoC();
    int Nv = inhomo_[v]->NoC();
    for (size_t u = 0; u < inhomo_.size(); u++) {
      int Nu                 = inhomo_[u]->NoE();
      C_.block(i, j, Nu, Nv) = inhomo_[v]->ModeMatrix(inhomo_[u]);
      i += Nu;
    }
  }
}

}  // namespace mss

#endif
