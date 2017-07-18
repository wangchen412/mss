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

#include "ConfigAssembly.h"

namespace mss {

template <typename T>
void ConfigAssembly<T>::Solve(const std::vector<Incident<T>*>& incident) {
  // Jacobi SVD:
  auto svd = C_.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::VectorXcd solution = svd.solve(inVect(incident));
  distSolution(solution);
}

template <typename T>
Eigen::VectorXcd ConfigAssembly<T>::inVect(
    const std::vector<Incident<T>*>& incident) {
  // The effect vector of incident wave along all the interfaces inside the
  // assembly.

  size_t n = 0, u = 0;
  for (auto& i : inhomo_) n += i->NoE();
  Eigen::VectorXcd rst(n);
  for (auto& i : inhomo_) {
    for (auto& j : incident) rst.segment(u, i->NoE()) += j->Effect(i->Node());
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
T ConfigAssembly<T>::Resultant(
    const CS* objCS, const Inhomogeneity<T>* in,
    const std::vector<Incident<T>*>& incident) const {
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
void ConfigAssembly<T>::add_inhomo() {
  add_fiber();
}
template <typename T>
void ConfigAssembly<T>::add_fiber() {
  add_configFiber();
  for (auto& i : input_.fiber)
    inhomo_.pushback(
        new Fiber<T>(FindID(configFiber_, i.configID), i.position));
}
template <typename T>
void ConfigAssembly<T>::add_configFiber() {
  for (auto& i : *input_.configFiber)
    configFiber_.pushback(new ConfigFiber<T>(i, this->Matrix()));
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
  for (size_t u = 0; u < inhomo_.size(); u++) {
    int i = 0, j = 0;
    for (size_t k = 0; k < u; k++) i += inhomo_[k]->NoE();
    int Nu = inhomo_[u]->NoE();
    for (size_t v = 0; v < inhomo_.size(); v++) {
      size_t Nv = inhomo_[v]->NoE();
      C_.block(i, j, Nu, Nv) = inhomo_[u]->ModeMatrix(inhomo_[v]);
      j += Nv;
    }
  }
}

}  // namespace mss
