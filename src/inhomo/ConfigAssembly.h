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

#include "../core/Input.h"
#include "ConfigFiber.h"
#include "Fiber.h"
#include "Inhomogeneity.h"

namespace mss {

template <typename T>
class ConfigAssembly {
 public:
  ConfigAssembly(const input::ConfigAssembly& input, const Matrix* matrix)
      : ID_(input.ID), matrix_(matrix), input_(input) {
    add_inhomo(input);
    allocate();
  }

  virtual ~ConfigAssembly() { delete_inhomo(); }

  const Eigen::MatrixXcd& TransMatrix() const;

  const double& CharLength() const { return height_ + width_; }
  size_t NoN() const;
  size_t NoE() const;
  size_t NoC() const;

  const double& Height() const { return height_; }
  const double& Width() const { return width_; }

  void Solve(const std::vector<Incident<T>*>& incident);

  Inhomogeneity<T>* InWhich(const CS* objCS) const;
  T Resultant(const CS* objCS, const Inhomogeneity<T>* inhomo,
              const std::vector<Incident<T>*>& incident) const;

 protected:
  std::vector<Inhomogeneity<T>*> inhomo_;
  std::vector<ConfigFiber<T>*> configFiber_;
  std::vector<ConfigAssembly<T>*> configAssembly_;
  const std::string ID_;
  const size_t P_;
  const double height_, width_;
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
  Eigen::VectorXcd inVect(const std::vector<Incident<T>*>& incident);
  void distSolution(const Eigen::VectorXcd& solution);
};

}  // namespace mss

#endif
