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
#include "Configuration.h"
#include "Fiber.h"
#include "Inhomogeneity.h"

namespace mss {

template <typename T>
class ConfigAssembly : public Configuration<T> {
 public:
  ConfigAssembly(const input::ConfigAssembly& input, const Matrix* matrix)
      : Configuration<T>(matrix), ID_(input.ID), input_(input) {
    add_inhomo(input);
    allocate();
  }

  virtual ~ConfigAssembly() {
    delete_inhomo();
  }

  const Eigen::MatrixXcd& TransMatrix() const override;

  const double& CharLength() const override { return height_ + width_; }
  int NoP() const override;
  int NoE() const override;
  int NoC() const override;

  const double& Height() const { return height_; }
  const double& Width() const { return width_; }

 protected:
  const input::ConfigAssembly& input_;
  std::vector<Inhomogeneity<T>*> inhomo_;
  std::vector<ConfigFiber<T>*> configFiber_;
  const std::string ID_;
  const int N_;
  const int P_;
  const double height_, width_;
  Eigen::MatrixXcd Q_;
  Eigen::MatrixXcd C_;

  void add_inhomo();
  void add_fiber();
  void add_configFiber();
  void delete_inhomo();
  void delete_configFiber();
  void allocate();
};

}  // namespace mss

#endif
