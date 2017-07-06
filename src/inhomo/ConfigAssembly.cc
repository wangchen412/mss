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
  for (auto& i : input_.configFiber)
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
  int noe = 0, noc = 0;
  for (auto& i : inhomo_) {
    noe += i->NoE();
    noc += i->NoC();
  }
  C_.resize(noe, noc);  // The combined transform matrix.
}

}  // namespace mss
