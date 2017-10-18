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

#ifndef MSS_AREA_H
#define MSS_AREA_H

#include "PointSet.h"

namespace mss {

namespace post {

template <typename T>
class Area : public PointSet<T> {
  using PointSet<T>::point_;

 public:
  Area(const Solution<T>* solution, const PosiVect& p1, const PosiVect& p2,
       size_t Nx, size_t Ny, const std::string& id = "1")
      : PointSet<T>(solution, "Area_" + id),
        p1_(p1),
        p2_(p2),
        Nx_(Nx),
        Ny_(Ny) {
    // Add points:
    PosiVect dx((p2 - p1).x / Nx, 0);
    PosiVect dy(0, (p2 - p1).y / Ny);
    point_.resize(Nx * Ny);

#ifdef NDEBUG
#pragma omp parallel for
#endif
    for (size_t j = 0; j < Ny; j++)
      for (size_t i = 0; i < Nx; i++)
        point_[j * Nx + i] = new Point<T>(solution, p1 + dx * i + dy * j);
  }

  std::string Shape() const override { return "Area"; }
  std::ostream& PrintParam(std::ostream& os) const override {
    return os << p1_ << "\t\t" << p2_ << "\t\t" << Nx_ << "\t\t" << Ny_
              << std::endl;
  }

 private:
  const PosiVect p1_, p2_;
  const double Nx_, Ny_;
};

typedef Area<IP> AreaIP;
typedef Area<AP> AreaAP;

}  // namespace post

}  // namespace mss

#endif
