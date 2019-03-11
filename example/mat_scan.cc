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

#include <fstream>
#include <sstream>
#include "../src/core/Solution.h"
#include "../src/post/Output.h"
#include "../src/post/check/Continuity.h"

using namespace mss;

void ReadBv(Eigen::VectorXcd& w, Eigen::VectorXcd& t) {
  std::ifstream file("bv.dat");
  std::string tmp;
  double x, y, ang;
  for (int i = 0; i < 1200; i++) {
    getline(file, tmp);
    std::stringstream ss(tmp);
    ss >> x >> y >> ang >> w(i) >> t(i);
  }
}

double BoundaryMismatch(const Eigen::VectorXcd& w, const Eigen::VectorXcd& t,
                        const Material& material) {
  Matrix m(material, 16576.243191120248 * 1.95);
  Boundary<AP, 4> b{500, {{-0.3, 0.3}, {0.3, -0.3}}, &m};

  return (b.MatrixH() * w - b.MatrixG() * t).norm();
}

void Scan(const Material& m1, const Material& m2, size_t n) {
  double dr_rho = (m2.Rho_comp() - m1.Rho_comp()).real() / n;
  double di_rho = (m2.Rho_comp() - m1.Rho_comp()).imag() / n;

  double dr_mu = (m2.Mu_comp() - m1.Mu_comp()).real() / n;
  double di_mu = (m2.Mu_comp() - m1.Mu_comp()).imag() / n;

  Eigen::VectorXcd w(1200), t(1200);
  ReadBv(w, t);

  double rst[n][n][n][n];

#ifdef NDEBUG
#pragma omp parallel for
#endif
  for (size_t jj = 0; jj < n * n; jj++) {
    size_t j = jj % n;
    size_t i = (jj - j) / n;
    for (size_t k = 0; k < n; k++)
      for (size_t l = 0; l < n; l++) {
        Material m({m1.Rho_comp().real() + i * dr_rho,
                    m1.Rho_comp().imag() + j * di_rho},
                   0,
                   {m1.Mu_comp().real() + k * dr_mu,
                    m1.Mu_comp().imag() + l * di_mu});
        rst[i][j][k][l] = BoundaryMismatch(w, t, m);
      }
  }

  std::ofstream file("mismatch.dat");
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      for (size_t k = 0; k < n; k++)
        for (size_t l = 0; l < n; l++)
          file << setMaxPrecision << m1.Rho_comp().real() + i * dr_rho << "\t"
               << m1.Rho_comp().imag() + j * di_rho << "\t"
               << m1.Mu_comp().real() + k * dr_mu << "\t"
               << m1.Mu_comp().imag() + l * di_mu << "\t" << rst[i][j][k][l]
               << std::endl;

  file.close();
}

int main(int argc, char* argv[]) {
  if (argc != 9) exit_error_msg({"Material parameters required."});
  dcomp rho1{atof(argv[1]), atof(argv[2])}, mu1{atof(argv[3]), atof(argv[4])};
  dcomp rho2{atof(argv[5]), atof(argv[6])}, mu2{atof(argv[7]), atof(argv[8])};

  Scan({rho1, mu1, mu1}, {rho2, mu2, mu2}, 16);

  return 0;
}
