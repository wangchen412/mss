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

// Incident base class.

#ifndef MSS_INCIDENT_H
#define MSS_INCIDENT_H

#include "../core/Matrix.h"
#include "../core/Modes.h"
#include "../core/State.h"
#include "../pre/Input.h"

namespace mss {

template <typename T>
class Incident {
 public:
  Incident(const Matrix& matrix, double amplitude = 1, double phase = 0)
      : amp_(amplitude),
        phase_(phase),
        material_(matrix.Material()),
        kl_(matrix.KL()),
        kt_(matrix.KT()) {}

  virtual ~Incident() {}

  // The normalization factor.
  T Norm() const;

  // The effect on the origin of the local CS.
  virtual T Effect(const CS* localCS) const = 0;

  // The effect vector along the boundary.
  VectorXcd EffectBv(const CSCPtrs& localCS) const {
    VectorXcd rst(localCS.size() * T::NumBv);
    for (size_t i = 0; i < localCS.size(); i++)
      rst.segment(i * T::NumBv, T::NumBv) = Effect(localCS[i]).Bv();
    return rst;
  }
  // The displacement vector along the boundary.
  VectorXcd EffectDv(const CSCPtrs& localCS) const {
    VectorXcd rst(localCS.size() * T::NumDv);
    for (size_t i = 0; i < localCS.size(); i++)
      rst.segment(i * T::NumDv, T::NumDv) = Effect(localCS[i]).Dv();
    return rst;
  }

  double Amplitude() const { return amp_; }
  double Phase() const { return phase_; }

 protected:
  double amp_, phase_;
  const Material& material_;
  const double kl_, kt_;
};

template <typename T>
using InciPtrs = std::vector<Incident<T>*>;
template <typename T>
using InciCPtrs = std::vector<const Incident<T>*>;

// ---------------------------------------------------------------------------
// Inline functions:

template <>
StateIP Incident<IP>::Norm() const {
  return StateIP(1, 1, material_.C(kl_, kl_, kl_)) * amp_;
}
template <>
StateAP Incident<AP>::Norm() const {
  return StateAP(1, material_.C(kt_, kt_)) * amp_;
}

// ---------------------------------------------------------------------------
// Plane P wave
class IncidentPlaneP : public Incident<IP> {
 public:
  IncidentPlaneP(const Matrix& matrix, double angle = 0, double amplitude = 1,
                 double phase = 0)
      : Incident<IP>(matrix, amplitude, phase),
        a_(angle),
        k_(cos(a_) * matrix.KL(), sin(a_) * matrix.KL()) {}
  IncidentPlaneP(const Matrix& matrix, const input::IncidentPlane& input)
      : IncidentPlaneP(matrix, input.angle, input.amplitude, input.phase) {}

  StateIP Effect(const CS* objCS) const override {
    return PlaneL<IP>(k_, objCS, material_) * amp_ * exp(ii * phase_);
  }
  double Angle() const { return a_; }

 private:
  const double a_;
  const Vector<double> k_;
};

// ---------------------------------------------------------------------------
// Plane SV wave
class IncidentPlaneSV : public Incident<IP> {
 public:
  IncidentPlaneSV(const Matrix& matrix, double angle = 0,
                  double amplitude = 1, double phase = 0)
      : Incident<IP>(matrix, amplitude, phase),
        a_(angle),
        k_(cos(a_) * matrix.KT(), sin(a_) * matrix.KT()) {}

  IncidentPlaneSV(const Matrix& matrix, const input::IncidentPlane& input)
      : IncidentPlaneSV(matrix, input.angle, input.amplitude, input.phase) {}

  StateIP Effect(const CS* objCS) const override {
    return PlaneT<IP>(k_, objCS, material_) * amp_ * exp(ii * phase_);
  }
  double Angle() const { return a_; }

 private:
  const double a_;
  const Vector<double> k_;
};

// ---------------------------------------------------------------------------
// Plane SH wave
class IncidentPlaneSH : public Incident<AP> {
 public:
  IncidentPlaneSH(const Matrix& matrix, double angle = 0,
                  double amplitude = 1, double phase = 0)
      : Incident<AP>(matrix, amplitude, phase),
        a_(angle),
        k_(cos(a_) * matrix.KT(), sin(a_) * matrix.KT()) {}
  IncidentPlaneSH(const Matrix& matrix, const input::IncidentPlane& input)
      : IncidentPlaneSH(matrix, input.angle, input.amplitude, input.phase) {}

  StateAP Effect(const CS* objCS) const override {
    return PlaneT<AP>(k_, objCS, material_) * amp_ * exp(ii * phase_);
  }
  double Angle() const { return a_; }

 private:
  const double a_;
  const Vector<double> k_;
};

class IncidentCylinSH : public Incident<AP> {
 public:
  IncidentCylinSH(const Matrix& matrix, const PosiVect& position,
                  double amplitude = 1, double phase = 0)
      : Incident<AP>(matrix, amplitude, phase), localCS_(position) {}

  StateAP Effect(const CS* objCS) const override {
    const PosiVect pc = objCS->PositionIn(&localCS_);
    const PosiVect p = pc.Polar();
    double r = p.x;
    const CS cs(pc, p.y, &localCS_);

    if (r < 0.1) return StateAP(DispAP(), StressAP(), objCS);

    DispAP w = ii / 4 * Hn(0, kt_ * r);
    dcomp gzr = -ii * kt_ / 4 * Hn(1, kt_ * r);
    dcomp gzt = 0;
    StressAP t = material_.C(gzr, gzt);

    return StateAP(w, t, &cs).in(objCS);
  }

 private:
  const CS localCS_;
};

// ---------------------------------------------------------------------------
// Input generator
template <typename T>
class IncidentGen {
 public:
  IncidentGen(const Matrix& matrix);
  Incident<T>* operator()(const input::IncidentPlane& input) {
    return funcMap[input.type](input);
  }

 private:
  typedef std::function<Incident<T>*(const input::IncidentPlane&)> funcType;
  std::map<std::string, funcType, ci_comp> funcMap;
  const Matrix& matrix_;
};

template <>
IncidentGen<IP>::IncidentGen(const Matrix& matrix) : matrix_(matrix) {
  funcMap["PlaneP"] = [this](const input::IncidentPlane& input) {
    return new IncidentPlaneP(matrix_, input);
  };
  funcMap["PlaneSV"] = [this](const input::IncidentPlane& input) {
    return new IncidentPlaneSV(matrix_, input);
  };
}

template <>
IncidentGen<AP>::IncidentGen(const Matrix& matrix) : matrix_(matrix) {
  funcMap["PlaneSH"] = [this](const input::IncidentPlane& input) {
    return new IncidentPlaneSH(matrix_, input);
  };
}

}  // namespace mss

#endif
