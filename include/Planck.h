// Copyright (c) 2000-2017, Texas Engineering Experiment Station (TEES), a
// component of the Texas A&M University System.
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are not permitted without specific prior written permission
// from TEES.

// If written permission is obtained for redistribution or further use, the
// following conditions must be met:

// 1) Redistributions of source code must retain the above copyright notice,
// this list of conditions and the disclaimer below.

// 2) Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions, and the disclaimer below in the documentation and/or
// other materials provided with the distribution.

// 3) Neither the name of TEES, the name of the Texas A&M University System, nor
// the names of its contributors may be used to endorse or promote products
// derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS AS IS
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#pragma once

#include <vector>
#include <cmath>
#include <cassert>
#include <limits>

#include "NDArray.h"
#include "Energy.h"
#include "Constants.h"

/*!
  \file Planck.h
  \brief Planck integration using:
  Infinite summation of Taylor expansion,
  Gaussian quadrature.

  All Planck integration functions should obtain B and dBdT simultaneously
  for all energy groups in the energy group set.

  Energy difference from infinite integral is put into the LAST energy
  group to conserve energy.
*/
class Planck
{
private:
  /// Accuracy parameter
  double m_accuracy;

  /// Points and weights are for integration by gaussian quadrature
  std::vector<long double> m_points;
  std::vector<long double> m_weights;

  /// Sets guassian quadrature for integration
  void gauss_quad_setup();

  bool equal(const double l, const double r, const int ulp = 2)
  {
    // The machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::fabs(l - r) <= std::numeric_limits<double>::epsilon() * std::fabs(l + r) * ulp ||
      std::fabs(l - r) < std::numeric_limits<double>::min(); // unless the result is subnormal
  }

public:

  Planck(double accuracy = std::numeric_limits<double>::epsilon());
  ~Planck() {}

  /// @brief returns the planck function B
  double get_B(double T, double E)
  {
    assert(T >= 0.0);
    assert(E >= 0.0);

    if(equal(T, 0.0))
      return 0.0;

    //                     2 E^3
    //    B(E, T) = ----------------------
    //              h^3 c^2 (e^(E/kT) - 1)

    const double h = Constants::PLANCK_CONSTANT;
    const double k = Constants::BOLTZMANN_CONSTANT;
    const double c = Constants::SPEED_OF_LIGHT;

    return 2.0 * std::pow(E, 3.0) * std::pow(h, -3.0) * std::pow(c, -2.0) /
      (exp(E/(k*T)) - 1.0);
  }

  /// @brief returns the partial of B with respect to temperature
  double get_dBdT(double T, double E)
  {
    assert(T >= 0.0);
    assert(E >= 0.0);

    if(equal(T, 0.0))
      return 0.0;

    //    dB(E, T)       2      E^4      e^(E/kT)
    //    -------- = ---------  ---  ----------------
    //       dT      h^3 c^2 k  T^2  (e^(E/kT) - 1)^2

    const double h = Constants::PLANCK_CONSTANT;
    const double k = Constants::BOLTZMANN_CONSTANT;
    const double c = Constants::SPEED_OF_LIGHT;

    return 2.0 * std::pow(h,-3.0) * std::pow(c,-2.0) * std::pow(k,-1.0) *
      std::pow(E, 4.0) * std::pow(T,-2.0) * exp(E/(k*T)) * std::pow(exp(E/(k*T)) - 1.0, -2.0);
}

  /// @brief returns Bg, integrated in energy from E_min to E_max
  double integrate_B(double T, double E_min, double E_max);

  /// @brief returns grey-case B
  double integrate_B_grey(double T);

  /// @brief returns (dB/dT)g, integrated in energy from E_min to E_max
  double integrate_dBdT(double T, double E_min, double E_max);

  /// @brief grey-case (dB/dT)g
  double integrate_dBdT_grey(double T);

  // This is a bit of a hack, but not too bad
  void get_Planck(double Te, std::vector<Group> &edisc,
		  ndarray::array<double, 1> &B, ndarray::array<double, 1> &dBdT);

  void get_Planck(double T, std::vector<Group> &edisc,
		  std::vector<double> &B, std::vector<double> &dBdT);
};
