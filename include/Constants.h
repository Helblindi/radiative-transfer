#pragma once

#include <cmath>

namespace Constants
{

static const double PLANCK_CONSTANT_JS     = 6.626083e-35;         // jk-sh
static const double PLANCK_CONSTANT        = 4.141895e-10;         // keV-sh
static const double BOLTZMANN_CONSTANT     = 1.0;                  // keV/keV
static const double BOLTZMANN_CONSTANT_JPK = 1.601558e-25;         // jk/keV
static const double SPEED_OF_LIGHT         = 299.79245800;         // cm/sh
static const double PI                     = 3.1415926546;
static const double FOUR_PI                = 4.0*PI;
static const double RADIATION_CONSTANT_A   = 1.3653104e-2;         // jk/(cm^3-keV^4)
static const double KELVIN2KEV             = 8.6173281e-8;         // keV/K
static const double NATURAL_LOG_2          = 0.6931471806;         // natural log of 2

using std::pow;
static const double RADIATION_CONSTANT_A_LONG = (8.0*pow(PI,5)*pow(BOLTZMANN_CONSTANT,4))/
  (15.0*pow(PLANCK_CONSTANT,3)*pow(SPEED_OF_LIGHT,3));   // keV/cm^3-KeV^4

}
